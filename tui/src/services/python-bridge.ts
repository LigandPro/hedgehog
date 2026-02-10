import { spawn, type ChildProcess } from 'child_process';
import { RpcClient, type NotificationHandler } from './rpc-client.js';
import { logger } from '../utils/logger.js';
import { fileURLToPath } from 'url';
import { dirname, resolve } from 'path';

// Find project root (where pyproject.toml is)
function findProjectRoot(): string {
  const __filename = fileURLToPath(import.meta.url);
  const __dirname = dirname(__filename);
  // From tui/dist/services/ go up to project root
  return resolve(__dirname, '..', '..', '..');
}

// Notification throttling to reduce UI flickering
class NotificationThrottler {
  private pendingNotifications: Map<string, { method: string; params: Record<string, unknown> }> = new Map();
  private timer: ReturnType<typeof setTimeout> | null = null;
  private handler: NotificationHandler | null = null;
  private throttleMs: number;

  constructor(throttleMs = 100) {
    this.throttleMs = throttleMs;
  }

  setHandler(handler: NotificationHandler | null): void {
    this.handler = handler;
  }

  push(method: string, params: Record<string, unknown>): void {
    // Immediate dispatch for critical notifications
    if (method === 'complete' || method === 'error' || method === 'stage_start' || method === 'stage_complete' || method === 'stage_error') {
      this.flush();
      this.handler?.(method, params);
      return;
    }

    // Throttle progress and log notifications - keep only latest per type
    const key = method === 'progress' ? `progress:${params?.stage}` : method;
    this.pendingNotifications.set(key, { method, params });

    if (!this.timer) {
      this.timer = setTimeout(() => this.flush(), this.throttleMs);
    }
  }

  private flush(): void {
    if (this.timer) {
      clearTimeout(this.timer);
      this.timer = null;
    }

    for (const { method, params } of this.pendingNotifications.values()) {
      this.handler?.(method, params);
    }
    this.pendingNotifications.clear();
  }

  clear(): void {
    if (this.timer) {
      clearTimeout(this.timer);
      this.timer = null;
    }
    this.pendingNotifications.clear();
  }
}

export class PythonBridge {
  private process: ChildProcess | null = null;
  private rpcClient: RpcClient | null = null;
  private isReady = false;
  private readyPromise: Promise<void> | null = null;
  private readyResolve: (() => void) | null = null;
  private readyReject: ((err: Error) => void) | null = null;
  private throttler = new NotificationThrottler(300);

  constructor(
    private command: string = 'uv',
    private args: string[] = ['run', 'python', '-m', 'hedgehog.tui_backend'],
    private cwd: string = findProjectRoot()
  ) {}

  async start(): Promise<void> {
    if (this.process) {
      return;
    }

    logger.info(`Starting Python backend in ${this.cwd}...`);
    logger.info(`Command: ${this.command} ${this.args.join(' ')}`);

    this.readyPromise = new Promise((resolve, reject) => {
      this.readyResolve = resolve;
      this.readyReject = reject;
    });

    this.process = spawn(this.command, this.args, {
      stdio: ['pipe', 'pipe', 'pipe'],
      cwd: this.cwd,
      env: {
        ...process.env,
        PYTHONUNBUFFERED: '1',
      },
    });

    this.rpcClient = new RpcClient((message) => {
      if (this.process?.stdin?.writable) {
        this.process.stdin.write(message);
      }
    });

    this.process.stdout?.on('data', (data: Buffer) => {
      const str = data.toString();
      logger.debug('Python stdout:', str);
      this.rpcClient?.handleData(str);
    });

    this.process.stderr?.on('data', (data: Buffer) => {
      const str = data.toString();
      logger.debug('Python stderr:', str);
      // Check for ready signal
      if (str.includes('HEDGEHOG_TUI_READY')) {
        this.isReady = true;
        this.readyResolve?.();
        this.readyResolve = null;
        this.readyReject = null;
      }
    });

    this.process.on('error', (error) => {
      logger.error('Python process error:', error);
      if (!this.isReady) {
        this.readyReject?.(error instanceof Error ? error : new Error(String(error)));
        this.readyResolve = null;
        this.readyReject = null;
      }
      this.cleanup();
    });

    this.process.on('exit', (code, signal) => {
      logger.info('Python process exited:', code, signal);
      if (!this.isReady) {
        this.readyReject?.(new Error(`Backend exited before ready (code=${code}, signal=${signal})`));
        this.readyResolve = null;
        this.readyReject = null;
      }
      this.cleanup();
    });

    // Wait for backend to be ready (with timeout)
    const timeout = new Promise<void>((_, reject) => {
      setTimeout(() => reject(new Error('Backend startup timeout')), 10000);
    });

    try {
      await Promise.race([this.readyPromise, timeout]);
      logger.info('Python backend ready');
    } catch (error) {
      logger.error('Failed to start Python backend:', error);
      this.stop();
      throw error;
    }
  }

  private cleanup(): void {
    this.throttler.clear();
    this.throttler.setHandler(null);
    this.rpcClient?.cancelPending();
    this.process = null;
    this.rpcClient = null;
    this.isReady = false;
    this.readyPromise = null;
    this.readyResolve = null;
    this.readyReject = null;
  }

  stop(): void {
    if (this.process) {
      this.rpcClient?.cancelPending();
      this.throttler.clear();
      this.process.kill();
      this.cleanup();
    }
  }

  onNotification(handler: NotificationHandler): () => void {
    if (!this.rpcClient) {
      return () => {};
    }
    // Set up throttled handler
    this.throttler.setHandler(handler);
    // Subscribe to raw notifications and push through throttler
    return this.rpcClient.onNotification((method, params) => {
      this.throttler.push(method, params);
    });
  }

  async call<T = unknown>(method: string, params?: Record<string, unknown>): Promise<T> {
    if (!this.rpcClient || !this.isReady) {
      throw new Error('Backend not ready');
    }
    return this.rpcClient.call<T>(method, params);
  }

  // Convenience methods
  async listFiles(path: string, extensions?: string[]): Promise<string[]> {
    return this.call<string[]>('list_files', { path, extensions });
  }

  async loadConfig(configType: string): Promise<Record<string, unknown>> {
    return this.call<Record<string, unknown>>('load_config', { config_type: configType });
  }

  async saveConfig(configType: string, data: Record<string, unknown>): Promise<boolean> {
    return this.call<boolean>('save_config', { config_type: configType, data });
  }

  async validateConfig(configType: string, data: Record<string, unknown>): Promise<{ valid: boolean; errors: string[] }> {
    return this.call('validate_config', { config_type: configType, data });
  }

  async startPipeline(stages: string[]): Promise<string> {
    return this.call<string>('start_pipeline', { stages });
  }

  async getProgress(jobId: string): Promise<Record<string, unknown>> {
    return this.call<Record<string, unknown>>('get_progress', { job_id: jobId });
  }

  async cancelPipeline(jobId: string): Promise<boolean> {
    return this.call<boolean>('cancel_pipeline', { job_id: jobId });
  }

  // History convenience methods
  async getJobHistory(limit = 50): Promise<Record<string, unknown>[]> {
    return this.call<Record<string, unknown>[]>('get_job_history', { limit });
  }

  async addJob(
    jobId: string,
    name: string | null,
    inputPath: string,
    outputPath: string,
    stages: string[]
  ): Promise<Record<string, unknown>> {
    return this.call<Record<string, unknown>>('add_job', {
      job_id: jobId,
      name,
      input_path: inputPath,
      output_path: outputPath,
      stages,
    });
  }

  async updateJob(
    jobId: string,
    status?: string,
    results?: Record<string, unknown>,
    error?: string
  ): Promise<Record<string, unknown> | null> {
    return this.call<Record<string, unknown> | null>('update_job', {
      job_id: jobId,
      status,
      results,
      error,
    });
  }

  async deleteJob(jobId: string): Promise<boolean> {
    return this.call<boolean>('delete_job', { job_id: jobId });
  }
}

// Singleton instance
let bridgeInstance: PythonBridge | null = null;

export function getBridge(): PythonBridge {
  if (!bridgeInstance) {
    const useUv = process.env.HEDGEHOG_PYTHON !== 'python';
    if (useUv) {
      bridgeInstance = new PythonBridge('uv', ['run', 'python', '-m', 'hedgehog.tui_backend']);
    } else {
      bridgeInstance = new PythonBridge('python', ['-m', 'hedgehog.tui_backend']);
    }
  }
  return bridgeInstance;
}

export function destroyBridge(): void {
  if (bridgeInstance) {
    bridgeInstance.stop();
    bridgeInstance = null;
  }
}

export default PythonBridge;
