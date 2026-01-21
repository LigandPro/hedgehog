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

export class PythonBridge {
  private process: ChildProcess | null = null;
  private rpcClient: RpcClient | null = null;
  private isReady = false;
  private readyPromise: Promise<void> | null = null;
  private readyResolve: (() => void) | null = null;

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

    this.readyPromise = new Promise((resolve) => {
      this.readyResolve = resolve;
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
      }
    });

    this.process.on('error', (error) => {
      logger.error('Python process error:', error);
      this.cleanup();
    });

    this.process.on('exit', (code, signal) => {
      logger.info('Python process exited:', code, signal);
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
    this.process = null;
    this.rpcClient = null;
    this.isReady = false;
  }

  stop(): void {
    if (this.process) {
      this.rpcClient?.cancelPending();
      this.process.kill();
      this.cleanup();
    }
  }

  onNotification(handler: NotificationHandler): () => void {
    if (!this.rpcClient) {
      return () => {};
    }
    return this.rpcClient.onNotification(handler);
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
