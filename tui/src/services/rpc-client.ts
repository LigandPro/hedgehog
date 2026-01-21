import type { RpcRequest, RpcResponse, RpcNotification } from '../types/index.js';
import { logger } from '../utils/logger.js';

export type NotificationHandler = (method: string, params: Record<string, unknown>) => void;

export class RpcClient {
  private requestId = 0;
  private pendingRequests: Map<number, {
    resolve: (value: unknown) => void;
    reject: (error: Error) => void;
  }> = new Map();
  private notificationHandlers: NotificationHandler[] = [];
  private buffer = '';

  constructor(
    private sendMessage: (message: string) => void,
  ) {}

  onNotification(handler: NotificationHandler): () => void {
    this.notificationHandlers.push(handler);
    return () => {
      const index = this.notificationHandlers.indexOf(handler);
      if (index !== -1) {
        this.notificationHandlers.splice(index, 1);
      }
    };
  }

  handleData(data: string): void {
    this.buffer += data;
    
    // Process complete JSON messages
    let newlineIndex: number;
    while ((newlineIndex = this.buffer.indexOf('\n')) !== -1) {
      const line = this.buffer.slice(0, newlineIndex).trim();
      this.buffer = this.buffer.slice(newlineIndex + 1);
      
      if (!line) continue;
      
      try {
        const message = JSON.parse(line);
        this.handleMessage(message);
      } catch (error) {
        logger.error('Failed to parse RPC message:', line, error);
      }
    }
  }

  private handleMessage(message: RpcResponse | RpcNotification): void {
    if ('id' in message && message.id !== undefined) {
      // This is a response
      const pending = this.pendingRequests.get(message.id);
      if (pending) {
        this.pendingRequests.delete(message.id);
        if ('error' in message && message.error) {
          pending.reject(new Error(message.error.message));
        } else {
          pending.resolve(message.result);
        }
      }
    } else if ('method' in message) {
      // This is a notification
      for (const handler of this.notificationHandlers) {
        try {
          handler(message.method, message.params || {});
        } catch (error) {
          logger.error('Notification handler error:', error);
        }
      }
    }
  }

  async call<T = unknown>(method: string, params?: Record<string, unknown>): Promise<T> {
    const id = ++this.requestId;
    const request: RpcRequest = {
      jsonrpc: '2.0',
      id,
      method,
      params,
    };

    return new Promise<T>((resolve, reject) => {
      this.pendingRequests.set(id, {
        resolve: resolve as (value: unknown) => void,
        reject,
      });

      const message = JSON.stringify(request) + '\n';
      logger.debug('RPC request:', message.trim());
      this.sendMessage(message);

      // Timeout after 60 seconds
      setTimeout(() => {
        if (this.pendingRequests.has(id)) {
          this.pendingRequests.delete(id);
          reject(new Error(`RPC call ${method} timed out`));
        }
      }, 60000);
    });
  }

  cancelPending(): void {
    for (const [id, pending] of this.pendingRequests) {
      pending.reject(new Error('Request cancelled'));
    }
    this.pendingRequests.clear();
  }
}

export default RpcClient;
