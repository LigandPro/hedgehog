import { appendFileSync, mkdirSync } from 'fs';
import { dirname } from 'path';
import { formatTimestamp } from './format.js';

const LOG_FILE = process.env.HEDGEHOG_TUI_LOG || '/tmp/hedgehog-tui.log';

export type LogLevel = 'debug' | 'info' | 'warn' | 'error';

function ensureLogDir() {
  try {
    mkdirSync(dirname(LOG_FILE), { recursive: true });
  } catch {
    // Directory might already exist
  }
}

export function log(level: LogLevel, message: string, ...args: unknown[]): void {
  if (process.env.HEDGEHOG_TUI_DEBUG !== '1' && level === 'debug') {
    return;
  }
  
  ensureLogDir();
  const timestamp = formatTimestamp();
  const formattedArgs = args.map(arg => 
    typeof arg === 'object' ? JSON.stringify(arg) : String(arg)
  ).join(' ');
  const logLine = `[${timestamp}] [${level.toUpperCase()}] ${message} ${formattedArgs}\n`;
  
  try {
    appendFileSync(LOG_FILE, logLine);
  } catch {
    // Silently fail if we can't write to log
  }
}

export const logger = {
  debug: (message: string, ...args: unknown[]) => log('debug', message, ...args),
  info: (message: string, ...args: unknown[]) => log('info', message, ...args),
  warn: (message: string, ...args: unknown[]) => log('warn', message, ...args),
  error: (message: string, ...args: unknown[]) => log('error', message, ...args),
};

export default logger;
