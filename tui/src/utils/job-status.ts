import type { JobStatus } from '../types/index.js';

export const STATUS_ICONS: Record<JobStatus, string> = {
  completed: '✓',
  running: '●',
  cancelled: '○',
  error: '✗',
};

export const STATUS_COLORS: Record<JobStatus, 'green' | 'yellow' | 'gray' | 'red'> = {
  completed: 'green',
  running: 'yellow',
  cancelled: 'gray',
  error: 'red',
};

export function getStatusIcon(status: JobStatus): string {
  return STATUS_ICONS[status];
}

export function getStatusColor(status: JobStatus): 'green' | 'yellow' | 'gray' | 'red' {
  return STATUS_COLORS[status];
}
