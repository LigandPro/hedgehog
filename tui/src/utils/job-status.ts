import type { JobStatus } from '../types/index.js';

export const STATUS_ICONS: Record<JobStatus, string> = {
  completed: '✓',
  running: '●',
  cancelled: '○',
  error: '✗',
};

export function getStatusIcon(status: JobStatus): string {
  return STATUS_ICONS[status];
}
