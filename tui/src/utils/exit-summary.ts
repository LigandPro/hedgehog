import type { JobHistoryRecord, Screen } from '../types/index.js';

const EXIT_LOGO = [
  '  █  █ █▀▀ █▀▀▄ █▀▀▀ █▀▀ █  █ █▀▀█ █▀▀▀',
  '  █▀▀█ █▀▀ █  █ █ ▀█ █▀▀ █▀▀█ █  █ █ ▀█',
  '  ▀  ▀ ▀▀▀ ▀▀▀  ▀▀▀▀ ▀▀▀ ▀  ▀ ▀▀▀▀ ▀▀▀▀',
];

interface ExitSummaryInput {
  currentJobId: string | null;
  selectedJobId: string | null;
  jobHistory: JobHistoryRecord[];
  screen: Screen;
}

export function buildExitSummary(input: ExitSummaryInput): string {
  const active = input.currentJobId
    ? input.jobHistory.find((job) => job.id === input.currentJobId)
    : undefined;
  const selected = input.selectedJobId
    ? input.jobHistory.find((job) => job.id === input.selectedJobId)
    : undefined;
  const recent = input.jobHistory[0];
  const sessionJob = active || selected || recent;

  const sessionTitle = sessionJob?.name
    || sessionJob?.id
    || `Last screen: ${input.screen}`;
  const continueCommand = sessionJob?.id
    ? `uv run hedgehog tui -s ${sessionJob.id}`
    : 'uv run hedgehog tui';

  return [
    '',
    ...EXIT_LOGO,
    '',
    `  Session   ${sessionTitle}`,
    `  Continue  ${continueCommand}`,
    '',
  ].join('\n');
}
