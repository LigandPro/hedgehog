import React, { memo } from 'react';
import { Box, Text } from 'ink';
import type { StageStatus } from '../types/index.js';
import { useTheme } from '../theme/context.js';

// Unicode block characters for smooth progress: ▏▎▍▌▋▊▉█
// Each represents 1/8 increments of a full block
const PROGRESS_BLOCKS = [' ', '▏', '▎', '▍', '▌', '▋', '▊', '▉', '█'];

interface ProgressBarProps {
  progress: number; // 0-100
  width?: number;
  status?: StageStatus;
  label?: string;
  showPercentage?: boolean;
}

export const ProgressBar = memo(function ProgressBar({
  progress,
  width = 24,
  status = 'pending',
  label,
  showPercentage = true,
}: ProgressBarProps): React.ReactElement {
  const { theme } = useTheme();
  // Calculate smooth progress with fractional blocks
  const totalUnits = width * 8; // 8 increments per character
  const filledUnits = Math.round((progress / 100) * totalUnits);

  const fullBlocks = Math.floor(filledUnits / 8);
  const partialBlock = filledUnits % 8;
  const emptyBlocks = width - fullBlocks - (partialBlock > 0 ? 1 : 0);

  const filled = '█'.repeat(fullBlocks);
  const partial = partialBlock > 0 ? PROGRESS_BLOCKS[partialBlock] : '';
  const empty = '░'.repeat(Math.max(0, emptyBlocks));

  const statusIcon = {
    pending: '○',
    running: '●',
    completed: '✓',
    error: '✗',
    skipped: '○',
  }[status];

  const barColor = {
    pending: theme.palette.textMuted,
    running: theme.palette.primary,
    completed: theme.palette.success,
    error: theme.palette.error,
    skipped: theme.palette.textMuted,
  }[status];

  return (
    <Box>
      <Text>{statusIcon} </Text>
      {label && <Text color={theme.palette.text}>{label.padEnd(16)} </Text>}
      <Text color={barColor}>{filled}{partial}</Text>
      <Text color={theme.palette.border}>{empty}</Text>
      {showPercentage && (
        <Text color={theme.palette.textMuted}> {progress.toString().padStart(3)}%</Text>
      )}
    </Box>
  );
});

export default ProgressBar;
