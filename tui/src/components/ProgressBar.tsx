import React, { memo } from 'react';
import { Box, Text } from 'ink';
import type { StageStatus } from '../types/index.js';
import { icons } from '../utils/colors.js';

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
    pending: icons.pending,
    running: icons.running,
    completed: icons.success,
    error: icons.error,
    skipped: icons.pending,
  }[status];

  const barColor = {
    pending: 'gray',
    running: 'cyan',
    completed: 'green',
    error: 'red',
    skipped: 'gray',
  }[status] as 'gray' | 'cyan' | 'green' | 'red';

  return (
    <Box>
      <Text>{statusIcon} </Text>
      {label && <Text>{label.padEnd(16)} </Text>}
      <Text color={barColor}>{filled}{partial}</Text>
      <Text color="gray">{empty}</Text>
      {showPercentage && (
        <Text dimColor> {progress.toString().padStart(3)}%</Text>
      )}
    </Box>
  );
});

export default ProgressBar;
