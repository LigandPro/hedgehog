import React from 'react';
import { Box, Text } from 'ink';
import type { StageStatus } from '../types/index.js';
import { icons } from '../utils/colors.js';

interface ProgressBarProps {
  progress: number; // 0-100
  width?: number;
  status?: StageStatus;
  label?: string;
  showPercentage?: boolean;
}

export function ProgressBar({
  progress,
  width = 24,
  status = 'pending',
  label,
  showPercentage = true,
}: ProgressBarProps): React.ReactElement {
  const filledWidth = Math.round((progress / 100) * width);
  const emptyWidth = width - filledWidth;
  
  const filled = '█'.repeat(filledWidth);
  const empty = '░'.repeat(emptyWidth);
  
  const statusIcon = {
    pending: icons.pending,
    running: icons.running,
    completed: icons.success,
    error: icons.error,
    skipped: icons.pending,
  }[status];

  const barColor = {
    pending: 'gray',
    running: 'yellow',
    completed: 'green',
    error: 'red',
    skipped: 'gray',
  }[status] as 'gray' | 'yellow' | 'green' | 'red';

  return (
    <Box>
      <Text>{statusIcon} </Text>
      {label && <Text>{label.padEnd(16)} </Text>}
      <Text color={barColor}>{filled}</Text>
      <Text color="gray">{empty}</Text>
      {showPercentage && (
        <Text dimColor> {progress.toString().padStart(3)}%</Text>
      )}
    </Box>
  );
}

export default ProgressBar;
