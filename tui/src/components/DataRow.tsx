import React from 'react';
import { Box, Text } from 'ink';
import { useTheme } from '../theme/context.js';

interface DataRowProps {
  label: string;
  value: string | number | boolean;
  labelWidth?: number;
  valueColor?: string;
}

export function DataRow({ 
  label, 
  value, 
  labelWidth = 20,
  valueColor,
}: DataRowProps): React.ReactElement {
  const { theme } = useTheme();
  const displayValue = typeof value === 'boolean' 
    ? (value ? 'Yes' : 'No')
    : String(value);

  return (
    <Box>
      <Text color={theme.palette.textMuted}>{label.padEnd(labelWidth)}</Text>
      <Text color={valueColor || theme.palette.text}>{displayValue}</Text>
    </Box>
  );
}

export default DataRow;
