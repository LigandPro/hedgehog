import React from 'react';
import { Box, Text } from 'ink';

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
  valueColor = 'white',
}: DataRowProps): React.ReactElement {
  const displayValue = typeof value === 'boolean' 
    ? (value ? 'Yes' : 'No')
    : String(value);

  return (
    <Box>
      <Text dimColor>{label.padEnd(labelWidth)}</Text>
      <Text color={valueColor as any}>{displayValue}</Text>
    </Box>
  );
}

export default DataRow;
