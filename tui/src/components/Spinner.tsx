import React from 'react';
import { Box, Text } from 'ink';
import InkSpinner from 'ink-spinner';

interface SpinnerProps {
  label?: string;
}

export function Spinner({ label }: SpinnerProps): React.ReactElement {
  return (
    <Box>
      <Text color="yellow">
        <InkSpinner type="dots" />
      </Text>
      {label && <Text> {label}</Text>}
    </Box>
  );
}

export default Spinner;
