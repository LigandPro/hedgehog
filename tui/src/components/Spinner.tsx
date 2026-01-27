import React, { memo } from 'react';
import { Box, Text } from 'ink';
import InkSpinner from 'ink-spinner';

interface SpinnerProps {
  label?: string;
}

export const Spinner = memo(function Spinner({ label }: SpinnerProps): React.ReactElement {
  return (
    <Box>
      <Text color="yellow">
        <InkSpinner type="dots" />
      </Text>
      {label && <Text> {label}</Text>}
    </Box>
  );
});

export default Spinner;
