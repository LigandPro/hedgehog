import React, { memo } from 'react';
import { Box, Text } from 'ink';
import InkSpinner from 'ink-spinner';
import { useTheme } from '../theme/context.js';

interface SpinnerProps {
  label?: string;
}

export const Spinner = memo(function Spinner({ label }: SpinnerProps): React.ReactElement {
  const { theme } = useTheme();
  return (
    <Box>
      <Text color={theme.palette.warning}>
        <InkSpinner type="dots" />
      </Text>
      {label && <Text color={theme.palette.text}> {label}</Text>}
    </Box>
  );
});

export default Spinner;
