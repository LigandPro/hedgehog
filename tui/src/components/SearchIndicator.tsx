import React from 'react';
import { Box, Text } from 'ink';
import { useTerminalSize } from '../hooks/useTerminalSize.js';
import { useTheme } from '../theme/context.js';
import { LineFill } from './LineFill.js';

interface SearchIndicatorProps {
  active: boolean;
  query: string;
}

export function SearchIndicator({ active, query }: SearchIndicatorProps): React.ReactElement | null {
  const { theme } = useTheme();
  const { width: terminalWidth } = useTerminalSize();
  if (!active) return null;

  const contentWidth = Math.max(8, terminalWidth - 2);

  return (
    <Box marginY={1} width={contentWidth}>
      <Text color={theme.palette.primary}>Search: </Text>
      <Text color={theme.palette.accent}>{query}</Text>
      <Text color={theme.palette.primary}>█</Text>
      <LineFill width={contentWidth} />
    </Box>
  );
}

export default SearchIndicator;
