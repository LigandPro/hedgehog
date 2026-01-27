import React from 'react';
import { Box, Text } from 'ink';

interface SearchIndicatorProps {
  active: boolean;
  query: string;
}

export function SearchIndicator({ active, query }: SearchIndicatorProps): React.ReactElement | null {
  if (!active) return null;

  return (
    <Box marginY={1}>
      <Text color="cyan">Search: </Text>
      <Text color="yellow">{query}</Text>
      <Text color="cyan">█</Text>
    </Box>
  );
}

export default SearchIndicator;
