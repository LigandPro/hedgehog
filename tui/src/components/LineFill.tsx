import React from 'react';
import { Box, Text } from 'ink';

/**
 * Render trailing spaces to fully overwrite the rest of the current line.
 *
 * Ink generally handles reflows, but when a line becomes shorter between renders
 * (for example, after toggling a view mode or when a value disappears),
 * terminals can show leftover characters at the end of the line.
 *
 * Place this as the LAST child in a horizontal Box with a fixed `width`.
 */
export function LineFill({ width }: { width: number }): React.ReactElement | null {
  if (width <= 0) return null;

  return (
    <Box flexGrow={1}>
      {/* Use plain truncation to avoid rendering an ellipsis character. */}
      <Text wrap="truncate">{' '.repeat(width)}</Text>
    </Box>
  );
}

export default LineFill;
