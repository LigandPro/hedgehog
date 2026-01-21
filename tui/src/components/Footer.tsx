import React from 'react';
import { Box, Text } from 'ink';

interface Shortcut {
  key: string;
  label: string;
}

interface FooterProps {
  shortcuts?: Shortcut[];
}

const defaultShortcuts: Shortcut[] = [
  { key: '↑↓', label: 'Navigate' },
  { key: 'Enter', label: 'Select' },
  { key: 'q', label: 'Quit' },
];

export function Footer({ shortcuts = defaultShortcuts }: FooterProps): React.ReactElement {
  return (
    <Box flexDirection="column" marginTop={1}>
      <Box>
        <Text color="gray">{'─'.repeat(60)}</Text>
      </Box>
      <Box gap={2}>
        {shortcuts.map((shortcut, index) => (
          <Box key={index}>
            <Text color="cyan">[{shortcut.key}]</Text>
            <Text dimColor> {shortcut.label}</Text>
          </Box>
        ))}
      </Box>
    </Box>
  );
}

export default Footer;
