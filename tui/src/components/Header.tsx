import React from 'react';
import { Box, Text } from 'ink';
import { colors, icons } from '../utils/colors.js';

interface HeaderProps {
  title?: string;
  subtitle?: string;
  showLogo?: boolean;
}

const ASCII_LOGO = `    __  __         __           __
   / / / /__  ____/ /___ ____  / /_  ____  ____
  / /_/ / _ \\/ __  / __ \`/ _ \\/ __ \\/ __ \\/ __ \\
 / __  /  __/ /_/ / /_/ /  __/ / / / /_/ / /_/ /
/_/ /_/\\___/\\__,_/\\__, /\\___/_/ /_/\\____/\\__, /
                 /____/                 /____/`;

export function Header({ title, subtitle, showLogo = false }: HeaderProps): React.ReactElement {
  return (
    <Box flexDirection="column" marginBottom={1}>
      {showLogo && (
        <Box flexDirection="column" marginBottom={1}>
          <Text color="cyan">{ASCII_LOGO}</Text>
        </Box>
      )}
      <Box>
        <Text>{icons.hedgehog} </Text>
        <Text color="cyan" bold>HEDGEHOG</Text>
        {title && (
          <>
            <Text color="gray"> │ </Text>
            <Text color="white">{title}</Text>
          </>
        )}
        <Box flexGrow={1} />
        <Text color="gray" dimColor>v1.0.0</Text>
      </Box>
      {subtitle && (
        <Box marginTop={0}>
          <Text dimColor>{subtitle}</Text>
        </Box>
      )}
      <Box marginTop={0}>
        <Text color="gray">{'─'.repeat(60)}</Text>
      </Box>
    </Box>
  );
}

export default Header;
