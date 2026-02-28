import React from 'react';
import { Box, Text } from 'ink';
import { useTheme } from '../theme/context.js';

interface SectionProps {
  title: string;
  children: React.ReactNode;
  collapsed?: boolean;
}

export function Section({ title, children, collapsed = false }: SectionProps): React.ReactElement {
  const { theme } = useTheme();
  return (
    <Box flexDirection="column" marginY={1}>
      <Box marginBottom={1}>
        <Text color={theme.palette.accent} bold>{title}</Text>
      </Box>
      {!collapsed && (
        <Box flexDirection="column" paddingLeft={2}>
          {children}
        </Box>
      )}
    </Box>
  );
}

export default Section;
