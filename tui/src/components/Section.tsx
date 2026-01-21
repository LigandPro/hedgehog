import React from 'react';
import { Box, Text } from 'ink';

interface SectionProps {
  title: string;
  children: React.ReactNode;
  collapsed?: boolean;
}

export function Section({ title, children, collapsed = false }: SectionProps): React.ReactElement {
  return (
    <Box flexDirection="column" marginY={1}>
      <Box marginBottom={1}>
        <Text color="cyan" bold>{title}</Text>
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
