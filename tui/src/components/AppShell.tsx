import React from 'react';
import { Box } from 'ink';

interface AppShellProps {
  header?: React.ReactNode;
  footer?: React.ReactNode;
  children: React.ReactNode;
  padding?: number;
  paddingX?: number;
  paddingY?: number;
  contentJustify?: 'flex-start' | 'center' | 'flex-end';
}

export function AppShell({
  header,
  footer,
  children,
  padding,
  paddingX,
  paddingY,
  contentJustify = 'flex-start',
}: AppShellProps): React.ReactElement {
  return (
    <Box flexDirection="column" flexGrow={1} padding={padding} paddingX={paddingX} paddingY={paddingY}>
      {header && <Box flexDirection="column">{header}</Box>}
      <Box flexDirection="column" flexGrow={1} justifyContent={contentJustify}>
        {children}
      </Box>
      {footer && <Box flexDirection="column">{footer}</Box>}
    </Box>
  );
}

export default AppShell;
