import React from 'react';
import { Box, Text, useInput } from 'ink';
import { useStore } from '../store/index.js';
import { useTerminalSize } from '../hooks/useTerminalSize.js';
import { useTheme } from '../theme/context.js';

export function ConfirmDialog(): React.ReactElement | null {
  const { width, height } = useTerminalSize();
  const { theme } = useTheme();
  const confirmDialog = useStore((state) => state.confirmDialog);
  const hideConfirm = useStore((state) => state.hideConfirm);

  useInput((input, key) => {
    if (!confirmDialog) return;

    if (input === 'y' || input === 'Y') {
      confirmDialog.onConfirm();
      // Clear dialog after confirm
      useStore.setState({ confirmDialog: null });
    } else if (input === 'n' || input === 'N' || key.escape) {
      hideConfirm();
    }
  }, { isActive: !!confirmDialog });

  if (!confirmDialog) {
    return null;
  }

  const { title, message, confirmLabel = 'Yes', cancelLabel = 'No' } = confirmDialog;

  return (
    <Box
      flexDirection="column"
      position="absolute"
      width={width}
      height={height}
    >
      <Box flexGrow={1} />
      <Box justifyContent="center">
        <Box
          flexDirection="column"
          width={Math.max(1, Math.min(74, width - 4))}
          borderStyle="double"
          borderColor={theme.palette.warning}
          paddingX={2}
          paddingY={1}
        >
          <Box justifyContent="center" marginBottom={1}>
            <Text color={theme.palette.warning} bold>{title}</Text>
          </Box>

          <Box justifyContent="center" marginBottom={1}>
            <Text color={theme.palette.text}>{message}</Text>
          </Box>

          <Box justifyContent="center" gap={4}>
            <Box>
              <Text color={theme.palette.primary} bold>[y]</Text>
              <Text color={theme.palette.text}> {confirmLabel}</Text>
            </Box>
            <Box>
              <Text color={theme.palette.primary} bold>[n]</Text>
              <Text color={theme.palette.text}> {cancelLabel}</Text>
            </Box>
          </Box>
        </Box>
      </Box>
      <Box flexGrow={1} />
    </Box>
  );
}

export default ConfirmDialog;
