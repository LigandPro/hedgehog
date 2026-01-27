import React from 'react';
import { Box, Text, useInput } from 'ink';
import { useStore } from '../store/index.js';

export function ConfirmDialog(): React.ReactElement | null {
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
      borderStyle="double"
      borderColor="yellow"
      paddingX={2}
      paddingY={1}
      marginY={1}
    >
      <Box justifyContent="center" marginBottom={1}>
        <Text color="yellow" bold>{title}</Text>
      </Box>

      <Box justifyContent="center" marginBottom={1}>
        <Text>{message}</Text>
      </Box>

      <Box justifyContent="center" gap={4}>
        <Box>
          <Text color="cyan" bold>[y]</Text>
          <Text> {confirmLabel}</Text>
        </Box>
        <Box>
          <Text color="cyan" bold>[n]</Text>
          <Text> {cancelLabel}</Text>
        </Box>
      </Box>
    </Box>
  );
}

export default ConfirmDialog;
