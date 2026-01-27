import React from 'react';
import { Box, Text } from 'ink';
import { useStore } from '../store/index.js';
import type { ToastType } from '../types/index.js';

const TOAST_ICONS: Record<ToastType, string> = {
  success: '✓',
  error: '✗',
  info: 'ℹ',
  warning: '⚠',
};

const TOAST_COLORS: Record<ToastType, string> = {
  success: 'green',
  error: 'red',
  info: 'cyan',
  warning: 'yellow',
};

export function ToastContainer(): React.ReactElement | null {
  const toasts = useStore((state) => state.toasts);

  if (toasts.length === 0) {
    return null;
  }

  return (
    <Box
      flexDirection="column"
      position="absolute"
      marginTop={-toasts.length - 1}
    >
      {toasts.map((toast) => (
        <Box
          key={toast.id}
          borderStyle="round"
          borderColor={TOAST_COLORS[toast.type]}
          paddingX={1}
          marginBottom={0}
        >
          <Text color={TOAST_COLORS[toast.type]} bold>
            {TOAST_ICONS[toast.type]}
          </Text>
          <Text> {toast.message}</Text>
        </Box>
      ))}
    </Box>
  );
}

export default ToastContainer;
