import React from 'react';
import { Box, Text } from 'ink';
import { useStore } from '../store/index.js';
import { useTheme } from '../theme/context.js';
import { useTerminalSize } from '../hooks/useTerminalSize.js';
import type { ToastType } from '../types/index.js';

const TOAST_ICONS: Record<ToastType, string> = {
  success: '✓',
  error: '✗',
  info: 'ℹ',
  warning: '⚠',
};

function truncate(value: string, maxWidth: number): string {
  if (maxWidth <= 0) return '';
  if (value.length <= maxWidth) return value;
  if (maxWidth <= 1) return '…';
  return `${value.slice(0, maxWidth - 1)}…`;
}

export function ToastContainer(): React.ReactElement | null {
  const { theme } = useTheme();
  const { width } = useTerminalSize();
  const toasts = useStore((state) => state.toasts);

  if (toasts.length === 0) {
    return null;
  }

  const toastColors: Record<ToastType, string> = {
    success: theme.palette.success,
    error: theme.palette.error,
    info: theme.palette.info,
    warning: theme.palette.warning,
  };
  const toast = toasts[toasts.length - 1];
  const messageWidth = Math.max(14, Math.min(58, width - 18));
  const message = truncate(toast.message, messageWidth).padEnd(messageWidth, ' ');

  return (
    <Box
      flexDirection="row"
      position="absolute"
      width={width}
      justifyContent="flex-end"
      marginTop={1}
      paddingRight={2}
    >
      <Box
        key={toast.id}
        borderStyle="single"
        borderColor={toastColors[toast.type]}
        paddingX={1}
      >
        <Text color={toastColors[toast.type]} backgroundColor={theme.palette.panel} bold>
          {TOAST_ICONS[toast.type]}
        </Text>
        <Text color={theme.palette.text} backgroundColor={theme.palette.panel}> {message}</Text>
      </Box>
    </Box>
  );
}

export default ToastContainer;
