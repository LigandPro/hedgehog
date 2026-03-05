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
  const { width, height } = useTerminalSize();
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
  const rightPadding = width >= 36 ? 2 : 1;
  const bottomOffset = height >= 24 ? 6 : height >= 16 ? 4 : 2;
  const safeBottomOffset = Math.min(bottomOffset, Math.max(0, Math.floor(height / 3)));
  const messageWidth = Math.max(1, Math.min(58, width - rightPadding - 6));
  const message = truncate(toast.message, messageWidth).padEnd(messageWidth, ' ');

  return (
    <Box
      flexDirection="column"
      position="absolute"
      width={width}
      height={height}
      justifyContent="flex-end"
      alignItems="flex-end"
      paddingRight={rightPadding}
      paddingBottom={safeBottomOffset}
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
