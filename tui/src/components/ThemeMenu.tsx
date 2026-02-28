import React, { useEffect, useMemo, useState } from 'react';
import { Box, Text, useInput } from 'ink';
import { useStore } from '../store/index.js';
import { useTerminalSize } from '../hooks/useTerminalSize.js';
import { useTheme } from '../theme/context.js';
import { THEMES } from '../theme/themes.js';
import { APP_VERSION } from '../constants/version.js';

const TEMPLATE_MENU_ITEMS = [
  { label: 'New Run', tone: 'primary' as const },
  { label: 'History', tone: 'text' as const },
  { label: 'Config', tone: 'accent' as const },
  { label: 'Reports', tone: 'info' as const },
];

function fitLine(value: string, width: number): string {
  if (width <= 0) return '';
  if (value.length >= width) {
    return value.slice(0, width);
  }
  return value.padEnd(width, ' ');
}

export function ThemeMenu(): React.ReactElement | null {
  const { width, height } = useTerminalSize();
  const { theme, themeId, themeIds, setThemeId, setPreviewThemeId, clearPreviewTheme } = useTheme();
  const showThemeMenu = useStore((state) => state.showThemeMenu);
  const setShowThemeMenu = useStore((state) => state.setShowThemeMenu);

  const [selectedIndex, setSelectedIndex] = useState(0);

  useEffect(() => {
    if (!showThemeMenu) return;
    const index = themeIds.indexOf(themeId);
    setSelectedIndex(index >= 0 ? index : 0);
  }, [showThemeMenu, themeId, themeIds]);

  useEffect(
    () => () => {
      clearPreviewTheme();
    },
    [clearPreviewTheme],
  );

  const close = () => {
    clearPreviewTheme();
    setShowThemeMenu(false);
  };

  const options = useMemo(
    () => themeIds.map((id) => ({ id, label: THEMES[id].name })),
    [themeIds],
  );

  const panelWidth = Math.max(1, Math.min(104, width - 2));
  const contentWidth = Math.max(1, panelWidth - 4);
  const useSidePreview = panelWidth >= 78;
  const previewColumnWidth = useSidePreview
    ? Math.max(30, Math.min(56, Math.floor(panelWidth * 0.5)))
    : contentWidth;
  const listColumnWidth = useSidePreview
    ? panelWidth - previewColumnWidth - 5
    : contentWidth;
  const reservedRows = useSidePreview ? 10 : 18;
  const maxVisibleOptions = Math.max(
    3,
    Math.min(options.length || 1, Math.max(3, height - reservedRows)),
  );
  const startIndex = Math.max(
    0,
    Math.min(
      selectedIndex - Math.floor(maxVisibleOptions / 2),
      Math.max(0, options.length - maxVisibleOptions),
    ),
  );
  const endIndex = Math.min(options.length, startIndex + maxVisibleOptions);
  const visibleOptions = options.slice(startIndex, endIndex);
  const hiddenAbove = startIndex;
  const hiddenBelow = options.length - endIndex;
  const previewOption = options[selectedIndex];
  const previewStatus = previewOption?.id === themeId ? 'Saved' : 'Preview';
  const previewTitle = 'HEDGEHOG';

  useEffect(() => {
    if (!showThemeMenu) return;
    const selected = options[selectedIndex];
    if (!selected) return;
    setPreviewThemeId(selected.id);
  }, [options, selectedIndex, setPreviewThemeId, showThemeMenu]);

  useInput(
    (input, key) => {
      if (!showThemeMenu) return;
      if (options.length === 0) return;

      if (key.upArrow || input === 'k') {
        setSelectedIndex((prev) => (prev <= 0 ? options.length - 1 : prev - 1));
        return;
      }
      if (key.downArrow || input === 'j') {
        setSelectedIndex((prev) => (prev >= options.length - 1 ? 0 : prev + 1));
        return;
      }
      if (key.escape) {
        close();
        return;
      }
      if (key.return) {
        const selected = options[selectedIndex];
        if (selected) {
          setThemeId(selected.id);
        }
        setShowThemeMenu(false);
      }
    },
    { isActive: showThemeMenu },
  );

  if (!showThemeMenu) {
    return null;
  }

  const activeBackground = theme.palette.panelStrong;
  const menuTone = (tone: (typeof TEMPLATE_MENU_ITEMS)[number]['tone']) => {
    switch (tone) {
      case 'primary':
        return theme.palette.primary;
      case 'accent':
        return theme.palette.accent;
      case 'info':
        return theme.palette.info;
      default:
        return theme.palette.text;
    }
  };

  return (
    <Box
      position="absolute"
      width={width}
      height={height}
      flexDirection="column"
    >
      <Box flexGrow={1} />
      <Box justifyContent="flex-end" paddingRight={1}>
        <Box
          width={panelWidth}
          paddingX={2}
          paddingY={1}
          borderStyle="single"
          borderColor={theme.palette.borderActive}
          flexDirection={useSidePreview ? 'row' : 'column'}
        >
          <Box flexDirection="column" width={listColumnWidth}>
            <Box marginBottom={1}>
              <Text color={theme.palette.text} bold>Theme Selector</Text>
            </Box>
            <Text color={theme.palette.textMuted}>Use Up/Down to preview. Enter saves.</Text>
            {hiddenAbove > 0 && (
              <Box marginTop={1} paddingX={1}>
                <Text color={theme.palette.textMuted}>... {hiddenAbove} themes above ...</Text>
              </Box>
            )}
            {visibleOptions.map((option, index) => {
              const absoluteIndex = startIndex + index;
              const active = selectedIndex === absoluteIndex;
              const saved = option.id === themeId;
              const rowBackground = active ? activeBackground : undefined;
              return (
                <Box key={option.id} paddingX={1}>
                  <Text
                    color={active ? theme.palette.text : theme.palette.textMuted}
                    backgroundColor={rowBackground}
                  >
                    {active ? '> ' : '  '}
                  </Text>
                  <Text
                    color={saved ? theme.palette.primary : theme.palette.text}
                    bold={saved}
                    backgroundColor={rowBackground}
                  >
                    {option.label}
                  </Text>
                  <Text color={theme.palette.textMuted} backgroundColor={rowBackground}>
                    {' '}({option.id})
                  </Text>
                  {saved && (
                    <Text
                      color={theme.palette.textMuted}
                      backgroundColor={rowBackground}
                    >
                      {' '}saved
                    </Text>
                  )}
                </Box>
              );
            })}
            {hiddenBelow > 0 && (
              <Box paddingX={1}>
                <Text color={theme.palette.textMuted}>... {hiddenBelow} themes below ...</Text>
              </Box>
            )}
            <Box marginTop={1}>
              <Text color={theme.palette.textMuted}>Enter save  Esc cancel</Text>
            </Box>
          </Box>
          {useSidePreview && (
            <Box width={1} marginX={1}>
              <Text color={theme.palette.border}>|</Text>
            </Box>
          )}
          <Box flexDirection="column" width={previewColumnWidth} marginTop={useSidePreview ? 0 : 1}>
            <Text color="white" bold>
              {fitLine(previewTitle, previewColumnWidth)}
            </Text>
            <Box marginTop={1}>
              <Text color={theme.palette.textMuted}>Version: </Text>
              <Text color={theme.palette.text}>{APP_VERSION}</Text>
            </Box>
            <Box>
              <Text color={theme.palette.textMuted}>Theme: </Text>
              <Text color={theme.palette.text}>{previewOption?.label ?? theme.name}</Text>
            </Box>
            <Box>
              <Text color={theme.palette.textMuted}>Mode: </Text>
              <Text color={previewStatus === 'Saved' ? theme.palette.success : theme.palette.warning}>{previewStatus}</Text>
            </Box>
            <Box marginTop={1}>
              <Text color={theme.palette.accent} bold>Template Menu</Text>
            </Box>
            {TEMPLATE_MENU_ITEMS.map((item, index) => (
              <Box key={item.label}>
                <Text color={index === 0 ? theme.palette.primary : theme.palette.textMuted}>
                  {index === 0 ? '> ' : '  '}
                </Text>
                <Text color={menuTone(item.tone)}>{item.label}</Text>
              </Box>
            ))}
            <Box marginTop={1}>
              <Text color={theme.palette.primary}>Primary</Text>
              <Text color={theme.palette.textMuted}> / </Text>
              <Text color={theme.palette.accent}>Accent</Text>
              <Text color={theme.palette.textMuted}> / </Text>
              <Text color={theme.palette.info}>Info</Text>
            </Box>
            <Box>
              <Text color={theme.palette.success}>Success</Text>
              <Text color={theme.palette.textMuted}> / </Text>
              <Text color={theme.palette.warning}>Warning</Text>
              <Text color={theme.palette.textMuted}> / </Text>
              <Text color={theme.palette.error}>Error</Text>
            </Box>
          </Box>
        </Box>
      </Box>
      <Box flexGrow={1} />
    </Box>
  );
}

export default ThemeMenu;
