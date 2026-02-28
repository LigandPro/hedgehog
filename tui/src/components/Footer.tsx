import React from 'react';
import { Box, Text } from 'ink';
import { useStore } from '../store/index.js';
import { SCREEN_SHORTCUTS } from '../screens/registry.js';
import { useTerminalSize } from '../hooks/useTerminalSize.js';
import { useTheme } from '../theme/context.js';
import { APP_VERSION } from '../constants/version.js';
import type { Screen, ScreenShortcut } from '../types/index.js';

// Breadcrumb paths for each screen
const SCREEN_BREADCRUMBS: Record<Screen, string[]> = {
  welcome: ['Home'],
  configMain: ['Home', 'Configuration'],
  configDescriptors: ['Home', 'Configuration', 'Descriptors'],
  configFilters: ['Home', 'Configuration', 'Filters'],
  configSynthesis: ['Home', 'Configuration', 'Synthesis'],
  configRetrosynthesis: ['Home', 'Configuration', 'Synthesis', 'Retrosynthesis'],
  configDocking: ['Home', 'Configuration', 'Docking'],
  pipelineRunner: ['Home', 'Pipeline'],
  history: ['Home', 'History'],
  results: ['Home', 'History', 'Results'],
  // Wizard screens
  wizardInputSelection: ['Home', 'Wizard', 'Input'],
  wizardStageSelection: ['Home', 'Wizard', 'Stages'],
  wizardStageOrder: ['Home', 'Wizard', 'Order'],
  wizardConfigMolPrep: ['Home', 'Wizard', 'Mol Prep'],
  wizardConfigDescriptors: ['Home', 'Wizard', 'Descriptors'],
  wizardConfigFilters: ['Home', 'Wizard', 'Filters'],
  wizardConfigSynthesis: ['Home', 'Wizard', 'Synthesis'],
  wizardConfigDocking: ['Home', 'Wizard', 'Docking'],
  wizardConfigDockingFilters: ['Home', 'Wizard', 'Docking Filters'],
  wizardReview: ['Home', 'Wizard', 'Review'],
};

function fitLineEnd(text: string, width: number): string {
  if (width <= 0) return '';
  if (text.length === width) return text;
  if (text.length < width) return text.padEnd(width, ' ');
  if (width <= 3) return text.slice(0, width);
  return `${text.slice(0, width - 3)}...`;
}

function fitBreadcrumbs(crumbs: string[], width: number): string {
  const full = crumbs.join(' > ');
  if (full.length <= width) return full.padEnd(width, ' ');
  if (width <= 3) return full.slice(0, width);

  const last = crumbs.at(-1) ?? '';
  const tail2 = crumbs.slice(-2).join(' > ');
  const compact2 = `... ${tail2}`;
  if (compact2.length <= width) return compact2.padEnd(width, ' ');

  const compact1 = `... ${last}`;
  return fitLineEnd(compact1, width);
}

function fitRight(text: string, width: number): string {
  if (width <= 0) return '';
  if (text.length >= width) return text.slice(0, width);
  return `${' '.repeat(width - text.length)}${text}`;
}

interface FooterProps {
  shortcuts?: ScreenShortcut[];
  overrideShortcuts?: boolean;
  showBreadcrumbs?: boolean;
}

interface ShortcutToken {
  key: string;
  label: string;
}

function getVisibleShortcutCount(tokens: ShortcutToken[], width: number): number {
  if (width <= 0) return 0;
  let used = 0;
  let count = 0;

  for (const token of tokens) {
    const tokenLen = token.key.length + 1 + token.label.length;
    const sepLen = count > 0 ? 2 : 0;
    if (used + sepLen + tokenLen > width) break;
    used += sepLen + tokenLen;
    count += 1;
  }

  return count;
}

export function Footer({
  shortcuts,
  overrideShortcuts = false,
  showBreadcrumbs = true,
}: FooterProps): React.ReactElement {
  const { theme } = useTheme();
  const screen = useStore((state) => state.screen);
  const debugMode = useStore((state) => state.debugMode);
  const isBackendReady = useStore((state) => state.isBackendReady);
  const isRunning = useStore((state) => state.isRunning);
  const searchActive = useStore((state) => state.searchActive);
  const searchQuery = useStore((state) => state.searchQuery);

  const { width: terminalWidth } = useTerminalSize();

  // Use provided shortcuts or fall back to screen-based shortcuts
  const displayShortcuts = overrideShortcuts && shortcuts
    ? shortcuts
    : (shortcuts || SCREEN_SHORTCUTS[screen] || []);

  const breadcrumbs = SCREEN_BREADCRUMBS[screen] || ['Home'];
  const separator = '─'.repeat(Math.max(0, terminalWidth));
  const shortcutTokens: ShortcutToken[] = displayShortcuts.map((shortcut) => ({
    key: shortcut.key,
    label: shortcut.label,
  }));
  if (
    ['configMain', 'configDescriptors', 'configDocking', 'configFilters', 'configSynthesis', 'configRetrosynthesis', 'history'].includes(screen)
    && !searchActive
  ) {
    shortcutTokens.push({ key: 'ctrl+f', label: 'search' });
  }
  shortcutTokens.push({ key: 'ctrl+p', label: 'commands' });
  const statusColor = isRunning
    ? theme.palette.warning
    : (isBackendReady ? theme.palette.success : theme.palette.textMuted);

  const rightStatus = ` ●${debugMode ? ' [D]' : ''}`;
  const rightWidth = rightStatus.length;
  const leftWidth = Math.max(0, terminalWidth - rightWidth);
  const visibleShortcutCount = getVisibleShortcutCount(shortcutTokens, leftWidth);
  const visibleShortcutTokens = shortcutTokens.slice(0, visibleShortcutCount);
  const hasOverflow = visibleShortcutCount < shortcutTokens.length;

  return (
    <Box flexDirection="column" marginTop={1}>
      {/* Search bar when active */}
      {searchActive && (
        <Box marginBottom={1} width={terminalWidth}>
          <Text color={theme.palette.primary}>
            {fitLineEnd(`Search: ${searchQuery}█ (Esc to cancel, Enter to confirm)`, terminalWidth)}
          </Text>
        </Box>
      )}

      {/* Breadcrumbs */}
      {showBreadcrumbs && !searchActive && (
        <Box marginBottom={0} width={terminalWidth}>
          <Text color={theme.palette.textMuted}>{fitBreadcrumbs(breadcrumbs, terminalWidth)}</Text>
        </Box>
      )}

      {/* Separator line */}
      <Box width={terminalWidth}>
        <Text color={theme.palette.border} wrap="truncate">{separator}</Text>
      </Box>

      <Box width={terminalWidth}>
        <Box width={leftWidth}>
          {visibleShortcutTokens.length > 0 ? (
            <>
              {visibleShortcutTokens.map((token, index) => (
                <React.Fragment key={`${token.key}-${token.label}-${index}`}>
                  {index > 0 && <Text color={theme.palette.textMuted}>  </Text>}
                  <Text color="white">{token.key}</Text>
                  <Text color={theme.palette.textMuted}> {token.label}</Text>
                </React.Fragment>
              ))}
              {hasOverflow && <Text color={theme.palette.textMuted}> ...</Text>}
            </>
          ) : (
            <Text color={theme.palette.textMuted}>{fitLineEnd('', leftWidth)}</Text>
          )}
        </Box>
        <Text color={statusColor}> ●</Text>
        {debugMode && <Text color={theme.palette.warning} bold> [D]</Text>}
      </Box>

      <Box width={terminalWidth} justifyContent="flex-end">
        <Text color={theme.palette.textMuted}>{fitRight(APP_VERSION, terminalWidth)}</Text>
      </Box>
    </Box>
  );
}

export default Footer;
