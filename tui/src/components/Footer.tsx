import React from 'react';
import { Box, Text } from 'ink';
import { useStore } from '../store/index.js';
import { SCREEN_SHORTCUTS } from '../App.js';
import { useTerminalSize } from '../hooks/useTerminalSize.js';
import type { Screen, ScreenShortcut } from '../types/index.js';

// Breadcrumb paths for each screen
const SCREEN_BREADCRUMBS: Record<Screen, string[]> = {
  welcome: ['Home'],
  configMain: ['Home', 'Configuration'],
  configDescriptors: ['Home', 'Configuration', 'Descriptors'],
  configFilters: ['Home', 'Configuration', 'Filters'],
  configSynthesis: ['Home', 'Configuration', 'Synthesis'],
  configDocking: ['Home', 'Configuration', 'Docking'],
  pipelineRunner: ['Home', 'Pipeline'],
  history: ['Home', 'History'],
  results: ['Home', 'History', 'Results'],
  // Wizard screens
  wizardInputSelection: ['Home', 'Wizard', 'Input'],
  wizardStageSelection: ['Home', 'Wizard', 'Stages'],
  wizardStageOrder: ['Home', 'Wizard', 'Order'],
  wizardConfigDescriptors: ['Home', 'Wizard', 'Descriptors'],
  wizardConfigFilters: ['Home', 'Wizard', 'Filters'],
  wizardConfigSynthesis: ['Home', 'Wizard', 'Synthesis'],
  wizardConfigDocking: ['Home', 'Wizard', 'Docking'],
  wizardReview: ['Home', 'Wizard', 'Review'],
};

interface FooterProps {
  shortcuts?: ScreenShortcut[];
  overrideShortcuts?: boolean;
  showBreadcrumbs?: boolean;
}

export function Footer({
  shortcuts,
  overrideShortcuts = false,
  showBreadcrumbs = true,
}: FooterProps): React.ReactElement {
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

  return (
    <Box flexDirection="column" marginTop={1}>
      {/* Search bar when active */}
      {searchActive && (
        <Box marginBottom={1}>
          <Text color="cyan">Search: </Text>
          <Text color="yellow">{searchQuery}</Text>
          <Text color="cyan">█</Text>
          <Text dimColor> (Esc to cancel, Enter to confirm)</Text>
        </Box>
      )}

      {/* Breadcrumbs */}
      {showBreadcrumbs && !searchActive && (
        <Box marginBottom={0}>
          {breadcrumbs.map((crumb, index) => (
            <React.Fragment key={index}>
              {index > 0 && <Text dimColor> {'>'} </Text>}
              <Text color={index === breadcrumbs.length - 1 ? 'cyan' : 'gray'}>
                {crumb}
              </Text>
            </React.Fragment>
          ))}
        </Box>
      )}

      {/* Separator line */}
      <Box>
        <Text color="gray">{'─'.repeat(terminalWidth - 2)}</Text>
      </Box>

      {/* Shortcuts as simple text */}
      <Text wrap="wrap">
        {displayShortcuts.map((shortcut, index) => (
          <Text key={index}>
            {shortcut.disabled ? (
              <Text color="gray" dimColor>[{shortcut.key}] {shortcut.label}</Text>
            ) : (
              <>
                <Text color="cyan" bold>[{shortcut.key}]</Text>
                <Text> {shortcut.label}</Text>
              </>
            )}
            <Text>  </Text>
          </Text>
        ))}

        {/* Search shortcut hint in list screens */}
        {['configMain', 'configDescriptors', 'configDocking', 'configFilters', 'configSynthesis', 'history'].includes(screen) && !searchActive && (
          <Text>
            <Text color="cyan" bold>[/]</Text>
            <Text> Search  </Text>
          </Text>
        )}

        {/* Status indicator */}
        {isRunning ? (
          <Text color="yellow">●</Text>
        ) : (
          <Text color={isBackendReady ? 'green' : 'red'}>
            {isBackendReady ? '●' : '○'}
          </Text>
        )}
        {debugMode && <Text color="yellow" bold> [D]</Text>}
      </Text>
    </Box>
  );
}

export default Footer;
