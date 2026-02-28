import React from 'react';
import { Box, Text, useInput } from 'ink';
import { useStore } from '../store/index.js';
import { SCREEN_SHORTCUTS, SCREEN_TITLES } from '../screens/registry.js';
import { useTerminalSize } from '../hooks/useTerminalSize.js';
import { useTheme } from '../theme/context.js';
import type { Screen, ScreenShortcut } from '../types/index.js';

// Extended help descriptions for each screen
const SCREEN_HELP_DESCRIPTIONS: Partial<Record<Screen, Record<string, string>>> = {
  welcome: {
    'n': 'Start a new pipeline run with current configuration',
    'c': 'Open configuration settings',
    'h': 'View job history and results',
    '→/Enter': 'Activate selected action or view job',
    'q': 'Exit the application',
    'Esc': 'Exit the application',
  },
  configMain: {
    '↑↓': 'Move between configuration fields',
    'e': 'Edit the selected field value',
    'Ctrl+F': 'Search/filter fields by name',
    's': 'Save all changes to disk',
    '←/Esc': 'Return to previous screen (prompts if unsaved)',
  },
  configDescriptors: {
    '↑↓': 'Move between descriptor settings',
    'e': 'Edit the selected setting',
    'b': 'Edit border min/max values',
    '1/2/3': 'Apply preset (Drug-like/Lead-like/Fragment)',
    's': 'Save configuration',
    '←/Esc': 'Return to configuration menu',
  },
  configFilters: {
    '↑↓': 'Navigate through filter options',
    'Space': 'Enable/disable selected filter',
    'r': 'Switch to rulesets view',
    'Ctrl+F': 'Search rulesets (in rulesets view)',
    'a/n': 'Select all/none rulesets',
    's': 'Save filter settings',
    '←/Esc': 'Return to settings / configuration menu',
  },
  configDocking: {
    '↑↓': 'Navigate docking parameters',
    'PgUp/Dn': 'Jump by page (16 items)',
    'e/Enter': 'Edit parameter or toggle boolean',
    'Ctrl+F': 'Search/filter parameters by name',
    's': 'Save docking configuration',
    '←/Esc': 'Return to welcome screen (prompts if unsaved)',
  },
  configSynthesis: {
    '↑↓': 'Navigate synthesis parameters',
    'e/Enter': 'Edit parameter or toggle boolean',
    's': 'Save synthesis configuration',
    '←/Esc': 'Return to welcome screen (prompts if unsaved)',
  },
  pipelineRunner: {
    'c': 'Cancel running pipeline (with confirmation)',
    '←/Esc': 'Return to welcome (pipeline continues in background)',
  },
  history: {
    '↑↓': 'Navigate job list',
    '→/Enter': 'View detailed results',
    'Ctrl+F': 'Search jobs by name or status',
    'd': 'Delete selected job (with confirmation)',
    '←/Esc': 'Return to welcome screen',
  },
  results: {
    'r': 'Re-run this job with same settings',
    '←/Esc': 'Return to history',
  },
  wizardInputSelection: {
    '↑↓': 'Move between input/output fields',
    'Space': 'Browse selected field (in file picker: search for files or select folder)',
    '→/e': 'Open selected path field directly in path edit mode',
    'Enter': 'Save and continue to stage selection',
    '←/Esc': 'Return to welcome screen',
  },
  wizardStageSelection: {
    '↑↓': 'Move between pipeline stages',
    'Space': 'Enable/disable selected stage',
    'c': 'Open configuration for selected stage',
    'Enter': 'Open detailed review screen',
    'r/→': 'Open detailed review screen',
    '←/Esc': 'Return to input selection',
  },
  wizardReview: {
    '↑↓': 'Move between selected stages',
    'e': 'Open config screen for focused stage',
    'r': 'Refresh preflight checks',
    'Tab': 'Switch summary/detailed stage view',
    'Enter': 'Run preflight and start pipeline',
    '←/Esc': 'Return to stage selection',
  },
};

// Global shortcuts always available
const GLOBAL_SHORTCUTS: ScreenShortcut[] = [
  { key: '?', label: 'Toggle help overlay' },
  { key: '/', label: 'Open command palette' },
  { key: 'Ctrl+F', label: 'Search/filter in list screens' },
];

export function HelpOverlay(): React.ReactElement | null {
  const { width, height } = useTerminalSize();
  const { theme } = useTheme();
  const showHelp = useStore((state) => state.showHelp);
  const setShowHelp = useStore((state) => state.setShowHelp);
  const screen = useStore((state) => state.screen);

  useInput(() => {
    if (showHelp) {
      setShowHelp(false);
    }
  }, { isActive: showHelp });

  if (!showHelp) {
    return null;
  }

  const screenShortcuts = SCREEN_SHORTCUTS[screen] || [];
  const helpDescriptions = SCREEN_HELP_DESCRIPTIONS[screen] || {};
  const title = SCREEN_TITLES[screen] || 'Help';

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
          width={Math.max(1, Math.min(86, width - 4))}
          borderStyle="double"
          borderColor={theme.palette.borderActive}
          paddingX={2}
          paddingY={1}
        >
          <Box justifyContent="center" marginBottom={1}>
            <Text color={theme.palette.primary} bold>Help: {title}</Text>
          </Box>

          <Box marginBottom={1}>
            <Text color={theme.palette.accent} bold>Screen Shortcuts</Text>
          </Box>

          {screenShortcuts.map((shortcut, index) => (
            <Box key={index}>
              <Text color={theme.palette.primary} bold>{shortcut.key.padEnd(12)}</Text>
              <Text color={theme.palette.text}>
                {helpDescriptions[shortcut.key] || shortcut.label}
              </Text>
            </Box>
          ))}

          <Box marginTop={1} marginBottom={1}>
            <Text color={theme.palette.accent} bold>Global Shortcuts</Text>
          </Box>

          {GLOBAL_SHORTCUTS.map((shortcut, index) => (
            <Box key={`global-${index}`}>
              <Text color={theme.palette.primary} bold>{shortcut.key.padEnd(12)}</Text>
              <Text color={theme.palette.text}>{shortcut.label}</Text>
            </Box>
          ))}

          <Box marginTop={1} justifyContent="center">
            <Text color={theme.palette.textMuted}>Press any key to close</Text>
          </Box>
        </Box>
      </Box>
      <Box flexGrow={1} />
    </Box>
  );
}

export default HelpOverlay;
