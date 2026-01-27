import React from 'react';
import { Box, Text, useInput } from 'ink';
import { useStore } from '../store/index.js';
import { SCREEN_SHORTCUTS, SCREEN_TITLES } from '../App.js';
import type { Screen, ScreenShortcut } from '../types/index.js';

// Extended help descriptions for each screen
const SCREEN_HELP_DESCRIPTIONS: Partial<Record<Screen, Record<string, string>>> = {
  welcome: {
    'n': 'Start a new pipeline run with current configuration',
    'c': 'Open configuration settings',
    'h': 'View job history and results',
    '→/Enter': 'Activate selected action or view job',
    '←': 'Exit the application',
    'q': 'Exit the application',
  },
  configMain: {
    '↑↓': 'Move between configuration fields',
    'e': 'Edit the selected field value',
    '/': 'Search/filter fields by name',
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
    '/': 'Search rulesets (in rulesets view)',
    'a/n': 'Select all/none rulesets',
    's': 'Save filter settings',
    '←/Esc': 'Return to settings / configuration menu',
  },
  configDocking: {
    '↑↓': 'Navigate docking parameters',
    'PgUp/Dn': 'Jump by page (16 items)',
    'e/Enter': 'Edit parameter or toggle boolean',
    '/': 'Search/filter parameters by name',
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
    '/': 'Search jobs by name or status',
    'd': 'Delete selected job (with confirmation)',
    '←/Esc': 'Return to welcome screen',
  },
  results: {
    'r': 'Re-run this job with same settings',
    '←/Esc': 'Return to history',
  },
};

// Global shortcuts always available
const GLOBAL_SHORTCUTS: ScreenShortcut[] = [
  { key: '?', label: 'Toggle help overlay' },
  { key: '/', label: 'Search/filter (in lists)' },
];

export function HelpOverlay(): React.ReactElement | null {
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
      borderStyle="double"
      borderColor="magenta"
      paddingX={2}
      paddingY={1}
      marginY={1}
    >
      <Box justifyContent="center" marginBottom={1}>
        <Text color="cyan" bold>Help: {title}</Text>
      </Box>

      <Box marginBottom={1}>
        <Text color="yellow" bold>Screen Shortcuts</Text>
      </Box>

      {screenShortcuts.map((shortcut, index) => (
        <Box key={index}>
          <Text color="cyan" bold>{shortcut.key.padEnd(12)}</Text>
          <Text>
            {helpDescriptions[shortcut.key] || shortcut.label}
          </Text>
        </Box>
      ))}

      <Box marginTop={1} marginBottom={1}>
        <Text color="yellow" bold>Global Shortcuts</Text>
      </Box>

      {GLOBAL_SHORTCUTS.map((shortcut, index) => (
        <Box key={`global-${index}`}>
          <Text color="cyan" bold>{shortcut.key.padEnd(12)}</Text>
          <Text>{shortcut.label}</Text>
        </Box>
      ))}

      <Box marginTop={1} justifyContent="center">
        <Text dimColor>Press any key to close</Text>
      </Box>
    </Box>
  );
}

export default HelpOverlay;
