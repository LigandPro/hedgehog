import type { Screen, ScreenShortcut } from '../types/index.js';

// Central screen metadata registry.
// Keep this file dependency-free (no imports from App/screens/components) to avoid cycles.

export const SCREEN_TITLES: Record<Screen, string> = {
  welcome: '',
  configMain: 'Main Configuration',
  configDescriptors: 'Descriptors Settings',
  configFilters: 'Structure Filters',
  configSynthesis: 'Synthesis Scoring',
  configRetrosynthesis: 'Retrosynthesis Config',
  configDocking: 'Docking Configuration',
  pipelineRunner: 'Pipeline Runner',
  history: 'Job History',
  results: 'Job Results',
  // Wizard screens
  wizardInputSelection: 'Pipeline Wizard',
  wizardStageSelection: 'Pipeline Wizard',
  wizardStageOrder: 'Pipeline Wizard',
  wizardConfigMolPrep: 'Pipeline Wizard',
  wizardConfigDescriptors: 'Pipeline Wizard',
  wizardConfigFilters: 'Pipeline Wizard',
  wizardConfigSynthesis: 'Pipeline Wizard',
  wizardConfigDocking: 'Pipeline Wizard',
  wizardConfigDockingFilters: 'Pipeline Wizard',
  wizardReview: 'Pipeline Wizard',
};

export const SCREEN_SHORTCUTS: Record<Screen, ScreenShortcut[]> = {
  welcome: [
    { key: 'n', label: 'New Run' },
    { key: 'h', label: 'History' },
    { key: 'c', label: 'Config' },
    { key: '→/Enter', label: 'Select' },
    { key: 'q', label: 'Quit' },
    { key: 'Esc', label: 'Quit' },
  ],
  configMain: [
    { key: '↑↓', label: 'Navigate' },
    { key: 'e/Enter', label: 'Edit' },
    { key: 's', label: 'Save' },
    { key: '←/Esc', label: 'Back' },
  ],
  configDescriptors: [
    { key: '↑↓', label: 'Navigate' },
    { key: 'e/Enter', label: 'Edit' },
    { key: 'b', label: 'Borders' },
    { key: 's', label: 'Save' },
    { key: '←/Esc', label: 'Back' },
  ],
  configFilters: [
    { key: '↑↓', label: 'Navigate' },
    { key: 'Space', label: 'Toggle' },
    { key: 'r', label: 'Rulesets' },
    { key: 's', label: 'Save' },
    { key: '←/Esc', label: 'Back' },
  ],
  configSynthesis: [
    { key: '↑↓', label: 'Navigate' },
    { key: 'e/Enter', label: 'Edit' },
    { key: 'r', label: 'Retrosynthesis' },
    { key: 's', label: 'Save' },
    { key: '←/Esc', label: 'Back' },
  ],
  configRetrosynthesis: [
    { key: '↑↓', label: 'Navigate' },
    { key: 'b/Enter', label: 'Browse' },
    { key: 'Ctrl+F', label: 'Search' },
    { key: 's', label: 'Save' },
    { key: '←/Esc', label: 'Back' },
  ],
  configDocking: [
    { key: '↑↓', label: 'Navigate' },
    { key: 'PgUp/Dn', label: 'Page' },
    { key: 'e/Enter', label: 'Edit' },
    { key: 's', label: 'Save' },
    { key: '←/Esc', label: 'Back' },
  ],
  pipelineRunner: [
    { key: 'c', label: 'Cancel' },
    { key: 'l', label: 'Show/Hide Log' },
    { key: '←/Esc', label: 'Back' },
  ],
  history: [
    { key: '↑↓', label: 'Navigate' },
    { key: '→/Enter', label: 'View' },
    { key: 'd', label: 'Delete' },
    { key: '←/Esc', label: 'Back' },
  ],
  results: [
    { key: 'r', label: 'Re-run' },
    { key: '←/Esc', label: 'Back' },
  ],
  // Wizard screens
  wizardInputSelection: [
    { key: '↑↓', label: 'Navigate' },
    { key: 'Space', label: 'Edit/Browse' },
    { key: '→/e', label: 'Edit path' },
    { key: 'Enter', label: 'Next' },
    { key: '←/Esc', label: 'Back' },
  ],
  wizardStageSelection: [
    { key: 'Space', label: 'Toggle' },
    { key: 'c', label: 'Configure stage' },
    { key: 'Enter', label: 'Review' },
    { key: 'r/→', label: 'Review' },
    { key: '←/Esc', label: 'Back' },
  ],
  wizardStageOrder: [
    { key: '→/Enter', label: 'Next' },
    { key: '←/Esc', label: 'Back' },
  ],
  wizardConfigMolPrep: [
    { key: '↑↓', label: 'Navigate' },
    { key: 'Space', label: 'Edit' },
    { key: 's/Enter', label: 'Save' },
    { key: '←/Esc', label: 'Back' },
  ],
  wizardConfigDescriptors: [
    { key: '↑↓', label: 'Navigate' },
    { key: 'Space', label: 'Edit' },
    { key: 's/Enter', label: 'Save' },
    { key: '←/Esc', label: 'Back' },
  ],
  wizardConfigFilters: [
    { key: '↑↓', label: 'Navigate' },
    { key: 'Space', label: 'Edit' },
    { key: 's/Enter', label: 'Save' },
    { key: '←/Esc', label: 'Back' },
  ],
  wizardConfigSynthesis: [
    { key: '↑↓', label: 'Navigate' },
    { key: 'Space', label: 'Edit' },
    { key: 'r', label: 'Retrosynthesis' },
    { key: 's/Enter', label: 'Save' },
    { key: '←/Esc', label: 'Back' },
  ],
  wizardConfigDocking: [
    { key: '↑↓', label: 'Navigate' },
    { key: 'Space', label: 'Edit' },
    { key: 's/Enter', label: 'Save' },
    { key: '←/Esc', label: 'Back' },
  ],
  wizardConfigDockingFilters: [
    { key: '↑↓', label: 'Navigate' },
    { key: 'Space', label: 'Edit' },
    { key: 's/Enter', label: 'Save' },
    { key: '←/Esc', label: 'Back' },
  ],
  wizardReview: [
    { key: '↑↓', label: 'Navigate' },
    { key: 'e', label: 'Edit stage' },
    { key: 'r', label: 'Refresh preflight' },
    { key: 'Tab', label: 'Summary/Detailed' },
    { key: 'Enter', label: 'Start' },
    { key: '←/Esc', label: 'Back' },
  ],
};
