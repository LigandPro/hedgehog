import React, { useEffect, useRef } from 'react';
import { Box, Text, useApp, useInput } from 'ink';
import { useStore } from './store/index.js';
import { getBridge, destroyBridge } from './services/python-bridge.js';
import { logger } from './utils/logger.js';
import type { Screen, ScreenShortcut, JobHistoryRecord } from './types/index.js';

// Screens
import { Welcome } from './screens/Welcome.js';
import { ConfigMain } from './screens/ConfigMain.js';
import { ConfigDescriptors } from './screens/ConfigDescriptors.js';
import { ConfigFilters } from './screens/ConfigFilters.js';
import { ConfigSynthesis } from './screens/ConfigSynthesis.js';
import { ConfigRetrosynthesis } from './screens/ConfigRetrosynthesis.js';
import { ConfigDocking } from './screens/ConfigDocking.js';
import { PipelineRunner } from './screens/PipelineRunner.js';
import { History } from './screens/History.js';
import { Results } from './screens/Results.js';

// Wizard screens
import {
  InputSelection,
  StageSelection,
  StageOrder,
  WizardConfigMolPrep,
  WizardConfigDescriptors,
  WizardConfigFilters,
  WizardConfigSynthesis,
  WizardConfigDocking,
  WizardConfigDockingFilters,
  ReviewRun,
} from './screens/wizard/index.js';

// Global overlays
import { ToastContainer } from './components/Toast.js';
import { ConfirmDialog } from './components/ConfirmDialog.js';
import { HelpOverlay } from './components/HelpOverlay.js';

// Back navigation map (Matcha pattern)
export const BACK_MAP: Record<Screen, Screen | null> = {
  welcome: null,
  configMain: 'welcome',
  configDescriptors: 'configMain',
  configFilters: 'configMain',
  configSynthesis: 'configMain',
  configRetrosynthesis: 'configSynthesis',
  configDocking: 'configMain',
  pipelineRunner: 'welcome',
  history: 'welcome',
  results: 'history',
  // Wizard screens - navigation handled internally
  wizardInputSelection: 'welcome',
  wizardStageSelection: 'wizardInputSelection',
  wizardStageOrder: 'wizardStageSelection',
  wizardConfigMolPrep: 'wizardStageSelection',
  wizardConfigDescriptors: 'wizardStageSelection',
  wizardConfigFilters: 'wizardStageSelection',
  wizardConfigSynthesis: 'wizardStageSelection',
  wizardConfigDocking: 'wizardStageSelection',
  wizardConfigDockingFilters: 'wizardStageSelection',
  wizardReview: 'wizardStageSelection',
};

// Screen titles map (Matcha pattern)
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

// Screen shortcuts map (Matcha pattern)
export const SCREEN_SHORTCUTS: Record<Screen, ScreenShortcut[]> = {
  welcome: [
    { key: 'n', label: 'New Run' },
    { key: 'h', label: 'History' },
    { key: '→/Enter', label: 'Select' },
    { key: 'q', label: 'Quit' },
  ],
  configMain: [
    { key: '↑↓', label: 'Navigate' },
    { key: 'e', label: 'Edit' },
    { key: 's', label: 'Save' },
    { key: '←/Esc', label: 'Back' },
  ],
  configDescriptors: [
    { key: '↑↓', label: 'Navigate' },
    { key: 'e', label: 'Edit' },
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
    { key: 'e', label: 'Edit' },
    { key: 'r', label: 'Retrosynthesis' },
    { key: 's', label: 'Save' },
    { key: '←/Esc', label: 'Back' },
  ],
  configRetrosynthesis: [
    { key: '↑↓', label: 'Navigate' },
    { key: 'b/Enter', label: 'Browse' },
    { key: '/', label: 'Search' },
    { key: 's', label: 'Save' },
    { key: '←/Esc', label: 'Back' },
  ],
  configDocking: [
    { key: '↑↓', label: 'Navigate' },
    { key: 'e', label: 'Edit' },
    { key: 's', label: 'Save' },
    { key: '←/Esc', label: 'Back' },
  ],
  pipelineRunner: [
    { key: 'c', label: 'Cancel' },
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
  // Wizard screens - shortcuts handled internally
  wizardInputSelection: [
    { key: '↑↓', label: 'Navigate' },
    { key: 'Enter', label: 'Browse' },
    { key: '→', label: 'Next' },
    { key: '←/Esc', label: 'Back' },
  ],
  wizardStageSelection: [
    { key: 'Space', label: 'Toggle' },
    { key: 'a/n', label: 'All/None' },
    { key: '→/Enter', label: 'Next' },
    { key: '←/Esc', label: 'Back' },
  ],
  wizardStageOrder: [
    { key: 'J/K', label: 'Move' },
    { key: 'Space', label: 'Toggle dep' },
    { key: '→/Enter', label: 'Next' },
    { key: '←/Esc', label: 'Back' },
  ],
  wizardConfigMolPrep: [
    { key: '↑↓', label: 'Navigate' },
    { key: 'Space', label: 'Edit' },
    { key: 's', label: 'Save' },
    { key: '←/→', label: 'Prev/Next' },
  ],
  wizardConfigDescriptors: [
    { key: '↑↓', label: 'Navigate' },
    { key: 'Enter', label: 'Edit' },
    { key: 's', label: 'Save' },
    { key: '1-3', label: 'Presets' },
    { key: '←/→', label: 'Prev/Next' },
  ],
  wizardConfigFilters: [
    { key: '↑↓', label: 'Navigate' },
    { key: 'Space', label: 'Toggle' },
    { key: 's', label: 'Save' },
    { key: '←/→', label: 'Prev/Next' },
  ],
  wizardConfigSynthesis: [
    { key: '↑↓', label: 'Navigate' },
    { key: 'Enter', label: 'Edit' },
    { key: 's', label: 'Save' },
    { key: '←/→', label: 'Prev/Next' },
  ],
  wizardConfigDocking: [
    { key: '↑↓', label: 'Navigate' },
    { key: 'Enter', label: 'Edit' },
    { key: 's', label: 'Save' },
    { key: '←/→', label: 'Prev/Next' },
  ],
  wizardConfigDockingFilters: [
    { key: '↑↓', label: 'Navigate' },
    { key: 'Enter', label: 'Edit' },
    { key: 's', label: 'Save' },
    { key: '←/→', label: 'Prev/Next' },
  ],
  wizardReview: [
    { key: '↑↓', label: 'Navigate' },
    { key: 'e', label: 'Edit stage' },
    { key: 'Enter', label: 'Start' },
    { key: '←/Esc', label: 'Back' },
  ],
};

const STAGE_NAMES: Record<string, string> = {
  mol_prep: 'Mol Prep',
  descriptors: 'Descriptors',
  struct_filters: 'Struct Filters',
  synthesis: 'Synthesis',
  docking: 'Docking',
  docking_filters: 'Docking Filters',
};

const PROGRESS_LOG_THRESHOLDS = [10, 25, 50, 75, 90, 100] as const;

function getStageName(stage: string): string {
  return STAGE_NAMES[stage] || stage;
}

// ScreenRouter component (Matcha pattern)
function ScreenRouter({ screen }: { screen: Screen }): React.ReactElement {
  // Each screen gets a unique key to force remount on screen change
  const content = (() => {
    switch (screen) {
      case 'welcome':
        return <Welcome />;
      case 'configMain':
        return <ConfigMain />;
      case 'configDescriptors':
        return <ConfigDescriptors />;
      case 'configFilters':
        return <ConfigFilters />;
      case 'configSynthesis':
        return <ConfigSynthesis />;
      case 'configRetrosynthesis':
        return <ConfigRetrosynthesis />;
      case 'configDocking':
        return <ConfigDocking />;
      case 'pipelineRunner':
        return <PipelineRunner />;
      case 'history':
        return <History />;
      case 'results':
        return <Results />;
      // Wizard screens
      case 'wizardInputSelection':
        return <InputSelection />;
      case 'wizardStageSelection':
        return <StageSelection />;
      case 'wizardStageOrder':
        return <StageOrder />;
      case 'wizardConfigMolPrep':
        return <WizardConfigMolPrep />;
      case 'wizardConfigDescriptors':
        return <WizardConfigDescriptors />;
      case 'wizardConfigFilters':
        return <WizardConfigFilters />;
      case 'wizardConfigSynthesis':
        return <WizardConfigSynthesis />;
      case 'wizardConfigDocking':
        return <WizardConfigDocking />;
      case 'wizardConfigDockingFilters':
        return <WizardConfigDockingFilters />;
      case 'wizardReview':
        return <ReviewRun />;
      default:
        return <Welcome />;
    }
  })();

  // Wrap in Box with key to force clean remount on screen change
  return <Box key={screen} flexDirection="column">{content}</Box>;
}

export function App(): React.ReactElement {
  useApp();
  const progressLogThresholdRef = useRef<Record<string, number>>({});
  const screen = useStore((state) => state.screen);
  const setBackendReady = useStore((state) => state.setBackendReady);
  const isBackendReady = useStore((state) => state.isBackendReady);
  const setJobHistory = useStore((state) => state.setJobHistory);
  const globalError = useStore((state) => state.globalError);
  const setGlobalError = useStore((state) => state.setGlobalError);
  const debugMode = useStore((state) => state.debugMode);
  const setDebugMode = useStore((state) => state.setDebugMode);
  const showHelp = useStore((state) => state.showHelp);
  const setShowHelp = useStore((state) => state.setShowHelp);
  const confirmDialog = useStore((state) => state.confirmDialog);
  const setRunning = useStore((state) => state.setRunning);
  const updateStage = useStore((state) => state.updateStage);
  const updatePipelineProgress = useStore((state) => state.updatePipelineProgress);
  const addLog = useStore((state) => state.addLog);
  const updateJobInHistory = useStore((state) => state.updateJobInHistory);
  const showToast = useStore((state) => state.showToast);
  const setSelectedJob = useStore((state) => state.setSelectedJob);
  const setScreen = useStore((state) => state.setScreen);

  // Initialize Python backend
  useEffect(() => {
    const initBackend = async () => {
      try {
        const bridge = getBridge();
        await bridge.start();
        setBackendReady(true);
        logger.info('Backend initialized');

        // Load job history from backend
        try {
          const history = await bridge.getJobHistory();
          setJobHistory(history as unknown as JobHistoryRecord[]);
          logger.info(`Loaded ${history.length} jobs from history`);
        } catch (historyError) {
          logger.warn('Failed to load job history:', historyError);
        }
      } catch (error) {
        logger.error('Failed to initialize backend:', error);
        setGlobalError('Backend connection failed');
      }
    };

    initBackend();

    return () => {
      destroyBridge();
    };
  }, []);

  // Global notification handler so pipeline can run in background across screens.
  useEffect(() => {
    if (!isBackendReady) return;
    const bridge = getBridge();
    const unsubscribe = bridge.onNotification((method, params) => {
      const state = useStore.getState();
      const stageOrder = state.wizard.stageOrder;
      const selectedStages = stageOrder.filter((s) => state.wizard.selectedStages.includes(s));

      if (method === 'progress') {
        const { stage, current, total, message } = params as any;
        const progress = total > 0 ? Math.round((current / total) * 100) : 0;
        updateStage(stage, {
          progress,
          status: progress === 100 ? 'completed' : 'running',
          message,
        });
        const stageIndex = selectedStages.indexOf(stage) + 1;
        updatePipelineProgress({
          currentStage: stage,
          stageIndex: stageIndex > 0 ? stageIndex : 0,
          totalStages: selectedStages.length,
          stageProgress: progress,
          latestMessage: message,
        });

        if (stage) {
          const lastLoggedThreshold = progressLogThresholdRef.current[stage] ?? 0;
          let reachedThreshold: number | null = null;

          for (const threshold of PROGRESS_LOG_THRESHOLDS) {
            if (progress >= threshold && threshold > lastLoggedThreshold) {
              reachedThreshold = threshold;
            }
          }

          if (reachedThreshold !== null) {
            progressLogThresholdRef.current[stage] = reachedThreshold;
            addLog({
              timestamp: new Date(),
              level: 'info',
              message: `Stage progress: ${getStageName(stage)} ${reachedThreshold}%`,
            });
          }
        }
        return;
      }

      if (method === 'log') {
        if (!debugMode) return;
        addLog({
          timestamp: new Date(),
          level: (params as any).level || 'info',
          message: (params as any).message,
        });
        return;
      }

      if (method === 'stage_start') {
        const stage = (params as any).stage;
        progressLogThresholdRef.current[stage] = 0;
        updateStage(stage, { status: 'running', progress: 0 });
        const stageIndex = selectedStages.indexOf(stage) + 1;
        updatePipelineProgress({
          currentStage: stage,
          stageIndex: stageIndex > 0 ? stageIndex : 0,
          totalStages: selectedStages.length,
          stageProgress: 0,
        });
        addLog({
          timestamp: new Date(),
          level: 'info',
          message: `Stage started: ${getStageName(stage)}`,
        });
        return;
      }

      if (method === 'stage_complete') {
        const stage = (params as any).stage;
        progressLogThresholdRef.current[stage] = 100;
        updateStage(stage, { status: 'completed', progress: 100 });
        addLog({
          timestamp: new Date(),
          level: 'info',
          message: `Stage completed: ${getStageName(stage)}`,
        });
        return;
      }

      if (method === 'stage_error') {
        const stage = (params as any).stage;
        const message = (params as any).message || 'Unknown error';
        updateStage(stage, { status: 'error', message });
        addLog({
          timestamp: new Date(),
          level: 'error',
          message: `Stage failed: ${getStageName(stage)} - ${message}`,
        });
        return;
      }

      if (method === 'complete') {
        setRunning(false);
        const currentJobId = state.currentJobId;
        if (currentJobId) {
          setSelectedJob(currentJobId);
        }
        showToast('success', 'Pipeline completed');
        if (currentJobId) {
          const results = (params as any).results || {};
          updateJobInHistory(currentJobId, {
            status: 'completed',
            endTime: new Date().toISOString(),
            results: {
              moleculesProcessed: results.molecules_processed || 0,
            },
          });
          bridge.updateJob(currentJobId, 'completed', {
            moleculesProcessed: results.molecules_processed || 0,
          });
        }
        if (state.screen === 'pipelineRunner') {
          setScreen('results');
        }
        progressLogThresholdRef.current = {};
        return;
      }

      if (method === 'error') {
        const message = (params as any).message;
        showToast('error', message);
        setRunning(false);
        progressLogThresholdRef.current = {};
        const currentJobId = state.currentJobId;
        if (currentJobId) {
          updateJobInHistory(currentJobId, {
            status: 'error',
            endTime: new Date().toISOString(),
            error: message,
          });
          bridge.updateJob(currentJobId, 'error', undefined, message);
        }
      }
    });

    return unsubscribe;
  }, [
    addLog,
    debugMode,
    isBackendReady,
    setRunning,
    setScreen,
    setSelectedJob,
    showToast,
    updateJobInHistory,
    updatePipelineProgress,
    updateStage,
  ]);

  // Global keyboard shortcuts (Matcha pattern)
  useInput((input, key) => {
    // Don't handle global shortcuts if confirm dialog is open
    if (confirmDialog) return;

    // '?' to toggle help overlay (repurposed from debug mode)
    if (input === '?') {
      setShowHelp(!showHelp);
      return;
    }

    // Ctrl+D to toggle debug mode (hidden shortcut)
    if (key.ctrl && input === 'd') {
      setDebugMode(!debugMode);
      return;
    }

    // Note: Navigation and 'q'/Escape are handled in individual screens
    // Pipeline can now run in background while user navigates
  });

  return (
    <Box flexDirection="column">
      {/* Global error banner with auto-dismiss (Matcha pattern) */}
      {globalError && (
        <Box marginBottom={1} paddingX={1}>
          <Text color="red" bold>Error: </Text>
          <Text color="red">{globalError}</Text>
        </Box>
      )}

      {/* Help overlay (takes priority over screen) */}
      {showHelp && <HelpOverlay />}

      {/* Confirmation dialog overlay */}
      {confirmDialog && <ConfirmDialog />}

      {/* Main screen content (hidden when overlays active) */}
      {!showHelp && !confirmDialog && <ScreenRouter screen={screen} />}

      {/* Toast notifications */}
      <ToastContainer />

      {/* Debug mode indicator (Matcha pattern) */}
      {debugMode && (
        <Box marginTop={1}>
          <Text color="yellow" bold>[DEBUG MODE]</Text>
          <Text dimColor> Screen: {screen}</Text>
        </Box>
      )}
    </Box>
  );
}

export default App;
