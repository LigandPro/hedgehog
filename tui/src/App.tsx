import React, { useEffect, useRef } from 'react';
import { Box, Text, useApp, useInput } from 'ink';
import { useStore } from './store/index.js';
import { getBridge, destroyBridge } from './services/python-bridge.js';
import { logger } from './utils/logger.js';
import { useTerminalSize } from './hooks/useTerminalSize.js';
import { useTheme } from './theme/context.js';
import { buildExitSummary } from './utils/exit-summary.js';
import type { Screen, JobHistoryRecord } from './types/index.js';

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
import { CommandMenu } from './components/CommandMenu.js';
import { ThemeMenu } from './components/ThemeMenu.js';

const STAGE_NAMES: Record<string, string> = {
  mol_prep: 'Mol Prep',
  descriptors: 'Descriptors',
  struct_filters: 'Struct Filters',
  synthesis: 'Synthesis',
  docking: 'Docking',
  docking_filters: 'Docking Filters',
};

const PROGRESS_LOG_THRESHOLDS = [10, 25, 50, 75, 90, 100] as const;
const REQUESTED_SESSION_ID = process.env.HEDGEHOG_TUI_SESSION?.trim() || null;

interface NotificationParams {
  stage?: string;
  current?: number;
  total?: number;
  message?: string;
  level?: string;
  results?: Record<string, unknown>;
}

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
  return <Box key={screen} flexDirection="column" flexGrow={1}>{content}</Box>;
}

export function App(): React.ReactElement {
  const { exit } = useApp();
  const { width, height } = useTerminalSize();
  const { theme } = useTheme();
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
  const showCommandMenu = useStore((state) => state.showCommandMenu);
  const setShowCommandMenu = useStore((state) => state.setShowCommandMenu);
  const setCommandMenuSeed = useStore((state) => state.setCommandMenuSeed);
  const showThemeMenu = useStore((state) => state.showThemeMenu);
  const setShowThemeMenu = useStore((state) => state.setShowThemeMenu);
  const confirmDialog = useStore((state) => state.confirmDialog);
  const hideConfirm = useStore((state) => state.hideConfirm);
  const setRunning = useStore((state) => state.setRunning);
  const updateStage = useStore((state) => state.updateStage);
  const updatePipelineProgress = useStore((state) => state.updatePipelineProgress);
  const addLog = useStore((state) => state.addLog);
  const updateJobInHistory = useStore((state) => state.updateJobInHistory);
  const showToast = useStore((state) => state.showToast);
  const setSelectedJob = useStore((state) => state.setSelectedJob);
  const setScreen = useStore((state) => state.setScreen);
  const currentJobId = useStore((state) => state.currentJobId);
  const selectedJobId = useStore((state) => state.selectedJobId);
  const jobHistory = useStore((state) => state.jobHistory);
  const setExitSummary = useStore((state) => state.setExitSummary);
  const inputLocked = useStore((state) => state.inputLocked);
  const searchActive = useStore((state) => state.searchActive);

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
          const typedHistory = history as unknown as JobHistoryRecord[];
          setJobHistory(typedHistory);
          logger.info(`Loaded ${typedHistory.length} jobs from history`);

          if (REQUESTED_SESSION_ID) {
            const requestedJob = typedHistory.find((job) => job.id === REQUESTED_SESSION_ID);
            if (requestedJob) {
              setSelectedJob(requestedJob.id);
              setScreen('results');
              logger.info(`Resumed session ${requestedJob.id}`);
            } else {
              logger.warn(`Requested session not found: ${REQUESTED_SESSION_ID}`);
            }
          }
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

    const handleProgress = (params: NotificationParams, selectedStages: string[]) => {
      const { stage, current, total, message } = params;
      if (!stage) return;
      const progress = (total ?? 0) > 0 ? Math.round(((current ?? 0) / (total ?? 1)) * 100) : 0;
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
    };

    const handleLog = (params: NotificationParams) => {
      if (!debugMode) return;
      addLog({
        timestamp: new Date(),
        level: (params.level as 'info' | 'warn' | 'error' | 'debug') || 'info',
        message: params.message ?? '',
      });
    };

    const handleStageStart = (params: NotificationParams, selectedStages: string[]) => {
      const stage = params.stage;
      if (!stage) return;
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
    };

    const handleStageComplete = (params: NotificationParams) => {
      const stage = params.stage;
      if (!stage) return;
      progressLogThresholdRef.current[stage] = 100;
      updateStage(stage, { status: 'completed', progress: 100 });
      addLog({
        timestamp: new Date(),
        level: 'info',
        message: `Stage completed: ${getStageName(stage)}`,
      });
    };

    const handleStageError = (params: NotificationParams) => {
      const stage = params.stage;
      if (!stage) return;
      const message = params.message || 'Unknown error';
      updateStage(stage, { status: 'error', message });
      addLog({
        timestamp: new Date(),
        level: 'error',
        message: `Stage failed: ${getStageName(stage)} - ${message}`,
      });
    };

    const handleComplete = (params: NotificationParams, state: ReturnType<typeof useStore.getState>) => {
      setRunning(false);
      const currentJobId = state.currentJobId;
      if (currentJobId) {
        setSelectedJob(currentJobId);
      }
      showToast('success', 'Pipeline completed');
      if (currentJobId) {
        const results = params.results || {};
        const moleculesProcessed = (results as Record<string, number>).molecules_processed || 0;
        updateJobInHistory(currentJobId, {
          status: 'completed',
          endTime: new Date().toISOString(),
          results: { moleculesProcessed },
        });
        bridge.updateJob(currentJobId, 'completed', { moleculesProcessed });
      }
      if (state.screen === 'pipelineRunner') {
        setScreen('results');
      }
      progressLogThresholdRef.current = {};
    };

    const handleError = (params: NotificationParams, state: ReturnType<typeof useStore.getState>) => {
      const message = params.message ?? 'Unknown error';
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
    };

    const unsubscribe = bridge.onNotification((method, rawParams) => {
      const state = useStore.getState();
      const stageOrder = state.wizard.stageOrder;
      const selectedStages = stageOrder.filter((s) => state.wizard.selectedStages.includes(s));
      const params = rawParams as NotificationParams;

      switch (method) {
        case 'progress':
          handleProgress(params, selectedStages);
          break;
        case 'log':
          handleLog(params);
          break;
        case 'stage_start':
          handleStageStart(params, selectedStages);
          break;
        case 'stage_complete':
          handleStageComplete(params);
          break;
        case 'stage_error':
          handleStageError(params);
          break;
        case 'complete':
          handleComplete(params, state);
          break;
        case 'error':
          handleError(params, state);
          break;
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
    if (key.ctrl && input === 'c') {
      if (confirmDialog) {
        hideConfirm();
        return;
      }
      if (showThemeMenu) {
        setShowThemeMenu(false);
        return;
      }
      if (showCommandMenu) {
        setShowCommandMenu(false);
        setCommandMenuSeed('/');
        return;
      }
      if (showHelp) {
        setShowHelp(false);
        return;
      }
      setExitSummary(buildExitSummary({
        currentJobId,
        selectedJobId,
        jobHistory,
        screen,
      }));
      exit();
      return;
    }

    if (key.ctrl && input === 'p') {
      if (showCommandMenu) {
        setShowCommandMenu(false);
        setCommandMenuSeed('/');
      } else if (!confirmDialog && !showThemeMenu && !showHelp && !inputLocked && !searchActive) {
        setCommandMenuSeed('/');
        setShowCommandMenu(true);
      }
      return;
    }

    // Don't handle global shortcuts if confirm dialog is open
    if (confirmDialog) return;
    if (showThemeMenu) return;
    if (showCommandMenu) return;
    if (showHelp) return;
    if (inputLocked || searchActive) return;

    if (input.startsWith('/')) {
      setCommandMenuSeed(input.toLowerCase());
      setShowCommandMenu(true);
      return;
    }

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

  const activeOverlay = confirmDialog
    ? 'confirm'
    : showThemeMenu
      ? 'theme'
      : showCommandMenu
        ? 'command'
        : showHelp
          ? 'help'
          : null;

  return (
    <Box
      flexDirection="column"
      width={width}
      height={height}
    >
      {/* Global error banner with auto-dismiss (Matcha pattern) */}
      {globalError && (
        <Box marginBottom={1} paddingX={1}>
          <Text color={theme.palette.error} bold>Error: </Text>
          <Text color={theme.palette.error}>{globalError}</Text>
        </Box>
      )}

      {/* Global overlays with deterministic precedence */}
      {activeOverlay === 'command' && <CommandMenu />}
      {activeOverlay === 'theme' && <ThemeMenu />}
      {activeOverlay === 'help' && <HelpOverlay />}
      {activeOverlay === 'confirm' && <ConfirmDialog />}

      {/* Main screen content (hidden when overlay is active) */}
      {activeOverlay === null && <ScreenRouter screen={screen} />}

      {/* Toast notifications */}
      <ToastContainer />

      {/* Debug mode indicator (Matcha pattern) */}
      {debugMode && (
        <Box marginTop={1}>
          <Text color={theme.palette.warning} bold>[DEBUG MODE]</Text>
          <Text color={theme.palette.textMuted}> Screen: {screen}</Text>
        </Box>
      )}
    </Box>
  );
}

export default App;
