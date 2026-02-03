import { create } from 'zustand';
import type {
  Screen,
  MainConfig,
  DescriptorsConfig,
  FiltersConfig,
  SynthesisConfig,
  RetrosynthesisConfig,
  DockingConfig,
  StageInfo,
  LogEntry,
  StageStatus,
  JobHistoryRecord,
  Toast,
  ToastType,
  ConfirmDialogConfig,
  WizardState,
  WizardStep,
  WizardStageConfig,
} from '../types/index.js';

interface Configs {
  main: MainConfig | null;
  descriptors: DescriptorsConfig | null;
  filters: FiltersConfig | null;
  synthesis: SynthesisConfig | null;
  retrosynthesis: RetrosynthesisConfig | null;
  docking: DockingConfig | null;
}

// Global pipeline progress for header display
interface PipelineProgress {
  currentStage: string;
  stageIndex: number;
  totalStages: number;
  stageProgress: number;
  latestMessage: string;
}

interface AppState {
  // Navigation
  screen: Screen;
  previousScreen: Screen | null;

  // Configs
  configs: Configs;
  configDirty: Record<string, boolean>;

  // Pipeline state
  isRunning: boolean;
  currentJobId: string | null;
  stages: Record<string, StageInfo>;
  logs: LogEntry[];

  // Global pipeline progress (for Header)
  pipelineProgress: PipelineProgress;
  elapsedSeconds: number;
  _elapsedTimer: ReturnType<typeof setInterval> | null;
  _startTime: number | null;

  // UI
  error: string | null;
  notification: string | null;
  isBackendReady: boolean;

  // Job history (Matcha pattern)
  jobHistory: JobHistoryRecord[];
  selectedJobId: string | null;

  // Debug mode (Matcha pattern)
  debugMode: boolean;

  // Global error with auto-dismiss (Matcha pattern)
  globalError: string | null;
  globalErrorTimer: ReturnType<typeof setTimeout> | null;

  // Toast notifications
  toasts: Toast[];

  // Confirmation dialog
  confirmDialog: ConfirmDialogConfig | null;

  // Help overlay
  showHelp: boolean;

  // Global search state
  searchQuery: string;
  searchActive: boolean;

  // Actions
  setScreen: (screen: Screen) => void;
  goBack: () => void;
  setConfig: <K extends keyof Configs>(type: K, data: Configs[K]) => void;
  updateConfig: <K extends keyof Configs>(type: K, data: Partial<NonNullable<Configs[K]>>) => void;
  markConfigDirty: (type: string, dirty: boolean) => void;
  setRunning: (running: boolean, jobId?: string) => void;
  updateStage: (stageName: string, update: Partial<StageInfo>) => void;
  updatePipelineProgress: (progress: Partial<PipelineProgress>) => void;
  addLog: (entry: LogEntry) => void;
  clearLogs: () => void;
  setError: (error: string | null) => void;
  setNotification: (notification: string | null) => void;
  setBackendReady: (ready: boolean) => void;
  reset: () => void;

  // Elapsed timer actions
  startElapsedTimer: () => void;
  stopElapsedTimer: () => void;

  // Job history actions (Matcha pattern)
  setJobHistory: (history: JobHistoryRecord[]) => void;
  addJobToHistory: (job: JobHistoryRecord) => void;
  updateJobInHistory: (jobId: string, update: Partial<JobHistoryRecord>) => void;
  removeJobFromHistory: (jobId: string) => void;
  setSelectedJob: (jobId: string | null) => void;

  // Debug mode actions
  setDebugMode: (enabled: boolean) => void;

  // Global error actions
  setGlobalError: (error: string | null, autoDismiss?: boolean) => void;

  // Toast actions
  showToast: (type: ToastType, message: string, duration?: number) => void;
  hideToast: (id: string) => void;
  clearToasts: () => void;

  // Confirm dialog actions
  showConfirm: (config: ConfirmDialogConfig) => void;
  hideConfirm: () => void;

  // Help overlay actions
  setShowHelp: (show: boolean) => void;

  // Search actions
  setSearchQuery: (query: string) => void;
  setSearchActive: (active: boolean) => void;
  clearSearch: () => void;

  // Wizard state
  wizard: WizardState;

  // Wizard actions
  setWizardStep: (step: WizardStep) => void;
  toggleWizardStage: (stage: string) => void;
  setWizardStageEnabled: (stage: string, enabled: boolean) => void;
  reorderWizardStage: (stage: string, direction: 'up' | 'down') => void;
  setWizardStageOrder: (order: string[]) => void;
  setWizardQuickParam: (stage: string, key: string, value: unknown) => void;
  applyWizardPreset: (stage: string, presetName: string, presetValues: Record<string, unknown>) => void;
  setWizardDependency: (key: keyof WizardState['dependencies'], value: boolean) => void;
  resetWizard: () => void;
  getWizardSelectedStagesInOrder: () => string[];
}

const initialStages: Record<string, StageInfo> = {
  descriptors: { name: 'descriptors', displayName: 'Descriptors', status: 'pending', progress: 0 },
  struct_filters: { name: 'struct_filters', displayName: 'Struct Filters', status: 'pending', progress: 0 },
  synthesis: { name: 'synthesis', displayName: 'Synthesis', status: 'pending', progress: 0 },
  docking: { name: 'docking', displayName: 'Docking', status: 'pending', progress: 0 },
};

const initialPipelineProgress: PipelineProgress = {
  currentStage: '',
  stageIndex: 0,
  totalStages: 0,
  stageProgress: 0,
  latestMessage: '',
};

const WIZARD_STAGES = ['descriptors', 'struct_filters', 'synthesis', 'docking'];

const initialWizardState: WizardState = {
  currentStep: 'stage-selection',
  selectedStages: [...WIZARD_STAGES],
  stageOrder: [...WIZARD_STAGES],
  stageConfigs: {
    descriptors: {
      enabled: true,
      order: 0,
      quickParams: {
        batch_size: 1000,
        filter_data: true,
        molWt_min: 200,
        molWt_max: 500,
        logP_min: -0.4,
        logP_max: 5.6,
        hbd_min: 0,
        hbd_max: 5,
        hba_min: 0,
        hba_max: 10,
        tpsa_min: 0,
        tpsa_max: 140,
      },
      preset: 'Drug-like',
    },
    struct_filters: {
      enabled: true,
      order: 1,
      quickParams: {
        run_before_descriptors: true,
        calculate_common_alerts: true,
        calculate_NIBR: true,
        calculate_lilly: true,
        filter_data: true,
      },
    },
    synthesis: {
      enabled: true,
      order: 2,
      quickParams: {
        filter_solved_only: true,
        sa_score_min: 0,
        sa_score_max: 4.5,
        syba_score_min: 0,
        syba_score_max: 'inf',
        ra_score_min: 0.5,
        ra_score_max: 1,
      },
    },
    docking: {
      enabled: true,
      order: 3,
      quickParams: {
        tools: 'smina',
        exhaustiveness: 8,
        num_modes: 9,
      },
    },
  },
  dependencies: {
    runFiltersBeforeDescriptors: true,
  },
};

export const useStore = create<AppState>((set, get) => ({
  // Initial state
  screen: 'welcome',
  previousScreen: null,
  configs: {
    main: null,
    descriptors: null,
    filters: null,
    synthesis: null,
    retrosynthesis: null,
    docking: null,
  },
  configDirty: {},
  isRunning: false,
  currentJobId: null,
  stages: { ...initialStages },
  logs: [],

  // Global pipeline progress
  pipelineProgress: { ...initialPipelineProgress },
  elapsedSeconds: 0,
  _elapsedTimer: null,
  _startTime: null,

  error: null,
  notification: null,
  isBackendReady: false,

  // Job history (Matcha pattern)
  jobHistory: [],
  selectedJobId: null,

  // Debug mode (Matcha pattern)
  debugMode: false,

  // Global error (Matcha pattern)
  globalError: null,
  globalErrorTimer: null,

  // Toast notifications
  toasts: [],

  // Confirmation dialog
  confirmDialog: null,

  // Help overlay
  showHelp: false,

  // Global search state
  searchQuery: '',
  searchActive: false,

  // Wizard state
  wizard: { ...initialWizardState },

  // Actions
  setScreen: (screen) => set((state) => ({
    screen,
    previousScreen: state.screen,
  })),

  goBack: () => set((state) => ({
    screen: state.previousScreen || 'welcome',
    previousScreen: null,
  })),

  setConfig: (type, data) => set((state) => ({
    configs: { ...state.configs, [type]: data },
    configDirty: { ...state.configDirty, [type]: false },
  })),

  updateConfig: (type, data) => set((state) => ({
    configs: {
      ...state.configs,
      [type]: state.configs[type] ? { ...state.configs[type], ...data } : data,
    },
    configDirty: { ...state.configDirty, [type]: true },
  })),

  markConfigDirty: (type, dirty) => set((state) => ({
    configDirty: { ...state.configDirty, [type]: dirty },
  })),

  setRunning: (running, jobId) => {
    if (running) {
      get().startElapsedTimer();
    } else {
      get().stopElapsedTimer();
    }
    set({
      isRunning: running,
      currentJobId: jobId ?? null,
      stages: running ? { ...initialStages } : get().stages,
      pipelineProgress: running ? { ...initialPipelineProgress } : get().pipelineProgress,
    });
  },

  updateStage: (stageName, update) => set((state) => {
    const current = state.stages[stageName];
    if (!current) return state;

    let changed = false;
    for (const key of Object.keys(update) as (keyof StageInfo)[]) {
      if (current[key] !== update[key]) {
        changed = true;
        break;
      }
    }

    if (!changed) return state;

    return {
      stages: {
        ...state.stages,
        [stageName]: { ...current, ...update },
      },
    };
  }),

  updatePipelineProgress: (progress) => set((state) => {
    let changed = false;
    for (const key of Object.keys(progress) as (keyof PipelineProgress)[]) {
      if (state.pipelineProgress[key] !== progress[key]) {
        changed = true;
        break;
      }
    }

    if (!changed) return state;

    return {
      pipelineProgress: { ...state.pipelineProgress, ...progress },
    };
  }),

  addLog: (entry) => set((state) => ({
    logs: [...state.logs.slice(-99), entry],
    pipelineProgress: {
      ...state.pipelineProgress,
      latestMessage: entry.message,
    },
  })),

  clearLogs: () => set({ logs: [] }),

  setError: (error) => set({ error }),

  setNotification: (notification) => set({ notification }),

  setBackendReady: (ready) => set({ isBackendReady: ready }),

  reset: () => {
    get().stopElapsedTimer();
    set({
      isRunning: false,
      currentJobId: null,
      stages: { ...initialStages },
      logs: [],
      error: null,
      pipelineProgress: { ...initialPipelineProgress },
      elapsedSeconds: 0,
    });
  },

  // Elapsed timer actions
  startElapsedTimer: () => {
    const state = get();
    if (state._elapsedTimer) {
      clearInterval(state._elapsedTimer);
    }
    const startTime = Date.now();
    const timer = setInterval(() => {
      set({ elapsedSeconds: Math.floor((Date.now() - startTime) / 1000) });
    }, 1000);
    set({ _elapsedTimer: timer, _startTime: startTime, elapsedSeconds: 0 });
  },

  stopElapsedTimer: () => {
    const state = get();
    if (state._elapsedTimer) {
      clearInterval(state._elapsedTimer);
      set({ _elapsedTimer: null, _startTime: null });
    }
  },

  // Job history actions (Matcha pattern)
  setJobHistory: (history) => set({ jobHistory: history }),

  addJobToHistory: (job) => set((state) => ({
    jobHistory: [job, ...state.jobHistory].slice(0, 50),
  })),

  updateJobInHistory: (jobId, update) => set((state) => ({
    jobHistory: state.jobHistory.map((j) =>
      j.id === jobId ? { ...j, ...update } : j
    ),
  })),

  removeJobFromHistory: (jobId) => set((state) => ({
    jobHistory: state.jobHistory.filter((j) => j.id !== jobId),
    selectedJobId: state.selectedJobId === jobId ? null : state.selectedJobId,
  })),

  setSelectedJob: (jobId) => set({ selectedJobId: jobId }),

  // Debug mode actions
  setDebugMode: (enabled) => set({ debugMode: enabled }),

  // Global error actions (Matcha pattern - auto-dismiss after 5s)
  setGlobalError: (error, autoDismiss = true) => {
    const state = get();
    if (state.globalErrorTimer) {
      clearTimeout(state.globalErrorTimer);
    }

    if (error && autoDismiss) {
      const timer = setTimeout(() => {
        set({ globalError: null, globalErrorTimer: null });
      }, 5000);
      set({ globalError: error, globalErrorTimer: timer });
    } else {
      set({ globalError: error, globalErrorTimer: null });
    }
  },

  // Toast actions
  showToast: (type, message, duration) => {
    const defaultDurations: Record<ToastType, number> = {
      success: 2000,
      error: 5000,
      info: 3000,
      warning: 4000,
    };
    const actualDuration = duration ?? defaultDurations[type];
    const id = `toast-${Date.now()}-${Math.random().toString(36).slice(2, 9)}`;
    const toast: Toast = { id, type, message, duration: actualDuration };

    set((state) => ({
      toasts: [...state.toasts.slice(-2), toast], // Keep max 3 toasts
    }));

    // Auto-dismiss
    setTimeout(() => {
      set((state) => ({
        toasts: state.toasts.filter((t) => t.id !== id),
      }));
    }, actualDuration);
  },

  hideToast: (id) => set((state) => ({
    toasts: state.toasts.filter((t) => t.id !== id),
  })),

  clearToasts: () => set({ toasts: [] }),

  // Confirm dialog actions
  showConfirm: (config) => set({ confirmDialog: config }),

  hideConfirm: () => {
    const state = get();
    if (state.confirmDialog?.onCancel) {
      state.confirmDialog.onCancel();
    }
    set({ confirmDialog: null });
  },

  // Help overlay actions
  setShowHelp: (show) => set({ showHelp: show }),

  // Search actions
  setSearchQuery: (query) => set({ searchQuery: query }),

  setSearchActive: (active) => set({
    searchActive: active,
    searchQuery: active ? get().searchQuery : '',
  }),

  clearSearch: () => set({ searchQuery: '', searchActive: false }),

  // Wizard actions
  setWizardStep: (step) => set((state) => ({
    wizard: { ...state.wizard, currentStep: step },
  })),

  toggleWizardStage: (stage) => set((state) => {
    const selected = state.wizard.selectedStages;
    const newSelected = selected.includes(stage)
      ? selected.filter((s) => s !== stage)
      : [...selected, stage];

    // Update stageConfigs enabled status
    const newConfigs = { ...state.wizard.stageConfigs };
    if (newConfigs[stage]) {
      newConfigs[stage] = { ...newConfigs[stage], enabled: newSelected.includes(stage) };
    }

    return {
      wizard: {
        ...state.wizard,
        selectedStages: newSelected,
        stageConfigs: newConfigs,
      },
    };
  }),

  setWizardStageEnabled: (stage, enabled) => set((state) => {
    const selected = state.wizard.selectedStages;
    const newSelected = enabled
      ? (selected.includes(stage) ? selected : [...selected, stage])
      : selected.filter((s) => s !== stage);

    const newConfigs = { ...state.wizard.stageConfigs };
    if (newConfigs[stage]) {
      newConfigs[stage] = { ...newConfigs[stage], enabled };
    }

    return {
      wizard: {
        ...state.wizard,
        selectedStages: newSelected,
        stageConfigs: newConfigs,
      },
    };
  }),

  reorderWizardStage: (stage, direction) => set((state) => {
    const order = [...state.wizard.stageOrder];
    const index = order.indexOf(stage);
    if (index === -1) return state;

    const newIndex = direction === 'up' ? index - 1 : index + 1;
    if (newIndex < 0 || newIndex >= order.length) return state;

    // Swap
    [order[index], order[newIndex]] = [order[newIndex], order[index]];

    // Update order in stageConfigs
    const newConfigs = { ...state.wizard.stageConfigs };
    order.forEach((s, i) => {
      if (newConfigs[s]) {
        newConfigs[s] = { ...newConfigs[s], order: i };
      }
    });

    return {
      wizard: {
        ...state.wizard,
        stageOrder: order,
        stageConfigs: newConfigs,
      },
    };
  }),

  setWizardStageOrder: (order) => set((state) => {
    const newConfigs = { ...state.wizard.stageConfigs };
    order.forEach((s, i) => {
      if (newConfigs[s]) {
        newConfigs[s] = { ...newConfigs[s], order: i };
      }
    });

    return {
      wizard: {
        ...state.wizard,
        stageOrder: order,
        stageConfigs: newConfigs,
      },
    };
  }),

  setWizardQuickParam: (stage, key, value) => set((state) => {
    const stageConfig = state.wizard.stageConfigs[stage];
    if (!stageConfig) return state;

    return {
      wizard: {
        ...state.wizard,
        stageConfigs: {
          ...state.wizard.stageConfigs,
          [stage]: {
            ...stageConfig,
            quickParams: {
              ...stageConfig.quickParams,
              [key]: value,
            },
          },
        },
      },
    };
  }),

  applyWizardPreset: (stage, presetName, presetValues) => set((state) => {
    const stageConfig = state.wizard.stageConfigs[stage];
    if (!stageConfig) return state;

    return {
      wizard: {
        ...state.wizard,
        stageConfigs: {
          ...state.wizard.stageConfigs,
          [stage]: {
            ...stageConfig,
            quickParams: {
              ...stageConfig.quickParams,
              ...presetValues,
            },
            preset: presetName,
          },
        },
      },
    };
  }),

  setWizardDependency: (key, value) => set((state) => ({
    wizard: {
      ...state.wizard,
      dependencies: {
        ...state.wizard.dependencies,
        [key]: value,
      },
    },
  })),

  resetWizard: () => set({
    wizard: { ...initialWizardState },
  }),

  getWizardSelectedStagesInOrder: () => {
    const state = get();
    const { selectedStages, stageOrder } = state.wizard;
    return stageOrder.filter((s) => selectedStages.includes(s));
  },
}));

export default useStore;
