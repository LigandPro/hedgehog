import { create } from 'zustand';
import type {
  Screen,
  MainConfig,
  DescriptorsConfig,
  FiltersConfig,
  SynthesisConfig,
  DockingConfig,
  StageInfo,
  LogEntry,
  StageStatus,
} from '../types/index.js';

interface Configs {
  main: MainConfig | null;
  descriptors: DescriptorsConfig | null;
  filters: FiltersConfig | null;
  synthesis: SynthesisConfig | null;
  docking: DockingConfig | null;
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

  // UI
  error: string | null;
  notification: string | null;
  isBackendReady: boolean;

  // Actions
  setScreen: (screen: Screen) => void;
  goBack: () => void;
  setConfig: <K extends keyof Configs>(type: K, data: Configs[K]) => void;
  updateConfig: <K extends keyof Configs>(type: K, data: Partial<NonNullable<Configs[K]>>) => void;
  markConfigDirty: (type: string, dirty: boolean) => void;
  setRunning: (running: boolean, jobId?: string) => void;
  updateStage: (stageName: string, update: Partial<StageInfo>) => void;
  addLog: (entry: LogEntry) => void;
  clearLogs: () => void;
  setError: (error: string | null) => void;
  setNotification: (notification: string | null) => void;
  setBackendReady: (ready: boolean) => void;
  reset: () => void;
}

const initialStages: Record<string, StageInfo> = {
  descriptors: { name: 'descriptors', displayName: 'Descriptors', status: 'pending', progress: 0 },
  struct_filters: { name: 'struct_filters', displayName: 'Struct Filters', status: 'pending', progress: 0 },
  synthesis: { name: 'synthesis', displayName: 'Synthesis', status: 'pending', progress: 0 },
  docking: { name: 'docking', displayName: 'Docking', status: 'pending', progress: 0 },
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
    docking: null,
  },
  configDirty: {},
  isRunning: false,
  currentJobId: null,
  stages: { ...initialStages },
  logs: [],
  error: null,
  notification: null,
  isBackendReady: false,

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

  setRunning: (running, jobId) => set({
    isRunning: running,
    currentJobId: jobId ?? null,
    stages: running ? { ...initialStages } : get().stages,
  }),

  updateStage: (stageName, update) => set((state) => ({
    stages: {
      ...state.stages,
      [stageName]: { ...state.stages[stageName], ...update },
    },
  })),

  addLog: (entry) => set((state) => ({
    logs: [...state.logs.slice(-99), entry],
  })),

  clearLogs: () => set({ logs: [] }),

  setError: (error) => set({ error }),

  setNotification: (notification) => set({ notification }),

  setBackendReady: (ready) => set({ isBackendReady: ready }),

  reset: () => set({
    isRunning: false,
    currentJobId: null,
    stages: { ...initialStages },
    logs: [],
    error: null,
  }),
}));

export default useStore;
