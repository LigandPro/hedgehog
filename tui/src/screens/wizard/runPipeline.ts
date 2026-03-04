import { getBridge } from '../../services/python-bridge.js';
import { useStore } from '../../store/index.js';
import type { JobHistoryRecord, PipelinePreflightResult, Screen } from '../../types/index.js';

type ConfigType = 'main' | 'mol_prep' | 'descriptors' | 'filters' | 'synthesis' | 'retrosynthesis' | 'docking' | 'docking_filters';

const STAGE_TO_CONFIG_TYPE: Record<string, Exclude<ConfigType, 'main'>> = {
  mol_prep: 'mol_prep',
  descriptors: 'descriptors',
  struct_filters: 'filters',
  synthesis: 'synthesis',
  docking: 'docking',
  docking_filters: 'docking_filters',
};

async function buildWizardConfigOverrides(
  selectedStages: string[]
): Promise<Partial<Record<ConfigType, Record<string, unknown>>>> {
  const bridge = getBridge();
  type StoreConfigs = ReturnType<typeof useStore.getState>['configs'];
  const overrides: Partial<Record<ConfigType, Record<string, unknown>>> = {};

  const getConfig = async (configType: ConfigType): Promise<Record<string, unknown> | null> => {
    const currentState = useStore.getState();
    const cached = currentState.configs[configType];
    if (cached) {
      return cached as unknown as Record<string, unknown>;
    }

    try {
      const loaded = await bridge.loadConfig(configType);
      useStore.getState().setConfig(configType as keyof StoreConfigs, loaded as never);
      return loaded as Record<string, unknown>;
    } catch (error) {
      const message = error instanceof Error ? error.message : String(error);
      useStore.getState().showToast(
        'warning',
        `Unable to load ${configType} config: ${message}`
      );
      return null;
    }
  };

  const main = await getConfig('main');
  if (main) {
    overrides.main = main;
  }

  for (const stage of selectedStages) {
    const configType = STAGE_TO_CONFIG_TYPE[stage];
    if (!configType) continue;
    const config = await getConfig(configType);
    if (config) {
      overrides[configType] = config;
    }
  }

  if (selectedStages.includes('synthesis')) {
    const retrosynthesis = await getConfig('retrosynthesis');
    if (retrosynthesis) {
      overrides.retrosynthesis = retrosynthesis;
    }
  }

  return overrides;
}

export interface PreflightCounters {
  errors: number;
  warnings: number;
  infos: number;
}

export interface RunWizardPipelineOptions {
  onPreflightErrorScreen?: Screen;
}

export interface RunWizardPipelineResult {
  started: boolean;
  preflight: PipelinePreflightResult | null;
  jobId?: string;
}

export function countPreflightChecks(preflight: PipelinePreflightResult | null): PreflightCounters {
  if (!preflight) {
    return { errors: 0, warnings: 0, infos: 0 };
  }

  const counters: PreflightCounters = { errors: 0, warnings: 0, infos: 0 };
  const allChecks = [
    ...preflight.checks,
    ...preflight.stage_reports.flatMap((report) => report.checks),
  ];

  for (const check of allChecks) {
    if (check.level === 'error') counters.errors += 1;
    else if (check.level === 'warning') counters.warnings += 1;
    else counters.infos += 1;
  }

  return counters;
}

export async function refreshWizardPreflight(
  stagesOverride?: string[]
): Promise<PipelinePreflightResult | null> {
  const state = useStore.getState();
  const selectedStages = stagesOverride ?? state.getWizardSelectedStagesInOrder();

  if (selectedStages.length === 0) {
    state.setWizardPreflight(null);
    return null;
  }

  try {
    const bridge = getBridge();
    const configOverrides = await buildWizardConfigOverrides(selectedStages);
    const preflight = await bridge.preflightPipeline(selectedStages, configOverrides);
    useStore.getState().setWizardPreflight(preflight);
    return preflight;
  } catch (error) {
    const message = error instanceof Error ? error.message : String(error);
    useStore.getState().showToast('error', `Preflight failed: ${message}`);
    return null;
  }
}

export async function runWizardPipeline(
  options: RunWizardPipelineOptions = {}
): Promise<RunWizardPipelineResult> {
  const state = useStore.getState();
  const selectedStages = state.getWizardSelectedStagesInOrder();

  if (selectedStages.length === 0) {
    state.showToast('warning', 'Select at least one stage');
    return { started: false, preflight: null };
  }

  const preflight = await refreshWizardPreflight(selectedStages);
  if (!preflight) {
    return { started: false, preflight: null };
  }

  const counters = countPreflightChecks(preflight);
  if (counters.errors > 0) {
    state.showToast('error', `Preflight has ${counters.errors} blocking error(s)`);
    if (options.onPreflightErrorScreen) {
      state.setScreen(options.onPreflightErrorScreen);
    }
    return { started: false, preflight };
  }

  if (counters.warnings > 0) {
    state.showToast('warning', `Preflight warnings: ${counters.warnings}`);
  }

  try {
    const bridge = getBridge();
    const latestState = useStore.getState();
    const configOverrides = await buildWizardConfigOverrides(selectedStages);
    const jobId = await bridge.startPipeline(selectedStages, configOverrides);
    const now = new Date();
    const mainConfig = latestState.configs.main;

    const jobRecord: JobHistoryRecord = {
      id: jobId,
      name: `Pipeline ${jobId}`,
      startTime: now.toISOString(),
      status: 'running',
      config: {
        inputPath: mainConfig?.generated_mols_path || '',
        outputPath: mainConfig?.folder_to_save || '',
        stages: selectedStages,
      },
    };

    latestState.addJobToHistory(jobRecord);
    latestState.updatePipelineProgress({
      currentStage: selectedStages[0],
      stageIndex: 1,
      totalStages: selectedStages.length,
      stageProgress: 0,
      latestMessage: 'Starting pipeline...',
    });

    latestState.setRunning(true, jobId);
    latestState.showToast('info', 'Pipeline started');
    latestState.setScreen('pipelineRunner');

    try {
      await bridge.addJob(
        jobId,
        null,
        mainConfig?.generated_mols_path || '',
        mainConfig?.folder_to_save || '',
        selectedStages
      );
    } catch (error) {
      const message = error instanceof Error ? error.message : String(error);
      useStore.getState().showToast('warning', `Failed to save job history: ${message}`);
    }

    return { started: true, preflight, jobId };
  } catch (error) {
    const message = error instanceof Error ? error.message : String(error);
    useStore.getState().showToast('error', `Failed to start: ${message}`);
    return { started: false, preflight };
  }
}
