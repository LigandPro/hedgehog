import { getBridge } from '../../services/python-bridge.js';
import { useStore } from '../../store/index.js';
import type { JobHistoryRecord, PipelinePreflightResult, Screen } from '../../types/index.js';

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
    const preflight = await bridge.preflightPipeline(selectedStages);
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
    const jobId = await bridge.startPipeline(selectedStages);
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
