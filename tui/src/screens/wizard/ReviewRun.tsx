import React, { useState, useMemo } from 'react';
import { Box, Text, useInput } from 'ink';
import { Header } from '../../components/Header.js';
import { Footer } from '../../components/Footer.js';
import { useStore } from '../../store/index.js';
import { getBridge } from '../../services/python-bridge.js';
import { useTerminalSize } from '../../hooks/useTerminalSize.js';
import type { Screen, JobHistoryRecord } from '../../types/index.js';

const STAGE_NAMES: Record<string, string> = {
  mol_prep: 'Mol Prep',
  descriptors: 'Descriptors',
  struct_filters: 'Struct Filters',
  synthesis: 'Synthesis',
  docking: 'Docking',
  docking_filters: 'Docking Filters',
};

interface StageSummary {
  name: string;
  displayName: string;
  summary: string;
  dependency?: string;
}

export function ReviewRun(): React.ReactElement {
  const setScreen = useStore((state) => state.setScreen);
  const wizard = useStore((state) => state.wizard);
  const config = useStore((state) => state.configs.main);
  const setRunning = useStore((state) => state.setRunning);
  const addJobToHistory = useStore((state) => state.addJobToHistory);
  const updatePipelineProgress = useStore((state) => state.updatePipelineProgress);
  const resetWizard = useStore((state) => state.resetWizard);
  const showToast = useStore((state) => state.showToast);
  const getWizardSelectedStagesInOrder = useStore((state) => state.getWizardSelectedStagesInOrder);

  const { width: terminalWidth } = useTerminalSize();

  const [focusedIndex, setFocusedIndex] = useState(0);
  const [starting, setStarting] = useState(false);

  const selectedStagesInOrder = useMemo(() => getWizardSelectedStagesInOrder(), [wizard.stageOrder, wizard.selectedStages]);
  const totalSteps = 2 + selectedStagesInOrder.length + 1;

  // Build stage summaries
  const stageSummaries = useMemo((): StageSummary[] => {
    return selectedStagesInOrder.map((stageName) => {
      const params = wizard.stageConfigs[stageName]?.quickParams || {};
      const preset = wizard.stageConfigs[stageName]?.preset;
      let summary = '';
      let dependency: string | undefined;

      switch (stageName) {
        case 'mol_prep':
          summary = 'Datamol standardization';
          break;
        case 'descriptors':
          summary = `Batch: ${params.batch_size}`;
          if (preset) {
            summary += ` | Preset: ${preset}`;
          }
          break;
        case 'struct_filters':
          const nibr = params.calculate_NIBR ? 'Yes' : 'No';
          const lilly = params.calculate_lilly ? 'Yes' : 'No';
          summary = `NIBR: ${nibr} | Lilly: ${lilly}`;
          break;
        case 'synthesis':
          const saMax = params.sa_score_max === 'inf' ? '∞' : params.sa_score_max;
          summary = `SA: ${params.sa_score_min}-${saMax} | RA: ${params.ra_score_min}-${params.ra_score_max}`;
          break;
        case 'docking':
          summary = `Tool: ${params.tools} | Exhaust: ${params.exhaustiveness} | Modes: ${params.num_modes}`;
          break;
        case 'docking_filters':
          summary = `Mode: ${params.aggregation_mode} | Clashes≤${params.max_clashes} | H-bonds≥${params.min_hbonds} | RMSD≤${params.max_rmsd_to_conformer}`;
          break;
        default:
          summary = 'Configured';
      }

      return {
        name: stageName,
        displayName: STAGE_NAMES[stageName] || stageName,
        summary,
        dependency,
      };
    });
  }, [selectedStagesInOrder, wizard.stageConfigs, wizard.selectedStages]);

  const goBack = () => {
    // Always go back to stage selection from review
    setScreen('wizardStageSelection');
  };

  const getStageScreen = (stageName: string): Screen => {
    if (stageName === 'mol_prep') return 'wizardConfigMolPrep';
    if (stageName === 'descriptors') return 'wizardConfigDescriptors';
    if (stageName === 'struct_filters') return 'wizardConfigFilters';
    if (stageName === 'synthesis') return 'wizardConfigSynthesis';
    if (stageName === 'docking') return 'wizardConfigDocking';
    if (stageName === 'docking_filters') return 'wizardConfigDockingFilters';
    return 'wizardReview';
  };

  const startPipeline = async () => {
    if (starting) return;

    setStarting(true);
    try {
      const bridge = getBridge();
      const jobId = await bridge.startPipeline(selectedStagesInOrder);
      const now = new Date();

      const jobRecord: JobHistoryRecord = {
        id: jobId,
        name: `Pipeline ${jobId}`,
        startTime: now.toISOString(),
        status: 'running',
        config: {
          inputPath: config?.generated_mols_path || '',
          outputPath: config?.folder_to_save || '',
          stages: selectedStagesInOrder,
        },
      };

      addJobToHistory(jobRecord);

      updatePipelineProgress({
        currentStage: selectedStagesInOrder[0],
        stageIndex: 1,
        totalStages: selectedStagesInOrder.length,
        stageProgress: 0,
        latestMessage: 'Starting pipeline...',
      });

      setRunning(true, jobId);
      showToast('info', 'Pipeline started');
      setScreen('pipelineRunner');

      // Best-effort persistence to backend history (pipeline is already running).
      try {
        await bridge.addJob(
          jobId,
          null,
          config?.generated_mols_path || '',
          config?.folder_to_save || '',
          selectedStagesInOrder
        );
      } catch (err) {
        showToast('warning', `Failed to save job history: ${err}`);
      }
    } catch (err) {
      showToast('error', `Failed to start: ${err}`);
      setStarting(false);
    }
  };

  useInput((input, key) => {
    if (starting) return;

    if (key.upArrow) {
      setFocusedIndex(Math.max(0, focusedIndex - 1));
    } else if (key.downArrow) {
      setFocusedIndex(Math.min(stageSummaries.length - 1, focusedIndex + 1));
    } else if (input === 'e' && stageSummaries[focusedIndex]) {
      // Edit stage
      setScreen(getStageScreen(stageSummaries[focusedIndex].name));
    } else if (key.return) {
      startPipeline();
    } else if (key.escape || key.leftArrow || input === 'q') {
      goBack();
    }
  });

  const shortcuts = [
    { key: '↑↓', label: 'Navigate' },
    { key: 'e', label: 'Edit stage' },
    { key: 'Enter', label: 'Start' },
    { key: '←/Esc', label: 'Back' },
  ];

  return (
    <Box flexDirection="column" padding={1}>
      <Header title="Pipeline Wizard" subtitle={`Review & Run (${totalSteps}/${totalSteps})`} />

      {/* Input/Output info */}
      {config && (
        <Box flexDirection="column" marginY={1}>
          <Box>
            <Text dimColor>Input:  </Text>
            <Text color="cyan">{config.generated_mols_path || '(not set)'}</Text>
          </Box>
          <Box>
            <Text dimColor>Output: </Text>
            <Text color="cyan">{config.folder_to_save || '(not set)'}</Text>
          </Box>
        </Box>
      )}

      <Box marginY={1}>
        <Text color="cyan" bold>Pipeline Stages</Text>
      </Box>

      <Box flexDirection="column" marginY={1}>
        {stageSummaries.map((stage, index) => {
          const isFocused = focusedIndex === index;

          return (
            <Box key={stage.name} flexDirection="column" marginBottom={1}>
              <Box>
                <Text dimColor>{`${index + 1}. `}</Text>
                <Text color={isFocused ? 'cyan' : 'white'}>
                  {isFocused ? '> ' : '  '}
                </Text>
                <Text color={isFocused ? 'white' : 'gray'} bold>
                  {stage.displayName}
                </Text>
                {stage.dependency && (
                  <Text dimColor> ({stage.dependency})</Text>
                )}
              </Box>
              <Box paddingLeft={5}>
                <Text dimColor>{stage.summary}</Text>
              </Box>
            </Box>
          );
        })}
      </Box>

      <Box marginY={1}>
        <Text color="gray">{'─'.repeat(terminalWidth - 2)}</Text>
      </Box>

      {starting ? (
        <Box marginY={1}>
          <Text color="yellow">Starting pipeline...</Text>
        </Box>
      ) : (
        <Box marginY={1}>
          <Text color="green" bold>Press Enter to start the pipeline</Text>
        </Box>
      )}

      <Footer shortcuts={shortcuts} />
    </Box>
  );
}

export default ReviewRun;
