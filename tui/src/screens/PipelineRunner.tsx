import React, { memo, useMemo } from 'react';
import { Box, Text, useInput } from 'ink';
import { Header } from '../components/Header.js';
import { Footer } from '../components/Footer.js';
import { ProgressBar } from '../components/ProgressBar.js';
import { useStore } from '../store/index.js';
import { getBridge } from '../services/python-bridge.js';
import { formatDuration } from '../utils/format.js';
import { useTerminalSize } from '../hooks/useTerminalSize.js';

const STAGE_NAMES: Record<string, string> = {
  descriptors: 'Descriptors',
  struct_filters: 'Struct Filters',
  synthesis: 'Synthesis',
  docking: 'Docking',
  docking_filters: 'Docking Filters',
};

// Overall progress display with current stage bar
const PipelineProgress = memo(function PipelineProgress({
  selectedStages,
}: {
  selectedStages: string[];
}): React.ReactElement {
  const stages = useStore((state) => state.stages);
  const pipelineProgress = useStore((state) => state.pipelineProgress);

  // Calculate completed stages
  const completedCount = selectedStages.filter(
    (stage) => stages[stage]?.status === 'completed'
  ).length;

  // Find current running stage
  const currentStage = selectedStages.find(
    (stage) => stages[stage]?.status === 'running'
  );
  const currentStageInfo = currentStage ? stages[currentStage] : null;

  // Overall progress percentage
  const overallProgress = Math.round(
    ((completedCount + (currentStageInfo?.progress || 0) / 100) / selectedStages.length) * 100
  );

  return (
    <Box flexDirection="column" marginY={1}>
      {/* Overall progress */}
      <Box marginBottom={1}>
        <Text color="cyan" bold>Overall: </Text>
        <Text color="cyan">{completedCount}/{selectedStages.length} stages</Text>
        <Text dimColor> ({overallProgress}%)</Text>
      </Box>

      {/* Overall progress bar */}
      <ProgressBar
        progress={overallProgress}
        status={completedCount === selectedStages.length ? 'completed' : 'running'}
        width={40}
        showPercentage={false}
      />

      {/* Stages list with status icons */}
      <Box flexDirection="column" marginTop={1}>
        {selectedStages.map((stage, index) => {
          const stageInfo = stages[stage];
          const status = stageInfo?.status || 'pending';
          const icon = status === 'completed' ? '✓' : status === 'running' ? '●' : status === 'error' ? '✗' : '○';
          const color = status === 'completed' ? 'green' : status === 'running' ? 'cyan' : status === 'error' ? 'red' : 'gray';

          return (
            <Box key={stage}>
              <Text dimColor>{`${index + 1}. `}</Text>
              <Text color={color}>{icon} </Text>
              <Text color={status === 'running' ? 'white' : 'gray'}>
                {STAGE_NAMES[stage] || stage}
              </Text>
              {status === 'running' && stageInfo?.progress !== undefined && (
                <Text color="cyan"> ({stageInfo.progress}%)</Text>
              )}
            </Box>
          );
        })}
      </Box>

      {/* Current stage progress bar removed to avoid duplicate progress display */}
    </Box>
  );
});

// Running status display
const RunningStatus = memo(function RunningStatus(): React.ReactElement | null {
  const isRunning = useStore((state) => state.isRunning);
  const elapsedSeconds = useStore((state) => state.elapsedSeconds);

  if (!isRunning) return null;

  return (
    <Box marginY={1}>
      <Text color="cyan">Elapsed: </Text>
      <Text color="cyan" bold>
        {formatDuration(elapsedSeconds)}
      </Text>
    </Box>
  );
});

export function PipelineRunner(): React.ReactElement {
  const setScreen = useStore((state) => state.setScreen);
  const isRunning = useStore((state) => state.isRunning);
  const currentJobId = useStore((state) => state.currentJobId);
  const setRunning = useStore((state) => state.setRunning);
  const addLog = useStore((state) => state.addLog);
  const isBackendReady = useStore((state) => state.isBackendReady);
  const config = useStore((state) => state.configs.main);
  const updateJobInHistory = useStore((state) => state.updateJobInHistory);
  const showConfirm = useStore((state) => state.showConfirm);
  const showToast = useStore((state) => state.showToast);
  const wizardStageOrder = useStore((state) => state.wizard.stageOrder);
  const wizardSelectedStages = useStore((state) => state.wizard.selectedStages);

  const { width: terminalWidth } = useTerminalSize();

  const selectedStages = useMemo(() => (
    wizardStageOrder.filter((stage) => wizardSelectedStages.includes(stage))
  ), [wizardStageOrder, wizardSelectedStages]);

  const doCancelPipeline = async () => {
    if (!currentJobId) return;
    try {
      const bridge = getBridge();
      await bridge.cancelPipeline(currentJobId);
      setRunning(false);

      updateJobInHistory(currentJobId, {
        status: 'cancelled',
        endTime: new Date().toISOString(),
      });
      await bridge.updateJob(currentJobId, 'cancelled');

      addLog({
        timestamp: new Date(),
        level: 'warn',
        message: 'Pipeline cancelled by user',
      });
      showToast('warning', 'Pipeline cancelled');
    } catch (err) {
      showToast('error', `Failed to cancel: ${err}`);
    }
  };

  const cancelPipeline = () => {
    showConfirm({
      title: 'Cancel Pipeline?',
      message: 'Stop the running pipeline? Progress will be lost.',
      confirmLabel: 'Cancel Pipeline',
      cancelLabel: 'Keep Running',
      onConfirm: doCancelPipeline,
    });
  };

  useInput((input, key) => {
    if (key.escape || key.leftArrow) {
      // Navigate back - pipeline continues in background
      setScreen('welcome');
    } else if (input === 'c' && isRunning) {
      cancelPipeline();
    }
  });

  const shortcuts = isRunning
    ? [
        { key: 'c', label: 'Cancel' },
        { key: '←/Esc', label: 'Back (keeps running)' },
      ]
    : [{ key: '←/Esc', label: 'Back' }];

  return (
    <Box flexDirection="column" padding={1}>
      <Header title="Pipeline Runner" />

      {!isBackendReady && (
        <Box marginY={1}>
          <Text color="yellow">Warning: Backend not connected</Text>
        </Box>
      )}

      {/* Input/Output info */}
      {config && (
        <Box flexDirection="column" marginY={1}>
          <Box>
            <Text dimColor>Input:  </Text>
            <Text color="cyan">{config.generated_mols_path}</Text>
          </Box>
          <Box>
            <Text dimColor>Output: </Text>
            <Text color="cyan">{config.folder_to_save}</Text>
          </Box>
        </Box>
      )}

      {/* Running status */}
      {isRunning ? (
        <>
          <RunningStatus />
          <PipelineProgress selectedStages={selectedStages} />
        </>
      ) : (
        <Box marginY={2}>
          <Text dimColor>Pipeline not running. Use [n] from Welcome screen to start a new run.</Text>
        </Box>
      )}

      {/* Separator */}
      <Box marginTop={1}>
        <Text color="gray">{'─'.repeat(terminalWidth - 2)}</Text>
      </Box>

      <Footer shortcuts={shortcuts} />
    </Box>
  );
}

export default PipelineRunner;
