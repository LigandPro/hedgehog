import React, { memo, useMemo } from 'react';
import { Box, Text, useInput } from 'ink';
import { Header } from '../components/Header.js';
import { Footer } from '../components/Footer.js';
import { ProgressBar } from '../components/ProgressBar.js';
import { useStore } from '../store/index.js';
import { getBridge } from '../services/python-bridge.js';
import { formatDuration } from '../utils/format.js';
import { useTerminalSize } from '../hooks/useTerminalSize.js';
import type { LogEntry } from '../types/index.js';

const STAGE_NAMES: Record<string, string> = {
  mol_prep: 'Mol Prep',
  descriptors: 'Descriptors',
  struct_filters: 'Struct Filters',
  synthesis: 'Synthesis',
  docking: 'Docking',
  docking_filters: 'Docking Filters',
};

const MIN_LOG_PANEL_TERMINAL_WIDTH = 120;

function formatLogTimestamp(timestamp: Date): string {
  const value = new Date(timestamp);
  const hours = String(value.getHours()).padStart(2, '0');
  const minutes = String(value.getMinutes()).padStart(2, '0');
  const seconds = String(value.getSeconds()).padStart(2, '0');
  return `${hours}:${minutes}:${seconds}`;
}

function trimText(value: string, maxWidth: number): string {
  if (value.length <= maxWidth) {
    return value;
  }
  if (maxWidth <= 1) {
    return value.slice(0, maxWidth);
  }
  return `${value.slice(0, maxWidth - 1)}…`;
}

function getLevelColor(level: LogEntry['level']): 'gray' | 'cyan' | 'yellow' | 'red' | 'magenta' {
  if (level === 'warn') return 'yellow';
  if (level === 'error') return 'red';
  if (level === 'debug') return 'magenta';
  return 'cyan';
}

function getMessageColor(level: LogEntry['level']): 'gray' | 'yellow' | 'red' | 'magenta' {
  if (level === 'warn') return 'yellow';
  if (level === 'error') return 'red';
  if (level === 'debug') return 'magenta';
  return 'gray';
}

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

const LiveLogPanel = memo(function LiveLogPanel({
  panelWidth,
  terminalHeight,
}: {
  panelWidth: number;
  terminalHeight: number;
}): React.ReactElement {
  const logs = useStore((state) => state.logs);

  const visibleLineCount = useMemo(() => (
    Math.max(6, Math.min(24, terminalHeight - 20))
  ), [terminalHeight]);

  const messageWidth = Math.max(12, panelWidth - 24);
  const visibleLogs = logs.slice(-visibleLineCount);

  return (
    <Box flexDirection="column" marginLeft={2} width={panelWidth}>
      <Text color="cyan" bold>Live Log</Text>
      <Text color="gray">{'─'.repeat(Math.max(10, panelWidth - 1))}</Text>
      {visibleLogs.length === 0 ? (
        <Text dimColor>No log events yet.</Text>
      ) : (
        visibleLogs.map((entry, index) => {
          const timestamp = formatLogTimestamp(entry.timestamp);
          const level = entry.level.toUpperCase();
          const message = trimText(entry.message, messageWidth);

          return (
            <Text key={`${entry.timestamp.getTime()}-${index}`}>
              <Text color="gray">[{timestamp}] </Text>
              <Text color={getLevelColor(entry.level)}>[{level}]</Text>
              <Text color={getMessageColor(entry.level)}> {message}</Text>
            </Text>
          );
        })
      )}
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
  const showPipelineLog = useStore((state) => state.showPipelineLog);
  const togglePipelineLog = useStore((state) => state.togglePipelineLog);

  const { width: terminalWidth, height: terminalHeight } = useTerminalSize();

  const selectedStages = useMemo(() => (
    wizardStageOrder.filter((stage) => wizardSelectedStages.includes(stage))
  ), [wizardStageOrder, wizardSelectedStages]);

  const canShowRightLogPanel = showPipelineLog && terminalWidth >= MIN_LOG_PANEL_TERMINAL_WIDTH;
  const logPanelWidth = Math.max(36, Math.min(52, Math.floor(terminalWidth * 0.35)));

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
    if (input === 'l') {
      togglePipelineLog();
    } else if (key.escape || key.leftArrow) {
      // Navigate back - pipeline continues in background
      setScreen('welcome');
    } else if (input === 'c' && isRunning) {
      cancelPipeline();
    }
  });

  const shortcuts = isRunning
    ? [
        { key: 'c', label: 'Cancel' },
        { key: 'l', label: showPipelineLog ? 'Hide Log' : 'Show Log' },
        { key: '←/Esc', label: 'Back (keeps running)' },
      ]
    : [
        { key: 'l', label: showPipelineLog ? 'Hide Log' : 'Show Log' },
        { key: '←/Esc', label: 'Back' },
      ];

  return (
    <Box flexDirection="column" padding={1}>
      <Header title="Pipeline Runner" />

      <Box flexDirection="row">
        <Box flexDirection="column" flexGrow={1}>
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

          {showPipelineLog && !canShowRightLogPanel && (
            <Box marginTop={1}>
              <Text color="yellow">
                Live Log panel requires terminal width {'>='} {MIN_LOG_PANEL_TERMINAL_WIDTH} columns.
              </Text>
            </Box>
          )}
        </Box>

        {canShowRightLogPanel && (
          <LiveLogPanel panelWidth={logPanelWidth} terminalHeight={terminalHeight} />
        )}
      </Box>

      {/* Separator */}
      <Box marginTop={1}>
        <Text color="gray">{'─'.repeat(terminalWidth - 2)}</Text>
      </Box>

      <Footer shortcuts={shortcuts} />
    </Box>
  );
}

export default PipelineRunner;
