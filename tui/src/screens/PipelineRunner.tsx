import React, { useState, useEffect } from 'react';
import { Box, Text, useInput } from 'ink';
import { Header } from '../components/Header.js';
import { Footer } from '../components/Footer.js';
import { Spinner } from '../components/Spinner.js';
import { ProgressBar } from '../components/ProgressBar.js';
import { useStore } from '../store/index.js';
import { getBridge } from '../services/python-bridge.js';
import { formatTimestamp } from '../utils/format.js';
import type { StageStatus, LogEntry } from '../types/index.js';

const STAGES = ['descriptors', 'struct_filters', 'synthesis', 'docking'];
const STAGE_NAMES: Record<string, string> = {
  descriptors: 'Descriptors',
  struct_filters: 'Struct Filters',
  synthesis: 'Synthesis',
  docking: 'Docking',
};

export function PipelineRunner(): React.ReactElement {
  const setScreen = useStore((state) => state.setScreen);
  const stages = useStore((state) => state.stages);
  const logs = useStore((state) => state.logs);
  const isRunning = useStore((state) => state.isRunning);
  const currentJobId = useStore((state) => state.currentJobId);
  const setRunning = useStore((state) => state.setRunning);
  const updateStage = useStore((state) => state.updateStage);
  const addLog = useStore((state) => state.addLog);
  const isBackendReady = useStore((state) => state.isBackendReady);
  const config = useStore((state) => state.configs.main);

  const [error, setError] = useState<string | null>(null);
  const [selectedStages, setSelectedStages] = useState<string[]>(STAGES);
  const [focusedIndex, setFocusedIndex] = useState(0);

  // Setup notification handler
  useEffect(() => {
    const bridge = getBridge();
    const unsubscribe = bridge.onNotification((method, params) => {
      if (method === 'progress') {
        const { stage, current, total, message } = params as any;
        const progress = total > 0 ? Math.round((current / total) * 100) : 0;
        updateStage(stage, {
          progress,
          status: progress === 100 ? 'completed' : 'running',
          message,
        });
      } else if (method === 'log') {
        addLog({
          timestamp: new Date(),
          level: (params as any).level || 'info',
          message: (params as any).message,
        });
      } else if (method === 'stage_start') {
        updateStage((params as any).stage, { status: 'running', progress: 0 });
      } else if (method === 'stage_complete') {
        updateStage((params as any).stage, { status: 'completed', progress: 100 });
      } else if (method === 'stage_error') {
        updateStage((params as any).stage, { status: 'error', message: (params as any).message });
      } else if (method === 'complete') {
        setRunning(false);
        addLog({
          timestamp: new Date(),
          level: 'info',
          message: 'Pipeline completed successfully!',
        });
      } else if (method === 'error') {
        setError((params as any).message);
        setRunning(false);
      }
    });

    return unsubscribe;
  }, []);

  const toggleStage = (stage: string) => {
    setSelectedStages(prev => 
      prev.includes(stage) 
        ? prev.filter(s => s !== stage)
        : [...prev, stage]
    );
  };

  const startPipeline = async () => {
    if (selectedStages.length === 0) {
      setError('Please select at least one stage');
      return;
    }
    try {
      setError(null);
      const bridge = getBridge();
      const jobId = await bridge.startPipeline(selectedStages);
      setRunning(true, jobId);
      addLog({
        timestamp: new Date(),
        level: 'info',
        message: \,
      });
    } catch (err) {
      setError(String(err));
    }
  };

  const cancelPipeline = async () => {
    if (!currentJobId) return;
    try {
      const bridge = getBridge();
      await bridge.cancelPipeline(currentJobId);
      setRunning(false);
      addLog({
        timestamp: new Date(),
        level: 'warn',
        message: 'Pipeline cancelled by user',
      });
    } catch (err) {
      setError(String(err));
    }
  };

  useInput((input, key) => {
    if (isRunning) {
      if (key.escape) {
        cancelPipeline();
      }
      return;
    }

    if (key.upArrow) {
      setFocusedIndex(prev => Math.max(0, prev - 1));
    } else if (key.downArrow) {
      setFocusedIndex(prev => Math.min(STAGES.length - 1, prev + 1));
    } else if (input === ' ') {
      toggleStage(STAGES[focusedIndex]);
    } else if (key.return) {
      startPipeline();
    } else if (key.escape || input === 'q') {
      setScreen('welcome');
    }
  });

  const visibleLogs = logs.slice(-6);

  const runningShortcuts = [
    { key: 'Esc', label: 'Cancel' },
  ];

  const idleShortcuts = [
    { key: 'Space', label: 'Toggle' },
    { key: 'Enter', label: 'Start' },
    { key: 'Esc', label: 'Back' },
  ];

  return (
    <Box flexDirection="column" padding={1}>
      <Header title="Pipeline Runner" />

      {error && (
        <Box marginY={1}>
          <Text color="red">Error: {error}</Text>
        </Box>
      )}

      {!isBackendReady && (
        <Box marginY={1}>
          <Text color="yellow">Warning: Backend not connected</Text>
        </Box>
      )}

      {/* Input info */}
      {config && (
        <Box marginY={1}>
          <Text dimColor>Input: </Text>
          <Text color="cyan">{config.generated_mols_path}</Text>
        </Box>
      )}

      {/* Stage selection */}
      {!isRunning && (
        <Box flexDirection="column" marginY={1}>
          <Text color="cyan" bold>Select Stages:</Text>
          <Box flexDirection="column" marginTop={1}>
            {STAGES.map((stage, index) => {
              const isSelected = selectedStages.includes(stage);
              const isFocused = focusedIndex === index;
              return (
                <Box key={stage}>
                  <Text color={isFocused ? 'cyan' : 'white'}>
                    {isFocused ? '>' : ' '} [{isSelected ? 'x' : ' '}] {STAGE_NAMES[stage]}
                  </Text>
                </Box>
              );
            })}
          </Box>
        </Box>
      )}

      {/* Stage progress */}
      <Box flexDirection="column" marginY={1}>
        <Text color="cyan" bold>Progress:</Text>
        <Box flexDirection="column" marginTop={1}>
          {STAGES.filter(s => selectedStages.includes(s)).map((stage) => {
            const stageInfo = stages[stage];
            return (
              <ProgressBar
                key={stage}
                progress={stageInfo?.progress || 0}
                status={stageInfo?.status || 'pending'}
                label={STAGE_NAMES[stage]}
                width={30}
              />
            );
          })}
        </Box>
      </Box>

      {/* Logs */}
      <Box flexDirection="column" marginY={1}>
        <Text color="gray">{'─'.repeat(60)}</Text>
        <Box flexDirection="column" height={6}>
          {visibleLogs.map((log, index) => (
            <Box key={index}>
              <Text dimColor>[{formatTimestamp(log.timestamp)}]</Text>
              <Text color={log.level === 'error' ? 'red' : log.level === 'warn' ? 'yellow' : 'white'}>
                {' '}{log.message}
              </Text>
            </Box>
          ))}
          {visibleLogs.length === 0 && (
            <Text dimColor>Press Space to toggle stages, Enter to start.</Text>
          )}
        </Box>
      </Box>

      {/* Status */}
      {isRunning && (
        <Box marginY={1}>
          <Spinner label="Pipeline running..." />
        </Box>
      )}

      <Footer shortcuts={isRunning ? runningShortcuts : idleShortcuts} />
    </Box>
  );
}

export default PipelineRunner;
