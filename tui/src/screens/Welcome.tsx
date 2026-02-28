import React, { useState } from 'react';
import { Box, Text, useInput, useApp } from 'ink';
import { Footer } from '../components/Footer.js';
import { useStore } from '../store/index.js';
import { useTerminalSize } from '../hooks/useTerminalSize.js';
import { useTheme } from '../theme/context.js';
import { formatDate } from '../utils/format.js';
import { getStatusIcon } from '../utils/job-status.js';
import { buildExitSummary } from '../utils/exit-summary.js';
import type { Screen } from '../types/index.js';

interface QuickAction {
  key: string;
  label: string;
  action: Screen | 'exit';
}

const QUICK_ACTIONS: QuickAction[] = [
  { key: 'n', label: 'New Pipeline Run', action: 'wizardInputSelection' },
  { key: 'h', label: 'View Job History', action: 'history' },
  { key: 'c', label: 'Open Configuration', action: 'configMain' },
  { key: 'q', label: 'Quit', action: 'exit' },
];

const WELCOME_LOGO = [
  '█  █ █▀▀ █▀▀▄ █▀▀▀ █▀▀ █  █ █▀▀█ █▀▀▀',
  '█▀▀█ █▀▀ █  █ █ ▀█ █▀▀ █▀▀█ █  █ █ ▀█',
  '▀  ▀ ▀▀▀ ▀▀▀  ▀▀▀▀ ▀▀▀ ▀  ▀ ▀▀▀▀ ▀▀▀▀',
];

export function Welcome(): React.ReactElement {
  const { exit } = useApp();
  const { width } = useTerminalSize();
  const { theme } = useTheme();
  const setScreen = useStore((state) => state.setScreen);
  const screen = useStore((state) => state.screen);
  const setSelectedJob = useStore((state) => state.setSelectedJob);
  const selectedJobId = useStore((state) => state.selectedJobId);
  const currentJobId = useStore((state) => state.currentJobId);
  const setExitSummary = useStore((state) => state.setExitSummary);
  const isBackendReady = useStore((state) => state.isBackendReady);
  const isRunning = useStore((state) => state.isRunning);
  const jobHistory = useStore((state) => state.jobHistory);
  const removeJobFromHistory = useStore((state) => state.removeJobFromHistory);
  const showConfirm = useStore((state) => state.showConfirm);
  const showToast = useStore((state) => state.showToast);
  const [focusedIndex, setFocusedIndex] = useState(0);
  const [section, setSection] = useState<'actions' | 'recent'>('actions');

  const recentJobs = jobHistory.slice(0, 5);

  const deleteSelectedJob = () => {
    if (section !== 'recent' || !recentJobs[focusedIndex]) return;

    const job = recentJobs[focusedIndex];
    showConfirm({
      title: 'Delete Job?',
      message: `Remove "${job.name || job.id}" from history?`,
      confirmLabel: 'Delete',
      cancelLabel: 'Cancel',
      onConfirm: () => {
        removeJobFromHistory(job.id);
        showToast('info', 'Job removed');
        // Adjust focus if needed
        if (focusedIndex >= recentJobs.length - 1 && focusedIndex > 0) {
          setFocusedIndex(focusedIndex - 1);
        }
        if (recentJobs.length <= 1) {
          setSection('actions');
          setFocusedIndex(0);
        }
      },
    });
  };

  const executeAction = (action: Screen | 'exit') => {
    if (action === 'exit') {
      setExitSummary(buildExitSummary({
        currentJobId,
        selectedJobId,
        jobHistory,
        screen,
      }));
      exit();
    } else {
      setScreen(action);
    }
  };

  const handleUpNavigation = () => {
    if (section === 'actions') {
      setFocusedIndex(Math.max(0, focusedIndex - 1));
    } else if (focusedIndex === 0) {
      setSection('actions');
      setFocusedIndex(QUICK_ACTIONS.length - 1);
    } else {
      setFocusedIndex(focusedIndex - 1);
    }
  };

  const handleDownNavigation = () => {
    if (section === 'actions') {
      if (focusedIndex === QUICK_ACTIONS.length - 1 && recentJobs.length > 0) {
        setSection('recent');
        setFocusedIndex(0);
      } else {
        setFocusedIndex(Math.min(QUICK_ACTIONS.length - 1, focusedIndex + 1));
      }
    } else {
      setFocusedIndex(Math.min(recentJobs.length - 1, focusedIndex + 1));
    }
  };

  const handleSelectAction = () => {
    if (section === 'actions') {
      executeAction(QUICK_ACTIONS[focusedIndex].action);
    } else if (section === 'recent' && recentJobs[focusedIndex]) {
      setSelectedJob(recentJobs[focusedIndex].id);
      setScreen('results');
    }
  };

  useInput((input, key) => {
    // Quick action shortcuts
    for (const item of QUICK_ACTIONS) {
      if (input === item.key) {
        executeAction(item.action);
        return;
      }
    }

    if (input === 'p' && isRunning) {
      setScreen('pipelineRunner');
      return;
    }

    if (input === 'd' && section === 'recent') {
      deleteSelectedJob();
      return;
    }

    if (key.upArrow) {
      handleUpNavigation();
    } else if (key.downArrow) {
      handleDownNavigation();
    } else if (key.return || key.rightArrow) {
      handleSelectAction();
    } else if (key.escape) {
      exit();
    }
  });

  return (
    <Box flexDirection="column" flexGrow={1} paddingY={1} minHeight={12}>
      <Box flexGrow={1} />

      <Box justifyContent="center">
        <Box flexDirection="column" width={Math.max(1, Math.min(86, width - 4))}>
          <Box marginBottom={1} justifyContent="center">
            <Box flexDirection="column">
              {WELCOME_LOGO.map((line, index) => (
                <Text key={index} color="white" bold>{line}</Text>
              ))}
            </Box>
          </Box>

          <Box flexDirection="column" paddingX={2} paddingY={1}>
            {!isBackendReady && (
              <Box marginBottom={1}>
                <Text color={theme.palette.warning}>⚠ Backend not connected. Some features may not work.</Text>
              </Box>
            )}

            {isRunning && (
              <Box marginBottom={1}>
                <Text color={theme.palette.warning} bold>● Pipeline Running</Text>
                <Text color={theme.palette.textMuted}>  Press </Text>
                <Text color={theme.palette.primary} bold>[p]</Text>
                <Text color={theme.palette.textMuted}> to view progress</Text>
              </Box>
            )}

            <Text color={theme.palette.accent} bold>Quick Actions</Text>
            <Box flexDirection="column" marginTop={1} marginBottom={1}>
              {QUICK_ACTIONS.map((item, index) => {
                const isSelected = section === 'actions' && focusedIndex === index;
                return (
                  <Box key={item.key} paddingX={1}>
                    <Text
                      color={theme.palette.textMuted}
                      backgroundColor={isSelected ? theme.palette.panelStrong : undefined}
                    >
                      {isSelected ? '› ' : '  '}
                    </Text>
                    <Text color={theme.palette.accent} bold backgroundColor={isSelected ? theme.palette.panelStrong : undefined}>
                      [{item.key}]
                    </Text>
                    <Text
                      color={isSelected ? theme.palette.text : theme.palette.textMuted}
                      backgroundColor={isSelected ? theme.palette.panelStrong : undefined}
                    >
                      {' '}{item.label}
                    </Text>
                  </Box>
                );
              })}
            </Box>

            {recentJobs.length > 0 && (
              <>
                <Text color={theme.palette.accent} bold>Recent Jobs</Text>
                <Box flexDirection="column" marginTop={1}>
                  {recentJobs.map((job, index) => {
                    const isSelected = section === 'recent' && focusedIndex === index;
                    const statusIcon = getStatusIcon(job.status);
                    const statusColor = job.status === 'completed'
                      ? theme.palette.success
                      : job.status === 'running'
                        ? theme.palette.warning
                        : job.status === 'error'
                          ? theme.palette.error
                          : theme.palette.textMuted;

                    return (
                      <Box key={job.id} paddingX={1}>
                        <Text
                          color={theme.palette.textMuted}
                          backgroundColor={isSelected ? theme.palette.panelStrong : undefined}
                        >
                          {isSelected ? '› ' : '  '}
                        </Text>
                        <Text color={statusColor} backgroundColor={isSelected ? theme.palette.panelStrong : undefined}>
                          {statusIcon}
                        </Text>
                        <Text
                          color={isSelected ? theme.palette.text : theme.palette.textMuted}
                          backgroundColor={isSelected ? theme.palette.panelStrong : undefined}
                        >
                          {' '}{job.name || job.id}
                        </Text>
                        <Text color={theme.palette.textMuted} backgroundColor={isSelected ? theme.palette.panelStrong : undefined}>
                          {' - '}{formatDate(job.startTime)}
                        </Text>
                      </Box>
                    );
                  })}
                </Box>
              </>
            )}
          </Box>

          <Box justifyContent="center" marginTop={1}>
            <Text color={theme.palette.text} bold>HEDGEHOG PIPELINE ENGINE</Text>
          </Box>
          <Box justifyContent="center">
            <Text color={theme.palette.textMuted}>Type </Text>
            <Text color="white">/</Text>
            <Text color={theme.palette.textMuted}> to open command palette</Text>
          </Box>
        </Box>
      </Box>

      <Box flexGrow={1} />
      <Footer showBreadcrumbs={false} />
    </Box>
  );
}

export default Welcome;
