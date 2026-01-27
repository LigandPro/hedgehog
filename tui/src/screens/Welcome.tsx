import React, { useState } from 'react';
import { Box, Text, useInput, useApp } from 'ink';
import { Header } from '../components/Header.js';
import { Footer } from '../components/Footer.js';
import { useStore } from '../store/index.js';
import { formatDate } from '../utils/format.js';
import { getStatusIcon, getStatusColor } from '../utils/job-status.js';
import type { Screen } from '../types/index.js';

interface QuickAction {
  key: string;
  label: string;
  action: Screen | 'exit';
}

const QUICK_ACTIONS: QuickAction[] = [
  { key: 'n', label: 'New Pipeline Run', action: 'wizardInputSelection' },
  { key: 'h', label: 'View Job History', action: 'history' },
  { key: 'q', label: 'Quit', action: 'exit' },
];

export function Welcome(): React.ReactElement {
  const { exit } = useApp();
  const setScreen = useStore((state) => state.setScreen);
  const setSelectedJob = useStore((state) => state.setSelectedJob);
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

  useInput((input, key) => {
    // Quick action shortcuts
    for (const item of QUICK_ACTIONS) {
      if (input === item.key) {
        if (item.action === 'exit') {
          exit();
        } else {
          setScreen(item.action);
        }
        return;
      }
    }

    // [p] Go to running pipeline
    if (input === 'p' && isRunning) {
      setScreen('pipelineRunner');
      return;
    }

    // Delete job with 'd' when in recent section
    if (input === 'd' && section === 'recent') {
      deleteSelectedJob();
      return;
    }

    // Navigation
    if (key.upArrow) {
      if (section === 'actions') {
        setFocusedIndex(Math.max(0, focusedIndex - 1));
      } else if (section === 'recent') {
        if (focusedIndex === 0) {
          setSection('actions');
          setFocusedIndex(QUICK_ACTIONS.length - 1);
        } else {
          setFocusedIndex(focusedIndex - 1);
        }
      }
    } else if (key.downArrow) {
      if (section === 'actions') {
        if (focusedIndex === QUICK_ACTIONS.length - 1 && recentJobs.length > 0) {
          setSection('recent');
          setFocusedIndex(0);
        } else {
          setFocusedIndex(Math.min(QUICK_ACTIONS.length - 1, focusedIndex + 1));
        }
      } else if (section === 'recent') {
        setFocusedIndex(Math.min(recentJobs.length - 1, focusedIndex + 1));
      }
    } else if (key.return || key.rightArrow) {
      if (section === 'actions') {
        const item = QUICK_ACTIONS[focusedIndex];
        if (item.action === 'exit') {
          exit();
        } else {
          setScreen(item.action);
        }
      } else if (section === 'recent' && recentJobs[focusedIndex]) {
        setSelectedJob(recentJobs[focusedIndex].id);
        setScreen('results');
      }
    } else if (key.escape) {
      exit();
    }
  });

  return (
    <Box flexDirection="column" padding={1}>
      <Header showLogo={true} />

      {!isBackendReady && (
        <Box marginY={1}>
          <Text color="yellow">⚠ Backend not connected. Some features may not work.</Text>
        </Box>
      )}

      {/* Running Pipeline indicator */}
      {isRunning && (
        <Box flexDirection="column" marginY={1}>
          <Box>
            <Text color="yellow" bold>● Pipeline Running</Text>
            <Text dimColor>  Press </Text>
            <Text color="cyan" bold>[p]</Text>
            <Text dimColor> to view progress</Text>
          </Box>
        </Box>
      )}

      {/* Quick Actions Section (Matcha pattern) */}
      <Box flexDirection="column" marginY={1}>
        <Text color="cyan" bold>Quick Actions</Text>
        <Box flexDirection="column" marginTop={1}>
          {QUICK_ACTIONS.map((item, index) => {
            const isSelected = section === 'actions' && focusedIndex === index;
            return (
              <Box key={item.key}>
                <Text color={isSelected ? 'cyan' : 'gray'}>{isSelected ? '▶ ' : '  '}</Text>
                <Text color="cyan" bold>[{item.key}]</Text>
                <Text color={isSelected ? 'white' : 'gray'}> {item.label}</Text>
              </Box>
            );
          })}
        </Box>
      </Box>

      {/* Recent Jobs Section (Matcha pattern) */}
      {recentJobs.length > 0 && (
        <Box flexDirection="column" marginY={1}>
          <Box>
            <Text color="cyan" bold>Recent Jobs</Text>
            {section === 'recent' && (
              <Text dimColor>  [Enter] View  [d] Delete</Text>
            )}
          </Box>
          <Box flexDirection="column" marginTop={1}>
            {recentJobs.map((job, index) => {
              const isSelected = section === 'recent' && focusedIndex === index;
              const statusIcon = getStatusIcon(job.status);
              const statusColor = getStatusColor(job.status);

              return (
                <Box key={job.id}>
                  <Text color={isSelected ? 'cyan' : 'gray'}>{isSelected ? '> ' : '  '}</Text>
                  <Text color={statusColor}>{statusIcon}</Text>
                  <Text color={isSelected ? 'white' : 'gray'}>
                    {' '}{job.name || job.id}
                  </Text>
                  <Text dimColor>
                    {' - '}{formatDate(job.startTime)}
                  </Text>
                </Box>
              );
            })}
          </Box>
        </Box>
      )}

      {/* Footer description (Matcha pattern) */}
      <Box marginTop={1}>
        <Text dimColor>Molecular design pipeline with filtering, synthesis scoring, and docking</Text>
      </Box>

      <Footer />
    </Box>
  );
}

export default Welcome;
