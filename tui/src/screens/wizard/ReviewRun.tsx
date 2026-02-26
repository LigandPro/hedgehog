import React, { useEffect, useMemo, useState } from 'react';
import { Box, Text, useInput } from 'ink';
import { Header } from '../../components/Header.js';
import { Footer } from '../../components/Footer.js';
import { useStore } from '../../store/index.js';
import { useTerminalSize } from '../../hooks/useTerminalSize.js';
import { countPreflightChecks, refreshWizardPreflight, runWizardPipeline } from './runPipeline.js';
import { getStageKeyParams, getStageMetadata, getStageScreen, getStageSummary } from './stageMetadata.js';

interface StageView {
  name: string;
  displayName: string;
  summary: string;
  whatItDoes: string;
  reads: string[];
  writes: string[];
  heavyLevel: 'low' | 'medium' | 'high';
  keyParams: string[];
}

function formatHeavyLevel(level: 'low' | 'medium' | 'high'): string {
  if (level === 'high') return 'high load';
  if (level === 'medium') return 'medium load';
  return 'low load';
}

function runtimeLabel(runtime: 'short' | 'medium' | 'long' | 'unknown'): string {
  if (runtime === 'short') return 'Short';
  if (runtime === 'medium') return 'Medium';
  if (runtime === 'long') return 'Long';
  return 'Unknown';
}

export function ReviewRun(): React.ReactElement {
  const setScreen = useStore((state) => state.setScreen);
  const wizard = useStore((state) => state.wizard);
  const config = useStore((state) => state.configs.main);
  const getWizardSelectedStagesInOrder = useStore((state) => state.getWizardSelectedStagesInOrder);

  const { width: terminalWidth } = useTerminalSize();

  const [focusedIndex, setFocusedIndex] = useState(0);
  const [starting, setStarting] = useState(false);
  const [refreshingPreflight, setRefreshingPreflight] = useState(false);
  const [viewMode, setViewMode] = useState<'summary' | 'detailed'>('detailed');

  const selectedStagesInOrder = useMemo(
    () => getWizardSelectedStagesInOrder(),
    [wizard.stageOrder, wizard.selectedStages, getWizardSelectedStagesInOrder]
  );

  const stageViews = useMemo((): StageView[] => {
    return selectedStagesInOrder.map((stageName) => {
      const params = wizard.stageConfigs[stageName]?.quickParams || {};
      const preset = wizard.stageConfigs[stageName]?.preset;
      const metadata = getStageMetadata(stageName);

      return {
        name: stageName,
        displayName: metadata?.title || stageName,
        summary: getStageSummary(stageName, params, preset),
        whatItDoes: metadata?.whatItDoes || 'Configured stage.',
        reads: metadata?.reads || ['-'],
        writes: metadata?.writes || ['-'],
        heavyLevel: metadata?.heavyLevel || 'medium',
        keyParams: getStageKeyParams(stageName, params),
      };
    });
  }, [selectedStagesInOrder, wizard.stageConfigs]);

  const preflight = wizard.preflight;
  const counters = useMemo(() => countPreflightChecks(preflight), [preflight]);

  const preflightItems = useMemo(() => {
    if (!preflight) return [];

    const globalChecks = preflight.checks.map((check) => ({
      scope: 'global',
      ...check,
    }));

    const stageChecks = preflight.stage_reports.flatMap((report) =>
      report.checks.map((check) => ({
        scope: report.stage,
        ...check,
      }))
    );

    return [...globalChecks, ...stageChecks];
  }, [preflight]);

  useEffect(() => {
    setFocusedIndex((value) => {
      if (stageViews.length === 0) return 0;
      return Math.min(value, stageViews.length - 1);
    });
  }, [stageViews.length]);

  useEffect(() => {
    let isActive = true;

    const run = async () => {
      setRefreshingPreflight(true);
      try {
        await refreshWizardPreflight(selectedStagesInOrder);
      } finally {
        if (isActive) {
          setRefreshingPreflight(false);
        }
      }
    };

    void run();

    return () => {
      isActive = false;
    };
  }, [selectedStagesInOrder.join(',')]);

  const refreshPreflight = async () => {
    if (refreshingPreflight) return;
    setRefreshingPreflight(true);
    try {
      await refreshWizardPreflight(selectedStagesInOrder);
    } finally {
      setRefreshingPreflight(false);
    }
  };

  const startPipeline = async () => {
    if (starting) return;

    setStarting(true);
    try {
      await runWizardPipeline({ onPreflightErrorScreen: 'wizardReview' });
    } finally {
      setStarting(false);
    }
  };

  const preflightStatus = !preflight
    ? { label: 'Not run', color: 'gray' as const }
    : preflight.valid
      ? { label: 'Pass', color: 'green' as const }
      : { label: 'Failed', color: 'red' as const };

  useInput((input, key) => {
    if (starting) return;

    if (key.upArrow) {
      if (stageViews.length === 0) return;
      setFocusedIndex(Math.max(0, focusedIndex - 1));
    } else if (key.downArrow) {
      if (stageViews.length === 0) return;
      setFocusedIndex(Math.min(stageViews.length - 1, focusedIndex + 1));
    } else if (key.tab) {
      setViewMode(viewMode === 'detailed' ? 'summary' : 'detailed');
    } else if (input === 'e' && stageViews[focusedIndex]) {
      setScreen(getStageScreen(stageViews[focusedIndex].name));
    } else if (input === 'r') {
      void refreshPreflight();
    } else if (key.return) {
      void startPipeline();
    } else if (key.escape || key.leftArrow || input === 'q') {
      setScreen('wizardStageSelection');
    }
  });

  const shortcuts = [
    { key: '↑↓', label: 'Navigate stages' },
    { key: 'e', label: 'Edit stage' },
    { key: 'r', label: 'Refresh preflight' },
    { key: 'Tab', label: 'Summary/Detailed' },
    { key: 'Enter', label: 'Start' },
    { key: '←/Esc', label: 'Back' },
  ];

  return (
    <Box flexDirection="column" padding={1}>
      <Header title="Pipeline Wizard" subtitle="Detailed Review (Optional)" />

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

      <Box flexDirection="column" marginY={1}>
        <Text color="cyan" bold>Preflight</Text>
        <Box>
          <Text dimColor>Status: </Text>
          <Text color={preflightStatus.color} bold>{preflightStatus.label}</Text>
          {refreshingPreflight && <Text color="yellow"> (refreshing...)</Text>}
        </Box>
        <Box>
          <Text dimColor>Molecules: </Text>
          <Text color="white">{preflight?.molecule_count?.toLocaleString() ?? 'Unknown'}</Text>
          <Text dimColor> | Runtime: </Text>
          <Text color="white">{runtimeLabel(preflight?.estimated_runtime || 'unknown')}</Text>
        </Box>
        <Box>
          <Text color="red">Errors: {counters.errors}</Text>
          <Text dimColor> | </Text>
          <Text color="yellow">Warnings: {counters.warnings}</Text>
          <Text dimColor> | </Text>
          <Text color="cyan">Info: {counters.infos}</Text>
        </Box>
      </Box>

      {preflightItems.length > 0 && (
        <Box flexDirection="column" marginBottom={1}>
          <Text color="cyan" bold>Checks</Text>
          {preflightItems.slice(0, 8).map((check, index) => {
            const marker = check.level === 'error' ? 'E' : check.level === 'warning' ? 'W' : 'I';
            const color = check.level === 'error' ? 'red' : check.level === 'warning' ? 'yellow' : 'cyan';
            const scope = check.scope === 'global' ? 'global' : check.scope;

            return (
              <Box key={`${check.code}-${index}`}>
                <Text color={color}>[{marker}]</Text>
                <Text dimColor> {scope}</Text>
                <Text dimColor> · </Text>
                <Text>{check.message}</Text>
              </Box>
            );
          })}
          {preflightItems.length > 8 && (
            <Text dimColor>... and {preflightItems.length - 8} more checks</Text>
          )}
        </Box>
      )}

      <Box marginBottom={1}>
        <Text color="cyan" bold>Pipeline Stages ({viewMode})</Text>
      </Box>

      <Box flexDirection="column" marginY={1}>
        {stageViews.map((stage, index) => {
          const isFocused = focusedIndex === index;

          return (
            <Box key={stage.name} flexDirection="column" marginBottom={1}>
              <Box>
                <Text dimColor>{`${index + 1}. `}</Text>
                <Text color={isFocused ? 'cyan' : 'white'}>{isFocused ? '> ' : '  '}</Text>
                <Text color={isFocused ? 'white' : 'gray'} bold>{stage.displayName}</Text>
                <Text dimColor> ({formatHeavyLevel(stage.heavyLevel)})</Text>
              </Box>

              <Box paddingLeft={5}>
                <Text dimColor>{stage.summary}</Text>
              </Box>

              {viewMode === 'detailed' && (
                <Box flexDirection="column" paddingLeft={5}>
                  <Text dimColor>Purpose: {stage.whatItDoes}</Text>
                  <Text dimColor>Reads: {stage.reads.join(', ')}</Text>
                  <Text dimColor>Writes: {stage.writes.join(', ')}</Text>
                  {stage.keyParams.length > 0 && (
                    <Text dimColor>Key params: {stage.keyParams.join(' | ')}</Text>
                  )}
                </Box>
              )}
            </Box>
          );
        })}
      </Box>

      <Box marginY={1}>
        <Text color="gray">{'─'.repeat(terminalWidth - 2)}</Text>
      </Box>

      {starting ? (
        <Box marginY={1}>
          <Text color="yellow">Running preflight and starting pipeline...</Text>
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
