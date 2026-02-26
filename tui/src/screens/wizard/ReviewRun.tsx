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

function truncateLine(value: string, maxLength: number): string {
  if (maxLength <= 0) return '';
  if (value.length <= maxLength) return value;
  if (maxLength <= 3) return '.'.repeat(maxLength);
  return `${value.slice(0, maxLength - 3)}...`;
}

export function ReviewRun(): React.ReactElement {
  const setScreen = useStore((state) => state.setScreen);
  const wizard = useStore((state) => state.wizard);
  const config = useStore((state) => state.configs.main);
  const getWizardSelectedStagesInOrder = useStore((state) => state.getWizardSelectedStagesInOrder);

  const { width: terminalWidth, height: terminalHeight } = useTerminalSize();

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

  const separatorWidth = Math.max(1, terminalWidth - 2);
  const maxPreflightItems = terminalHeight < 34 ? 3 : terminalHeight < 42 ? 5 : 8;
  const showTwoColumns = terminalWidth >= 92;
  const stageColumnWidth = showTwoColumns
    ? Math.max(26, Math.min(36, Math.floor(terminalWidth * 0.33)))
    : Math.max(20, terminalWidth - 2);
  const stageLabelWidth = Math.max(12, stageColumnWidth - 8);
  const focusedStage = stageViews[focusedIndex];

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
      <Header title="Pipeline Wizard" subtitle="Detailed Review" />

      {config && (
        <Box flexDirection="column" marginY={1}>
          <Box>
            <Text dimColor>Input:  </Text>
            <Text color="cyan" wrap="truncate-middle">{config.generated_mols_path || '(not set)'}</Text>
          </Box>
          <Box>
            <Text dimColor>Output: </Text>
            <Text color="cyan" wrap="truncate-middle">{config.folder_to_save || '(not set)'}</Text>
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
          {preflightItems.slice(0, maxPreflightItems).map((check, index) => {
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
          {preflightItems.length > maxPreflightItems && (
            <Text dimColor>... and {preflightItems.length - maxPreflightItems} more checks</Text>
          )}
        </Box>
      )}

      <Box marginBottom={1}>
        <Text color="cyan" bold>Pipeline Stages ({viewMode})</Text>
      </Box>

      <Box flexDirection={showTwoColumns ? 'row' : 'column'} marginY={1}>
        <Box
          flexDirection="column"
          width={stageColumnWidth}
          flexShrink={0}
          paddingRight={showTwoColumns ? 2 : 0}
          marginBottom={showTwoColumns ? 0 : 1}
        >
          <Text color="cyan" bold>Stages</Text>
          {stageViews.map((stage, index) => {
            const isFocused = focusedIndex === index;
            const label = truncateLine(stage.displayName, stageLabelWidth);

            return (
              <Box key={stage.name}>
                <Text color={isFocused ? 'cyan' : 'gray'}>
                  {isFocused ? '> ' : '  '}
                </Text>
                <Text color={isFocused ? 'white' : 'gray'}>{`${index + 1}. ${label}`}</Text>
              </Box>
            );
          })}
        </Box>

        <Box flexDirection="column" flexGrow={1} flexShrink={1}>
          {focusedStage ? (
            <>
              <Box>
                <Text color="white" bold>{focusedStage.displayName}</Text>
                <Text dimColor> ({formatHeavyLevel(focusedStage.heavyLevel)})</Text>
              </Box>
              <Text dimColor wrap="wrap">{focusedStage.summary}</Text>

              {viewMode === 'detailed' && (
                <Box flexDirection="column">
                  <Text dimColor wrap="wrap">Purpose: {focusedStage.whatItDoes}</Text>
                  <Text dimColor wrap="wrap">Reads: {focusedStage.reads.join(', ')}</Text>
                  <Text dimColor wrap="wrap">Writes: {focusedStage.writes.join(', ')}</Text>
                  {focusedStage.keyParams.length > 0 && (
                    <Text dimColor wrap="wrap">Key params: {focusedStage.keyParams.join(' | ')}</Text>
                  )}
                </Box>
              )}
            </>
          ) : (
            <Text dimColor>No stages selected</Text>
          )}
        </Box>
      </Box>

      <Box marginY={1}>
        <Text color="gray">{'─'.repeat(separatorWidth)}</Text>
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
