import React, { useEffect, useMemo, useState } from 'react';
import { Box, Text, useInput } from 'ink';
import { Header } from '../../components/Header.js';
import { Footer } from '../../components/Footer.js';
import { useStore } from '../../store/index.js';
import { useTerminalSize } from '../../hooks/useTerminalSize.js';
import { useTheme } from '../../theme/context.js';
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

function truncateMiddle(value: string, maxLength: number): string {
  if (maxLength <= 0) return '';
  if (value.length <= maxLength) return value;
  if (maxLength <= 3) return '.'.repeat(maxLength);
  const keep = maxLength - 3;
  const left = Math.ceil(keep / 2);
  const right = Math.floor(keep / 2);
  return `${value.slice(0, left)}...${value.slice(value.length - right)}`;
}

export function ReviewRun(): React.ReactElement {
  const { theme } = useTheme();
  const setScreen = useStore((state) => state.setScreen);
  const goBack = useStore((state) => state.goBack);
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
    ? { label: 'Not run', color: theme.palette.textMuted }
    : preflight.valid
      ? { label: 'Pass', color: theme.palette.success }
      : { label: 'Failed', color: theme.palette.error };

  const separatorWidth = Math.max(1, terminalWidth - 2);
  const maxPreflightItems = terminalHeight < 34 ? 3 : terminalHeight < 42 ? 5 : 8;
  const showTwoColumns = terminalWidth >= 92;
  const stageColumnWidth = showTwoColumns
    ? Math.max(26, Math.min(36, Math.floor(terminalWidth * 0.33)))
    : Math.max(20, terminalWidth - 2);
  const detailsColumnWidth = showTwoColumns
    ? Math.max(20, separatorWidth - stageColumnWidth - 2)
    : separatorWidth;
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
    } else if (key.escape || key.leftArrow) {
      goBack();
    }
  });

  return (
    <Box flexDirection="column" flexGrow={1} padding={1}>
      <Header title="Pipeline Wizard" subtitle="Detailed Review" />

      {config && (
        <Box flexDirection="column" marginY={1}>
          {(() => {
            const inputLabel = 'Input:  ';
            const outputLabel = 'Output: ';
            const inputValue = config.generated_mols_path || '(not set)';
            const outputValue = config.folder_to_save || '(not set)';

            const inputMax = Math.max(0, separatorWidth - inputLabel.length);
            const outputMax = Math.max(0, separatorWidth - outputLabel.length);
            const inputShown = truncateMiddle(inputValue, inputMax);
            const outputShown = truncateMiddle(outputValue, outputMax);

            const inputPad = Math.max(0, separatorWidth - (inputLabel.length + inputShown.length));
            const outputPad = Math.max(0, separatorWidth - (outputLabel.length + outputShown.length));

            return (
              <>
                <Box width={separatorWidth}>
                  <Text color={theme.palette.textMuted}>{inputLabel}</Text>
                  <Text color={theme.palette.primary}>{inputShown}</Text>
                  {inputPad > 0 && <Text>{' '.repeat(inputPad)}</Text>}
                </Box>
                <Box width={separatorWidth}>
                  <Text color={theme.palette.textMuted}>{outputLabel}</Text>
                  <Text color={theme.palette.primary}>{outputShown}</Text>
                  {outputPad > 0 && <Text>{' '.repeat(outputPad)}</Text>}
                </Box>
              </>
            );
          })()}
        </Box>
      )}

      <Box flexDirection="column" marginY={1}>
        <Text color={theme.palette.accent} bold>Preflight</Text>
        <Box width={separatorWidth}>
          <Text color={theme.palette.textMuted}>Status: </Text>
          <Text color={preflightStatus.color} bold>{preflightStatus.label}</Text>
          {refreshingPreflight && <Text color={theme.palette.warning}> (refreshing...)</Text>}
          {(() => {
            const suffix = refreshingPreflight ? ' (refreshing...)' : '';
            const used = 'Status: '.length + preflightStatus.label.length + suffix.length;
            const pad = Math.max(0, separatorWidth - used);
            return pad > 0 ? <Text>{' '.repeat(pad)}</Text> : null;
          })()}
        </Box>
        <Box>
          <Text color={theme.palette.textMuted}>Molecules: </Text>
          <Text color={theme.palette.text}>{preflight?.molecule_count?.toLocaleString() ?? 'Unknown'}</Text>
          <Text color={theme.palette.textMuted}> | Runtime: </Text>
          <Text color={theme.palette.text}>{runtimeLabel(preflight?.estimated_runtime || 'unknown')}</Text>
        </Box>
        <Box>
          <Text color={theme.palette.error}>Errors: {counters.errors}</Text>
          <Text color={theme.palette.textMuted}> | </Text>
          <Text color={theme.palette.warning}>Warnings: {counters.warnings}</Text>
          <Text color={theme.palette.textMuted}> | </Text>
          <Text color={theme.palette.info}>Info: {counters.infos}</Text>
        </Box>
      </Box>

      {preflightItems.length > 0 && (
        <Box flexDirection="column" marginBottom={1}>
          <Text color={theme.palette.accent} bold>Checks</Text>
          {preflightItems.slice(0, maxPreflightItems).map((check, index) => {
            const marker = check.level === 'error' ? 'E' : check.level === 'warning' ? 'W' : 'I';
            const color = check.level === 'error' ? theme.palette.error : check.level === 'warning' ? theme.palette.warning : theme.palette.info;
            const scope = check.scope === 'global' ? 'global' : check.scope;
            const prefix = `[${marker}] ${scope} · `;
            const messageMax = Math.max(0, separatorWidth - prefix.length);
            const messageShown = truncateLine(check.message, messageMax);
            const pad = Math.max(0, separatorWidth - (prefix.length + messageShown.length));

            return (
              <Box key={`${check.code}-${index}`} width={separatorWidth}>
                <Text color={color}>[{marker}]</Text>
                <Text color={theme.palette.textMuted}> {scope}</Text>
                <Text color={theme.palette.textMuted}> · </Text>
                <Text color={theme.palette.text} wrap="truncate-end">{messageShown}</Text>
                {pad > 0 && <Text>{' '.repeat(pad)}</Text>}
              </Box>
            );
          })}
          {preflightItems.length > maxPreflightItems && (
            <Box width={separatorWidth}>
              {(() => {
                const line = `... and ${preflightItems.length - maxPreflightItems} more checks`;
                const shown = truncateLine(line, separatorWidth);
                const pad = Math.max(0, separatorWidth - shown.length);
                return (
                  <>
                    <Text color={theme.palette.textMuted} wrap="truncate-end">{shown}</Text>
                    {pad > 0 && <Text>{' '.repeat(pad)}</Text>}
                  </>
                );
              })()}
            </Box>
          )}
        </Box>
      )}

      <Box marginBottom={1} width={separatorWidth}>
        {(() => {
          const line = `Pipeline Stages (${viewMode})`;
          const shown = truncateLine(line, separatorWidth);
          const pad = Math.max(0, separatorWidth - shown.length);
          return (
            <>
              <Text color={theme.palette.accent} bold wrap="truncate-end">{shown}</Text>
              {pad > 0 && <Text>{' '.repeat(pad)}</Text>}
            </>
          );
        })()}
      </Box>

      <Box flexDirection={showTwoColumns ? 'row' : 'column'} marginY={1} width={separatorWidth}>
        <Box
          flexDirection="column"
          width={stageColumnWidth}
          flexShrink={0}
          paddingRight={showTwoColumns ? 2 : 0}
          marginBottom={showTwoColumns ? 0 : 1}
        >
          <Text color={theme.palette.accent} bold>Stages</Text>
          {stageViews.map((stage, index) => {
            const isFocused = focusedIndex === index;
            const label = truncateLine(stage.displayName, stageLabelWidth);
            const prefix = `${index + 1}. `;
            const rowText = `${prefix}${label}`;
            const pad = Math.max(0, stageColumnWidth - (2 + rowText.length));

            return (
              <Box key={stage.name} width={stageColumnWidth}>
                <Text color={isFocused ? theme.palette.primary : theme.palette.textMuted}>
                  {isFocused ? '> ' : '  '}
                </Text>
                <Text color={isFocused ? theme.palette.text : theme.palette.textMuted} wrap="truncate-end">{rowText}</Text>
                {pad > 0 && <Text>{' '.repeat(pad)}</Text>}
              </Box>
            );
          })}
        </Box>

        <Box flexDirection="column" width={detailsColumnWidth} flexGrow={1} flexShrink={1}>
          {focusedStage ? (
            <>
              <Box>
                <Text color={theme.palette.text} bold>{focusedStage.displayName}</Text>
                <Text color={theme.palette.textMuted}> ({formatHeavyLevel(focusedStage.heavyLevel)})</Text>
              </Box>
              <Text color={theme.palette.textMuted} wrap="truncate-end">{focusedStage.summary}</Text>

              {viewMode === 'detailed' && (
                <Box flexDirection="column">
                  <Text color={theme.palette.textMuted} wrap="truncate-end">Purpose: {focusedStage.whatItDoes}</Text>
                  <Text color={theme.palette.textMuted} wrap="truncate-end">Reads: {focusedStage.reads.join(', ')}</Text>
                  <Text color={theme.palette.textMuted} wrap="truncate-end">Writes: {focusedStage.writes.join(', ')}</Text>
                  {focusedStage.keyParams.length > 0 && (
                    <Text color={theme.palette.textMuted} wrap="truncate-end">Key params: {focusedStage.keyParams.join(' | ')}</Text>
                  )}
                </Box>
              )}
            </>
          ) : (
            <Text color={theme.palette.textMuted}>No stages selected</Text>
          )}
        </Box>
      </Box>

      <Box marginY={1}>
        <Text color={theme.palette.border}>{'─'.repeat(separatorWidth)}</Text>
      </Box>

      {starting ? (
        <Box marginY={1} width={separatorWidth}>
          {(() => {
            const line = 'Running preflight and starting pipeline...';
            const shown = truncateLine(line, separatorWidth);
            const pad = Math.max(0, separatorWidth - shown.length);
            return (
              <>
                <Text color={theme.palette.warning} wrap="truncate-end">{shown}</Text>
                {pad > 0 && <Text>{' '.repeat(pad)}</Text>}
              </>
            );
          })()}
        </Box>
      ) : (
        <Box marginY={1} width={separatorWidth}>
          {(() => {
            const line = 'Press Enter to start the pipeline';
            const shown = truncateLine(line, separatorWidth);
            const pad = Math.max(0, separatorWidth - shown.length);
            return (
              <>
                <Text color={theme.palette.success} bold wrap="truncate-end">{shown}</Text>
                {pad > 0 && <Text>{' '.repeat(pad)}</Text>}
              </>
            );
          })()}
        </Box>
      )}

      <Footer />
    </Box>
  );
}

export default ReviewRun;
