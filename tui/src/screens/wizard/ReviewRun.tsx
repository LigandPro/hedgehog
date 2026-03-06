import React, { useEffect, useMemo, useState } from 'react';
import { Box, Text, useInput } from 'ink';
import { Header } from '../../components/Header.js';
import { Footer } from '../../components/Footer.js';
import { AppShell } from '../../components/AppShell.js';
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
    if (preflight && !preflight.valid) return;

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
  const isPreflightBlocked = Boolean(preflight && !preflight.valid);

  const separatorWidth = Math.max(1, terminalWidth - 2);
  const showTwoColumns = terminalWidth >= 80;
  const columnGap = showTwoColumns ? 2 : 0;
  const stageColumnWidth = showTwoColumns
    ? Math.max(26, Math.min(36, Math.floor(separatorWidth * 0.34)))
    : separatorWidth;
  const detailsColumnWidth = showTwoColumns
    ? Math.max(20, separatorWidth - stageColumnWidth - columnGap)
    : separatorWidth;
  const stageLabelWidth = Math.max(12, stageColumnWidth - 8);
  const focusedStage = stageViews[focusedIndex];

  const hasChecks = preflightItems.length > 0;
  const shellReservedRows = 12;
  const bodyRowsBudget = Math.max(8, terminalHeight - shellReservedRows);
  const fixedRows = (config ? 2 : 0) + 2 + (hasChecks ? 1 : 0) + 1 + 2 + 1;
  const variableRows = Math.max(2, bodyRowsBudget - fixedRows);

  const preferredStageRows = showTwoColumns ? 7 : 8;
  let stageRowsBudget = Math.max(2, Math.min(variableRows, preferredStageRows));
  let checksRowsBudget = hasChecks ? Math.max(0, variableRows - stageRowsBudget) : 0;
  if (hasChecks && checksRowsBudget === 0) {
    stageRowsBudget = Math.max(1, stageRowsBudget - 1);
    checksRowsBudget = 1;
  }

  const canRenderChecksOverflow = checksRowsBudget > 1;
  const preflightNeedsOverflow = hasChecks && canRenderChecksOverflow && preflightItems.length > checksRowsBudget;
  const visiblePreflightItems = !hasChecks || checksRowsBudget <= 0
    ? []
    : preflightNeedsOverflow
      ? preflightItems.slice(0, Math.max(1, checksRowsBudget - 1))
      : preflightItems.slice(0, checksRowsBudget);
  const hiddenPreflightItems = preflightItems.length - visiblePreflightItems.length;

  const stageListRows = showTwoColumns
    ? Math.max(1, stageRowsBudget - 1)
    : Math.max(1, Math.min(stageViews.length, Math.max(1, Math.floor((stageRowsBudget - 2) / 2))));
  const stageWindowStart = stageViews.length <= stageListRows
    ? 0
    : Math.min(
      Math.max(0, focusedIndex - Math.floor(stageListRows / 2)),
      stageViews.length - stageListRows
    );
  const visibleStages = stageViews.slice(stageWindowStart, stageWindowStart + stageListRows);

  const detailLineBudget = showTwoColumns
    ? Math.max(1, stageRowsBudget - 1)
    : Math.max(1, stageRowsBudget - stageListRows - 2);
  const detailLines = focusedStage
    ? [
      focusedStage.summary,
      ...(viewMode === 'detailed'
        ? [
          `Purpose: ${focusedStage.whatItDoes}`,
          `Reads: ${focusedStage.reads.join(', ')}`,
          `Writes: ${focusedStage.writes.join(', ')}`,
          ...(focusedStage.keyParams.length > 0
            ? [`Key params: ${focusedStage.keyParams.join(' | ')}`]
            : []),
        ]
        : []),
    ]
    : ['No stages selected'];

  let visibleDetailLines = detailLines.slice(0, detailLineBudget);
  const hiddenDetailLines = detailLines.length - visibleDetailLines.length;
  if (hiddenDetailLines > 0 && visibleDetailLines.length > 0) {
    const suffix = `... and ${hiddenDetailLines} more line${hiddenDetailLines === 1 ? '' : 's'}`;
    visibleDetailLines = [...visibleDetailLines.slice(0, -1), suffix];
  }

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
      if (!isPreflightBlocked) {
        void startPipeline();
      }
    } else if (key.escape || key.leftArrow) {
      goBack();
    }
  });

  return (
    <AppShell
      padding={1}
      header={<Header title="Pipeline Wizard" subtitle="Detailed Review" />}
      footer={<Footer />}
    >

      {config && (
        <Box flexDirection="column">
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

      <Box flexDirection="column">
        <Box width={separatorWidth}>
          <Text color={theme.palette.accent} bold>Preflight: </Text>
          <Text color={preflightStatus.color} bold>{preflightStatus.label}</Text>
          {refreshingPreflight && <Text color={theme.palette.warning}> (refreshing...)</Text>}
          {(() => {
            const suffix = refreshingPreflight ? ' (refreshing...)' : '';
            const used = 'Preflight: '.length + preflightStatus.label.length + suffix.length;
            const pad = Math.max(0, separatorWidth - used);
            return pad > 0 ? <Text>{' '.repeat(pad)}</Text> : null;
          })()}
        </Box>
        <Box width={separatorWidth}>
          {(() => {
            const line = `Molecules: ${preflight?.molecule_count?.toLocaleString() ?? 'Unknown'} | Runtime: ${runtimeLabel(preflight?.estimated_runtime || 'unknown')} | E:${counters.errors} W:${counters.warnings} I:${counters.infos}`;
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
      </Box>

      {hasChecks && (
        <Box flexDirection="column">
          <Box width={separatorWidth}>
            {(() => {
              const line = `Checks (${preflightItems.length})`;
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
          {visiblePreflightItems.map((check, index) => {
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
          {canRenderChecksOverflow && hiddenPreflightItems > 0 && (
            <Box width={separatorWidth}>
              {(() => {
                const line = `... and ${hiddenPreflightItems} more checks`;
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

      <Box width={separatorWidth}>
        <Text color={theme.palette.border}>{'─'.repeat(separatorWidth)}</Text>
      </Box>

      <Box width={separatorWidth}>
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

      <Box flexDirection={showTwoColumns ? 'row' : 'column'} width={separatorWidth}>
        <Box
          flexDirection="column"
          width={stageColumnWidth}
          flexShrink={0}
          paddingRight={showTwoColumns ? columnGap : 0}
        >
          <Text color={theme.palette.accent} bold>Stages</Text>
          {visibleStages.map((stage, localIndex) => {
            const index = stageWindowStart + localIndex;
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
            <Box width={detailsColumnWidth}>
              <Text color={theme.palette.text} bold wrap="truncate-end">{truncateLine(focusedStage.displayName, detailsColumnWidth)}</Text>
              <Text color={theme.palette.textMuted}> ({formatHeavyLevel(focusedStage.heavyLevel)})</Text>
            </Box>
          ) : (
            <Text color={theme.palette.textMuted}>No stages selected</Text>
          )}
          {visibleDetailLines.map((line, index) => {
            const shown = truncateLine(line, detailsColumnWidth);
            const pad = Math.max(0, detailsColumnWidth - shown.length);
            return (
              <Box key={`detail-${index}`} width={detailsColumnWidth}>
                <Text color={theme.palette.textMuted} wrap="truncate-end">{shown}</Text>
                {pad > 0 && <Text>{' '.repeat(pad)}</Text>}
              </Box>
            );
          })}
        </Box>
      </Box>

      <Box width={separatorWidth}>
        <Text color={theme.palette.border}>{'─'.repeat(separatorWidth)}</Text>
      </Box>

      {starting ? (
        <Box width={separatorWidth}>
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
        <Box width={separatorWidth}>
          {(() => {
            const line = isPreflightBlocked
              ? 'Preflight failed. Fix errors and press r to refresh.'
              : 'Press Enter to start the pipeline';
            const shown = truncateLine(line, separatorWidth);
            const pad = Math.max(0, separatorWidth - shown.length);
            return (
              <>
                <Text
                  color={isPreflightBlocked ? theme.palette.error : theme.palette.success}
                  bold
                  wrap="truncate-end"
                >
                  {shown}
                </Text>
                {pad > 0 && <Text>{' '.repeat(pad)}</Text>}
              </>
            );
          })()}
        </Box>
      )}
    </AppShell>
  );
}

export default ReviewRun;
