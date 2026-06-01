import React, { useEffect, useState } from 'react';
import { Box, Text, useInput } from 'ink';
import { Header } from '../components/Header.js';
import { Footer } from '../components/Footer.js';
import { AppShell } from '../components/AppShell.js';
import { useStore } from '../store/index.js';
import { useTheme } from '../theme/context.js';
import { formatTimestamp, formatDuration } from '../utils/format.js';
import { createReadStream, existsSync } from 'fs';
import { access } from 'fs/promises';
import { createInterface } from 'readline';
import { dirname, isAbsolute, join, resolve } from 'path';

export function Results(): React.ReactElement {
  const { theme } = useTheme();
  const setScreen = useStore((state) => state.setScreen);
  const goBack = useStore((state) => state.goBack);
  const selectedJobId = useStore((state) => state.selectedJobId);
  const jobHistory = useStore((state) => state.jobHistory);

  const job = jobHistory.find((j) => j.id === selectedJobId);
  const [stageStats, setStageStats] = useState<Array<{
    label: string;
    count: number | null;
    percent: number | null;
  }> | null>(null);
  const [stageStatsError, setStageStatsError] = useState<string | null>(null);

  const findProjectRoot = (): string => {
    let current = process.cwd();
    while (true) {
      if (existsSync(join(current, 'pyproject.toml'))) {
        return current;
      }
      const parent = dirname(current);
      if (parent === current) {
        return process.cwd();
      }
      current = parent;
    }
  };

  useInput((input, key) => {
    if (key.escape || key.leftArrow) {
      goBack();
    } else if (input === 'r' && job) {
      // Re-run with same config - go to pipeline runner
      setScreen('pipelineRunner');
    }
  });

  useEffect(() => {
    if (!job) {
      setStageStats(null);
      setStageStatsError(null);
      return;
    }

    let cancelled = false;

    const countCsvRows = async (filePath: string): Promise<number> => {
      await access(filePath);
      return new Promise((resolveCount, reject) => {
        let lineCount = 0;
        const stream = createReadStream(filePath);
        const rl = createInterface({ input: stream, crlfDelay: Infinity });

        rl.on('line', (line) => {
          if (line.trim().length > 0) {
            lineCount += 1;
          }
        });
        rl.on('close', () => {
          resolveCount(Math.max(0, lineCount - 1));
        });
        rl.on('error', reject);
        stream.on('error', reject);
      });
    };

    const findExistingPath = async (basePath: string, candidates: string[]): Promise<string | null> => {
      for (const relPath of candidates) {
        const fullPath = resolve(basePath, relPath);
        try {
          await access(fullPath);
          return fullPath;
        } catch {
          // Try next candidate
        }
      }
      return null;
    };

    const loadStageStats = async () => {
      try {
        const outputPath = job.config.outputPath;
        if (!outputPath) {
          setStageStats(null);
          return;
        }

        const projectRoot = findProjectRoot();
        const basePath = isAbsolute(outputPath) ? outputPath : resolve(projectRoot, outputPath);
        let initialCount = job.results?.moleculesProcessed ?? 0;

        if (!initialCount && job.config.inputPath) {
          const inputPath = isAbsolute(job.config.inputPath)
            ? job.config.inputPath
            : resolve(projectRoot, job.config.inputPath);
          try {
            initialCount = await countCsvRows(inputPath);
          } catch {
            initialCount = 0;
          }
        }

        const stageDefinitions = [
          {
            key: 'mol_prep',
            label: 'Mol Prep',
            candidates: [
              'stages/01_mol_prep/filtered_molecules.csv',
            ],
          },
          {
            key: 'descriptors',
            label: 'Descriptors',
            candidates: [
              'stages/02_descriptors_initial/filtered/filtered_molecules.csv',
              'Descriptors/passDescriptorsSMILES.csv',
            ],
          },
          {
            key: 'struct_filters',
            label: 'Struct Filters',
            candidates: [
              'stages/03_structural_filters_post/filtered_molecules.csv',
              'StructFilters/passStructFiltersSMILES.csv',
            ],
          },
          {
            key: 'synthesis',
            label: 'Synthesis',
            candidates: [
              'stages/04_synthesis/filtered_molecules.csv',
              'Synthesis/passSynthesisSMILES.csv',
            ],
          },
          {
            key: 'docking',
            label: 'Docking',
            candidates: [
              'stages/05_docking/ligands.csv',
            ],
          },
        ];

        const enabledStages = new Set(job.config.stages || []);
        const stats: Array<{ label: string; count: number | null; percent: number | null }> = [];

        for (const stage of stageDefinitions) {
          if (!enabledStages.has(stage.key)) continue;
          const path = await findExistingPath(basePath, stage.candidates);
          if (!path) {
            stats.push({ label: stage.label, count: null, percent: null });
            continue;
          }
          let count = 0;
          try {
            count = await countCsvRows(path);
          } catch {
            stats.push({ label: stage.label, count: null, percent: null });
            continue;
          }
          const percent = initialCount > 0 ? Math.round((count / initialCount) * 1000) / 10 : null;
          stats.push({ label: stage.label, count, percent });
        }

        if (!cancelled) {
          setStageStats(stats);
          setStageStatsError(null);
        }
      } catch (err) {
        if (!cancelled) {
          setStageStats(null);
          setStageStatsError(String(err));
        }
      }
    };

    loadStageStats();

    return () => {
      cancelled = true;
    };
  }, [job]);

  if (!job) {
    return (
      <AppShell
        padding={1}
        header={<Header title="Job Results" />}
        footer={<Footer />}
      >
        <Box marginY={2}>
          <Text color={theme.palette.error}>Job not found</Text>
        </Box>
      </AppShell>
    );
  }

  const statusColor = job.status === 'completed'
    ? theme.palette.success
    : job.status === 'running'
      ? theme.palette.warning
      : job.status === 'error'
        ? theme.palette.error
        : theme.palette.textMuted;

  const duration = job.endTime
    ? formatDuration((new Date(job.endTime).getTime() - new Date(job.startTime).getTime()) / 1000)
    : 'running...';

  return (
    <AppShell
      padding={1}
      header={<Header title={`Results: ${job.name || job.id}`} />}
      footer={<Footer />}
    >
      <Box flexDirection="column" marginY={1}>
        {/* Status */}
        <Box>
          <Text color={theme.palette.textMuted}>Status: </Text>
          <Text color={statusColor} bold>{job.status.toUpperCase()}</Text>
        </Box>

        {/* Timing */}
        <Box marginTop={1}>
          <Text color={theme.palette.accent} bold>Timing</Text>
        </Box>
        <Box marginLeft={2} flexDirection="column">
          <Text color={theme.palette.textMuted}>Started: {formatTimestamp(new Date(job.startTime))}</Text>
          {job.endTime && <Text color={theme.palette.textMuted}>Ended: {formatTimestamp(new Date(job.endTime))}</Text>}
          <Text color={theme.palette.textMuted}>Duration: {duration}</Text>
        </Box>

        {/* Config */}
        <Box marginTop={1}>
          <Text color={theme.palette.accent} bold>Configuration</Text>
        </Box>
        <Box marginLeft={2} flexDirection="column">
          <Text color={theme.palette.textMuted}>Input: {job.config.inputPath || '(not set)'}</Text>
          <Text color={theme.palette.textMuted}>Output: {job.config.outputPath || '(not set)'}</Text>
          <Text color={theme.palette.textMuted}>Stages: {job.config.stages.join(', ') || '(none)'}</Text>
        </Box>

        {/* Results */}
        {job.results && (
          <>
            <Box marginTop={1}>
              <Text color={theme.palette.accent} bold>Results</Text>
            </Box>
            <Box marginLeft={2} flexDirection="column">
              <Text color={theme.palette.textMuted}>Molecules processed: {job.results.moleculesProcessed}</Text>
              {job.results.moleculesFiltered !== undefined && (
                <Text color={theme.palette.textMuted}>Molecules filtered: {job.results.moleculesFiltered}</Text>
              )}
              {job.results.dockingHits !== undefined && (
                <Text color={theme.palette.textMuted}>Docking hits: {job.results.dockingHits}</Text>
              )}
            </Box>
          </>
        )}

        {/* Stage stats */}
        {stageStats && stageStats.length > 0 && (
          <>
            <Box marginTop={1}>
              <Text color={theme.palette.accent} bold>Stage Pass Counts</Text>
            </Box>
            <Box marginLeft={2} flexDirection="column">
              {stageStats.map((stat) => (
                <Text key={stat.label} color={theme.palette.textMuted}>
                  {stat.label}: {stat.count !== null ? stat.count : 'n/a'}
                  {stat.percent !== null ? ` (${stat.percent}%)` : ''}
                </Text>
              ))}
            </Box>
          </>
        )}

        {stageStatsError && (
          <Box marginTop={1}>
            <Text color={theme.palette.warning}>Stage stats unavailable: {stageStatsError}</Text>
          </Box>
        )}

        {/* Error */}
        {job.error && (
          <>
            <Box marginTop={1}>
              <Text color={theme.palette.error} bold>Error</Text>
            </Box>
            <Box marginLeft={2}>
              <Text color={theme.palette.error}>{job.error}</Text>
            </Box>
          </>
        )}
      </Box>
    </AppShell>
  );
}

export default Results;
