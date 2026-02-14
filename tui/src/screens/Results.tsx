import React, { useEffect, useState } from 'react';
import { Box, Text, useInput } from 'ink';
import { Header } from '../components/Header.js';
import { Footer } from '../components/Footer.js';
import { useStore } from '../store/index.js';
import { formatTimestamp, formatDuration } from '../utils/format.js';
import { getStatusColor } from '../utils/job-status.js';
import { createReadStream, existsSync } from 'fs';
import { access } from 'fs/promises';
import { createInterface } from 'readline';
import { dirname, isAbsolute, join, resolve } from 'path';

export function Results(): React.ReactElement {
  const setScreen = useStore((state) => state.setScreen);
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
      setScreen('history');
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
              'stages/00_mol_prep/filtered_molecules.csv',
            ],
          },
          {
            key: 'descriptors',
            label: 'Descriptors',
            candidates: [
              'stages/01_descriptors_initial/filtered/filtered_molecules.csv',
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
      <Box flexDirection="column" padding={1}>
        <Header title="Job Results" />
        <Box marginY={2}>
          <Text color="red">Job not found</Text>
        </Box>
        <Footer
          shortcuts={[{ key: 'Esc', label: 'Back' }]}
          overrideShortcuts
        />
      </Box>
    );
  }

  const statusColor = getStatusColor(job.status);

  const duration = job.endTime
    ? formatDuration((new Date(job.endTime).getTime() - new Date(job.startTime).getTime()) / 1000)
    : 'running...';

  return (
    <Box flexDirection="column" padding={1}>
      <Header title={`Results: ${job.name || job.id}`} />

      <Box flexDirection="column" marginY={1}>
        {/* Status */}
        <Box>
          <Text dimColor>Status: </Text>
          <Text color={statusColor} bold>{job.status.toUpperCase()}</Text>
        </Box>

        {/* Timing */}
        <Box marginTop={1}>
          <Text color="cyan" bold>Timing</Text>
        </Box>
        <Box marginLeft={2} flexDirection="column">
          <Text dimColor>Started: {formatTimestamp(new Date(job.startTime))}</Text>
          {job.endTime && <Text dimColor>Ended: {formatTimestamp(new Date(job.endTime))}</Text>}
          <Text dimColor>Duration: {duration}</Text>
        </Box>

        {/* Config */}
        <Box marginTop={1}>
          <Text color="cyan" bold>Configuration</Text>
        </Box>
        <Box marginLeft={2} flexDirection="column">
          <Text dimColor>Input: {job.config.inputPath || '(not set)'}</Text>
          <Text dimColor>Output: {job.config.outputPath || '(not set)'}</Text>
          <Text dimColor>Stages: {job.config.stages.join(', ') || '(none)'}</Text>
        </Box>

        {/* Results */}
        {job.results && (
          <>
            <Box marginTop={1}>
              <Text color="cyan" bold>Results</Text>
            </Box>
            <Box marginLeft={2} flexDirection="column">
              <Text dimColor>Molecules processed: {job.results.moleculesProcessed}</Text>
              {job.results.moleculesFiltered !== undefined && (
                <Text dimColor>Molecules filtered: {job.results.moleculesFiltered}</Text>
              )}
              {job.results.dockingHits !== undefined && (
                <Text dimColor>Docking hits: {job.results.dockingHits}</Text>
              )}
            </Box>
          </>
        )}

        {/* Stage stats */}
        {stageStats && stageStats.length > 0 && (
          <>
            <Box marginTop={1}>
              <Text color="cyan" bold>Stage Pass Counts</Text>
            </Box>
            <Box marginLeft={2} flexDirection="column">
              {stageStats.map((stat) => (
                <Text key={stat.label} dimColor>
                  {stat.label}: {stat.count !== null ? stat.count : 'n/a'}
                  {stat.percent !== null ? ` (${stat.percent}%)` : ''}
                </Text>
              ))}
            </Box>
          </>
        )}

        {stageStatsError && (
          <Box marginTop={1}>
            <Text color="yellow">Stage stats unavailable: {stageStatsError}</Text>
          </Box>
        )}

        {/* Error */}
        {job.error && (
          <>
            <Box marginTop={1}>
              <Text color="red" bold>Error</Text>
            </Box>
            <Box marginLeft={2}>
              <Text color="red">{job.error}</Text>
            </Box>
          </>
        )}
      </Box>

      <Footer />
    </Box>
  );
}

export default Results;
