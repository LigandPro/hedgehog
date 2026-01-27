import React, { useState, useEffect } from 'react';
import { Box, Text, useInput } from 'ink';
import { Header } from '../components/Header.js';
import { Footer } from '../components/Footer.js';
import { SearchIndicator } from '../components/SearchIndicator.js';
import { useStore } from '../store/index.js';
import { getBridge } from '../services/python-bridge.js';
import { useSearch } from '../hooks/useSearch.js';
import { formatTimestamp, formatDuration } from '../utils/format.js';
import { getStatusIcon, getStatusColor } from '../utils/job-status.js';

export function History(): React.ReactElement {
  const setScreen = useStore((state) => state.setScreen);
  const jobHistory = useStore((state) => state.jobHistory);
  const setSelectedJob = useStore((state) => state.setSelectedJob);
  const removeJobFromHistory = useStore((state) => state.removeJobFromHistory);
  const showConfirm = useStore((state) => state.showConfirm);
  const showToast = useStore((state) => state.showToast);
  const setSearchActive = useStore((state) => state.setSearchActive);
  const setSearchQuery = useStore((state) => state.setSearchQuery);

  const [selectedIndex, setSelectedIndex] = useState(0);

  // Search functionality
  const {
    searchQuery,
    searchActive,
    filteredItems: filteredJobs,
    handleSearchInput,
    highlightMatch,
  } = useSearch({
    items: jobHistory,
    searchFields: (job) => [job.name || job.id, job.status, job.config?.inputPath || ''],
  });

  // Sync search state with global store for Footer display
  useEffect(() => {
    setSearchActive(searchActive);
    setSearchQuery(searchQuery);
    return () => {
      setSearchActive(false);
      setSearchQuery('');
    };
  }, [searchActive, searchQuery, setSearchActive, setSearchQuery]);

  // Reset selection when search changes
  useEffect(() => {
    setSelectedIndex(0);
  }, [searchQuery]);

  const handleDelete = (jobId: string, jobName: string) => {
    showConfirm({
      title: 'Delete Job?',
      message: `Permanently remove "${jobName}" from history?`,
      confirmLabel: 'Delete',
      cancelLabel: 'Keep',
      onConfirm: () => {
        removeJobFromHistory(jobId);
        getBridge().deleteJob(jobId).catch(() => {});
        if (selectedIndex >= filteredJobs.length - 1) {
          setSelectedIndex(Math.max(0, filteredJobs.length - 2));
        }
      },
    });
  };

  useInput((input, key) => {
    // Handle search input first
    if (handleSearchInput(input, key)) {
      return;
    }

    if (key.upArrow) {
      setSelectedIndex(Math.max(0, selectedIndex - 1));
    } else if (key.downArrow) {
      setSelectedIndex(Math.min(filteredJobs.length - 1, selectedIndex + 1));
    } else if ((key.return || key.rightArrow) && filteredJobs[selectedIndex]) {
      setSelectedJob(filteredJobs[selectedIndex].id);
      setScreen('results');
    } else if (input === 'd' && filteredJobs[selectedIndex]) {
      const job = filteredJobs[selectedIndex];
      handleDelete(job.id, job.name || `Job ${job.id}`);
    } else if (key.escape || key.leftArrow) {
      setScreen('welcome');
    }
  });

  return (
    <Box flexDirection="column" padding={1}>
      <Header title="Job History" />

      <SearchIndicator active={searchActive} query={searchQuery} />

      {filteredJobs.length === 0 ? (
        <Box marginY={2}>
          {jobHistory.length === 0 ? (
            <Text dimColor>No job history yet. Run a pipeline to see it here.</Text>
          ) : (
            <Text dimColor>No jobs match "{searchQuery}"</Text>
          )}
        </Box>
      ) : (
        <Box flexDirection="column" marginY={1}>
          {filteredJobs.map((job, index) => {
            const isSelected = selectedIndex === index;
            const statusIcon = getStatusIcon(job.status);
            const statusColor = getStatusColor(job.status);

            const duration = job.endTime
              ? formatDuration((new Date(job.endTime).getTime() - new Date(job.startTime).getTime()) / 1000)
              : 'running...';

            const jobName = job.name || `Job ${job.id}`;
            const nameParts = searchQuery ? highlightMatch(jobName) : [{ text: jobName, highlighted: false }];

            return (
              <Box key={job.id} flexDirection="column">
                <Box>
                  <Text color={isSelected ? 'cyan' : 'gray'}>{isSelected ? '▶ ' : '  '}</Text>
                  <Text color={statusColor}>{statusIcon}</Text>
                  <Text> </Text>
                  {nameParts.map((part, pi) => (
                    <Text
                      key={pi}
                      color={isSelected ? (part.highlighted ? 'yellow' : 'white') : (part.highlighted ? 'yellow' : 'gray')}
                      bold={part.highlighted}
                    >
                      {part.text}
                    </Text>
                  ))}
                </Box>
                {isSelected && (
                  <Box marginLeft={4} flexDirection="column">
                    <Text dimColor>Started: {formatTimestamp(new Date(job.startTime))}</Text>
                    <Text dimColor>Duration: {duration}</Text>
                    {job.results && (
                      <Text dimColor>Processed: {job.results.moleculesProcessed} molecules</Text>
                    )}
                  </Box>
                )}
              </Box>
            );
          })}
        </Box>
      )}

      {/* Stats */}
      {filteredJobs.length > 0 && (
        <Box>
          <Text dimColor>
            {filteredJobs.length} job{filteredJobs.length !== 1 ? 's' : ''}
            {searchQuery && ` (filtered from ${jobHistory.length})`}
          </Text>
        </Box>
      )}

      <Footer />
    </Box>
  );
}

export default History;
