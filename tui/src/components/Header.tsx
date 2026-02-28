import React, { memo, useEffect, useState } from 'react';
import { Box, Text } from 'ink';
import { useStore } from '../store/index.js';
import { useTerminalSize } from '../hooks/useTerminalSize.js';
import { useTheme } from '../theme/context.js';
import { selectAsciiLogo } from '../constants/ascii-logos.js';

interface HeaderProps {
  title?: string;
  subtitle?: string;
  showLogo?: boolean;
}

function fitLine(value: string, width: number): string {
  if (width <= 0) return '';
  if (value.length >= width) {
    return value.slice(0, width);
  }
  return value.padEnd(width, ' ');
}

const STAGE_NAMES: Record<string, string> = {
  mol_prep: 'Mol Prep',
  descriptors: 'Descriptors',
  struct_filters: 'Filters',
  synthesis: 'Synthesis',
  docking: 'Docking',
  docking_filters: 'Docking Filters',
};

const SPINNER_FRAMES = ['⠋', '⠙', '⠹', '⠸'];

// Memoized running indicator to prevent unnecessary re-renders
// NOTE: Elapsed time is intentionally NOT shown here to prevent 1-second re-renders
// on all screens. Elapsed time is only shown in PipelineRunner screen.
const RunningIndicator = memo(function RunningIndicator(): React.ReactElement | null {
  const { theme } = useTheme();
  const isRunning = useStore((state) => state.isRunning);
  const stages = useStore((state) => state.stages);
  const selectedStages = useStore((state) =>
    state.wizard.stageOrder.filter((stage) => state.wizard.selectedStages.includes(stage))
  );
  const currentStage = useStore((state) => state.pipelineProgress.currentStage);
  const stageIndex = useStore((state) => state.pipelineProgress.stageIndex);
  const totalStages = useStore((state) => state.pipelineProgress.totalStages);
  const stageProgress = useStore((state) => state.pipelineProgress.stageProgress);
  const [frameIndex, setFrameIndex] = useState(0);

  useEffect(() => {
    if (!isRunning) {
      setFrameIndex(0);
      return;
    }

    const timer = setInterval(() => {
      setFrameIndex((idx) => (idx + 1) % SPINNER_FRAMES.length);
    }, 500);

    return () => clearInterval(timer);
  }, [isRunning]);

  if (!isRunning) return null;

  const completedVisibleStages = selectedStages.filter(
    (stage) => stages[stage]?.status === 'completed'
  ).length;
  const isFinalizing = selectedStages.length > 0 && completedVisibleStages === selectedStages.length;

  const stageName = STAGE_NAMES[currentStage] || currentStage;
  const progressText = totalStages > 0
    ? `${stageIndex}/${totalStages}`
    : '';
  const percentText = stageProgress > 0
    ? ` (${stageProgress}%)`
    : '';

  return (
    <Box marginLeft={1}>
      <Text color={theme.palette.textMuted}> | </Text>
      <Text color={theme.palette.warning}>{SPINNER_FRAMES[frameIndex]}</Text>
      <Text color={theme.palette.warning} bold> {isFinalizing ? 'Finalizing' : 'Running'}</Text>
      {isFinalizing ? (
        <>
          <Text color={theme.palette.textMuted}>: </Text>
          <Text color={theme.palette.text}>Writing outputs and report</Text>
        </>
      ) : stageName && (
        <>
          <Text color={theme.palette.textMuted}>: </Text>
          <Text color={theme.palette.text}>{progressText} {stageName}{percentText}</Text>
        </>
      )}
    </Box>
  );
});

export const Header = memo(function Header({ title, subtitle, showLogo = false }: HeaderProps): React.ReactElement {
  const { width: terminalWidth } = useTerminalSize();
  const { theme } = useTheme();
  const contentWidth = Math.max(8, terminalWidth - 2);
  const titleLine = `HEDGEHOG PIPELINE ENGINE${title ? ` | ${title}` : ''}`;

  // Select the largest logo that fits, or none if too narrow
  const selectedLogo = showLogo
    ? selectAsciiLogo(contentWidth)
    : null;

  return (
    <Box flexDirection="column" marginBottom={1}>
      {selectedLogo && (
        <Box flexDirection="column" marginBottom={1}>
          {selectedLogo.logo.split('\n').map((line, index) => (
            <Text
              key={`logo-${index}`}
              color="white"
            >
              {fitLine(line, contentWidth)}
            </Text>
          ))}
        </Box>
      )}
      <Box>
        {/* Keep the header title readable and consistent across themes. */}
        <Text color={theme.palette.text} bold>{fitLine(titleLine, contentWidth)}</Text>
      </Box>
      <RunningIndicator />
      {subtitle && (
        <Box marginTop={0}>
          <Text color={theme.palette.textMuted}>{fitLine(subtitle, contentWidth)}</Text>
        </Box>
      )}
      <Box marginTop={0}>
        <Text color={theme.palette.border}>{'─'.repeat(contentWidth)}</Text>
      </Box>
    </Box>
  );
});

export default Header;
