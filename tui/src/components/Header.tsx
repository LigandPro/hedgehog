import React, { memo, useEffect, useState } from 'react';
import { Box, Text } from 'ink';
import { useStore } from '../store/index.js';
import { useTerminalSize } from '../hooks/useTerminalSize.js';

interface HeaderProps {
  title?: string;
  subtitle?: string;
  showLogo?: boolean;
}

// Progressive ASCII logos - from full HEDGEHOG down to HEDGE
// Each has a width property for easy selection based on terminal width
const ASCII_LOGOS = [
  {
    width: 68,
    logo: `██╗  ██╗███████╗██████╗  ██████╗ ███████╗██╗  ██╗ ██████╗  ██████╗
██║  ██║██╔════╝██╔══██╗██╔════╝ ██╔════╝██║  ██║██╔═══██╗██╔════╝
███████║█████╗  ██║  ██║██║  ███╗█████╗  ███████║██║   ██║██║  ███╗
██╔══██║██╔══╝  ██║  ██║██║   ██║██╔══╝  ██╔══██║██║   ██║██║   ██║
██║  ██║███████╗██████╔╝╚██████╔╝███████╗██║  ██║╚██████╔╝╚██████╔╝
╚═╝  ╚═╝╚══════╝╚═════╝  ╚═════╝ ╚══════╝╚═╝  ╚═╝ ╚═════╝  ╚═════╝`,
  },
  {
    width: 59,
    logo: `██╗  ██╗███████╗██████╗  ██████╗ ███████╗██╗  ██╗ ██████╗
██║  ██║██╔════╝██╔══██╗██╔════╝ ██╔════╝██║  ██║██╔═══██╗
███████║█████╗  ██║  ██║██║  ███╗█████╗  ███████║██║   ██║
██╔══██║██╔══╝  ██║  ██║██║   ██║██╔══╝  ██╔══██║██║   ██║
██║  ██║███████╗██████╔╝╚██████╔╝███████╗██║  ██║╚██████╔╝
╚═╝  ╚═╝╚══════╝╚═════╝  ╚═════╝ ╚══════╝╚═╝  ╚═╝ ╚═════╝`,
  },
  {
    width: 49,
    logo: `██╗  ██╗███████╗██████╗  ██████╗ ███████╗██╗  ██╗
██║  ██║██╔════╝██╔══██╗██╔════╝ ██╔════╝██║  ██║
███████║█████╗  ██║  ██║██║  ███╗█████╗  ███████║
██╔══██║██╔══╝  ██║  ██║██║   ██║██╔══╝  ██╔══██║
██║  ██║███████╗██████╔╝╚██████╔╝███████╗██║  ██║
╚═╝  ╚═╝╚══════╝╚═════╝  ╚═════╝ ╚══════╝╚═╝  ╚═╝`,
  },
  {
    width: 40,
    logo: `██╗  ██╗███████╗██████╗  ██████╗ ███████╗
██║  ██║██╔════╝██╔══██╗██╔════╝ ██╔════╝
███████║█████╗  ██║  ██║██║  ███╗█████╗
██╔══██║██╔══╝  ██║  ██║██║   ██║██╔══╝
██║  ██║███████╗██████╔╝╚██████╔╝███████╗
╚═╝  ╚═╝╚══════╝╚═════╝  ╚═════╝ ╚══════╝`,
  },
];

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
  const isRunning = useStore((state) => state.isRunning);
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

  const stageName = STAGE_NAMES[currentStage] || currentStage;
  const progressText = totalStages > 0
    ? `${stageIndex}/${totalStages}`
    : '';
  const percentText = stageProgress > 0
    ? ` (${stageProgress}%)`
    : '';

  return (
    <Box marginLeft={1}>
      <Text color="gray"> | </Text>
      <Text color="yellow">{SPINNER_FRAMES[frameIndex]}</Text>
      <Text color="yellow" bold> Running</Text>
      {stageName && (
        <>
          <Text color="gray">: </Text>
          <Text color="white">{progressText} {stageName}{percentText}</Text>
        </>
      )}
    </Box>
  );
});

export const Header = memo(function Header({ title, subtitle, showLogo = true }: HeaderProps): React.ReactElement {
  const { width: terminalWidth } = useTerminalSize();

  // Select the largest logo that fits, or none if too narrow
  const selectedLogo = showLogo
    ? ASCII_LOGOS.find((l) => terminalWidth >= l.width)
    : null;

  return (
    <Box flexDirection="column" marginBottom={1}>
      {selectedLogo && (
        <Box flexDirection="column" marginBottom={1}>
          <Text color="cyan">{selectedLogo.logo}</Text>
        </Box>
      )}
      <Box>
        <Text color="cyan" bold>HEDGEHOG PIPELINE ENGINE</Text>
        <Text color="gray"> v0.1.0</Text>
        {title && (
          <>
            <Text color="gray"> | </Text>
            <Text color="white">{title}</Text>
          </>
        )}
        <RunningIndicator />
      </Box>
      {subtitle && (
        <Box marginTop={0}>
          <Text dimColor>{subtitle}</Text>
        </Box>
      )}
      <Box marginTop={0}>
        <Text color="gray">{'─'.repeat(terminalWidth - 2)}</Text>
      </Box>
    </Box>
  );
});

export default Header;
