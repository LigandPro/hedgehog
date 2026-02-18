import React, { useMemo, useState } from 'react';
import { Box, Text, useInput } from 'ink';
import { Header } from '../../components/Header.js';
import { Footer } from '../../components/Footer.js';
import { useStore } from '../../store/index.js';
import { runWizardPipeline } from './runPipeline.js';
import { STAGE_METADATA, WIZARD_STAGE_ORDER } from './stageMetadata.js';

export function StageSelection(): React.ReactElement {
  const setScreen = useStore((state) => state.setScreen);
  const wizard = useStore((state) => state.wizard);
  const toggleWizardStage = useStore((state) => state.toggleWizardStage);
  const showToast = useStore((state) => state.showToast);

  const [focusedIndex, setFocusedIndex] = useState(0);
  const [starting, setStarting] = useState(false);

  const stages = useMemo(
    () => WIZARD_STAGE_ORDER.map((key) => ({
      key,
      name: STAGE_METADATA[key].title,
      description: STAGE_METADATA[key].shortDescription,
      configScreen: STAGE_METADATA[key].configScreen,
    })),
    []
  );

  const selectedCount = wizard.selectedStages.length;

  const openConfig = () => {
    const stage = stages[focusedIndex];
    if (!stage?.configScreen) {
      showToast('warning', `No config screen for stage: ${stage?.name || 'unknown'}`);
      return;
    }

    if (!wizard.selectedStages.includes(stage.key)) {
      toggleWizardStage(stage.key);
    }
    setScreen(stage.configScreen);
  };

  const openReview = () => {
    if (selectedCount === 0) {
      showToast('warning', 'Select at least one stage');
      return;
    }
    setScreen('wizardReview');
  };

  const startFast = async () => {
    if (starting) return;
    if (selectedCount === 0) {
      showToast('warning', 'Select at least one stage');
      return;
    }

    setStarting(true);
    try {
      await runWizardPipeline({ onPreflightErrorScreen: 'wizardReview' });
    } finally {
      setStarting(false);
    }
  };

  useInput((input, key) => {
    if (starting) return;

    if (key.upArrow) {
      setFocusedIndex(Math.max(0, focusedIndex - 1));
    } else if (key.downArrow) {
      setFocusedIndex(Math.min(stages.length - 1, focusedIndex + 1));
    } else if (input === ' ') {
      if (stages[focusedIndex].key === 'mol_prep') {
        showToast('info', 'Mol Prep is a prerequisite stage and is usually kept enabled');
      } else {
        toggleWizardStage(stages[focusedIndex].key);
      }
    } else if (input === 'c') {
      openConfig();
    } else if (input === 'r' || key.rightArrow) {
      openReview();
    } else if (key.return) {
      void startFast();
    } else if (key.escape || key.leftArrow || input === 'q') {
      setScreen('wizardInputSelection');
    }
  });

  const shortcuts = [
    { key: 'Space', label: 'Toggle' },
    { key: 'c', label: 'Configure' },
    { key: 'Enter', label: 'Fast start' },
    { key: 'r/→', label: 'Detailed review' },
    { key: '←/Esc', label: 'Back' },
  ];

  return (
    <Box flexDirection="column" padding={1}>
      <Header title="Pipeline Wizard" subtitle="Step 2: Select Stages" />

      <Box flexDirection="column" marginY={1}>
        <Text color="cyan" bold>Select stages to include in pipeline:</Text>
      </Box>

      <Box flexDirection="column" marginY={1}>
        {stages.map((stage, index) => {
          const isSelected = wizard.selectedStages.includes(stage.key);
          const isFocused = focusedIndex === index;

          return (
            <Box key={stage.key}>
              <Text color={isFocused ? 'cyan' : 'white'}>
                {isFocused ? '▸ ' : '  '}
              </Text>
              <Text color={isSelected ? 'green' : 'gray'}>
                [{isSelected ? '✓' : ' '}]
              </Text>
              <Text color={isFocused ? 'white' : 'gray'}> {stage.name.padEnd(16)}</Text>
              <Text dimColor>{stage.description}</Text>
            </Box>
          );
        })}
      </Box>

      <Box marginY={1}>
        <Text color="cyan">{selectedCount}</Text>
        <Text dimColor> stage{selectedCount !== 1 ? 's' : ''} selected</Text>
      </Box>

      {starting && (
        <Box marginBottom={1}>
          <Text color="yellow">Running preflight and starting pipeline...</Text>
        </Box>
      )}

      <Footer shortcuts={shortcuts} />
    </Box>
  );
}

export default StageSelection;
