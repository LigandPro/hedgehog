import React, { useState } from 'react';
import { Box, Text, useInput } from 'ink';
import { Header } from '../../components/Header.js';
import { Footer } from '../../components/Footer.js';
import { useStore } from '../../store/index.js';
import type { Screen } from '../../types/index.js';

const STAGES = [
  { key: 'descriptors', name: 'Descriptors', description: 'Calculate molecular descriptors', configScreen: 'wizardConfigDescriptors' as Screen },
  { key: 'struct_filters', name: 'Struct Filters', description: 'Apply structural filters', configScreen: 'wizardConfigFilters' as Screen },
  { key: 'synthesis', name: 'Synthesis', description: 'Score synthesizability', configScreen: 'wizardConfigSynthesis' as Screen },
  { key: 'docking', name: 'Docking', description: 'Run molecular docking', configScreen: 'wizardConfigDocking' as Screen },
  { key: 'docking_filters', name: 'Docking Filters', description: 'Filter docking poses and interactions', configScreen: 'wizardConfigDockingFilters' as Screen },
];

export function StageSelection(): React.ReactElement {
  const setScreen = useStore((state) => state.setScreen);
  const wizard = useStore((state) => state.wizard);
  const toggleWizardStage = useStore((state) => state.toggleWizardStage);
  const showToast = useStore((state) => state.showToast);

  const [focusedIndex, setFocusedIndex] = useState(0);

  const selectedCount = wizard.selectedStages.length;

  const openConfig = () => {
    const stage = STAGES[focusedIndex];
    // Enable the stage if not already enabled
    if (!wizard.selectedStages.includes(stage.key)) {
      toggleWizardStage(stage.key);
    }
    setScreen(stage.configScreen);
  };

  const goNext = () => {
    if (selectedCount === 0) {
      showToast('warning', 'Select at least one stage');
      return;
    }
    setScreen('wizardReview');
  };

  useInput((input, key) => {
    if (key.upArrow) {
      setFocusedIndex(Math.max(0, focusedIndex - 1));
    } else if (key.downArrow) {
      setFocusedIndex(Math.min(STAGES.length - 1, focusedIndex + 1));
    } else if (input === ' ') {
      // Space - toggle stage
      toggleWizardStage(STAGES[focusedIndex].key);
    } else if (input === 'c') {
      // 'c' - open config
      openConfig();
    } else if (key.return || key.rightArrow) {
      // Enter or Right arrow - go to Review
      goNext();
    } else if (key.escape || key.leftArrow || input === 'q') {
      setScreen('wizardInputSelection');
    }
  });

  const shortcuts = [
    { key: 'Space', label: 'Toggle' },
    { key: 'c', label: 'Configure' },
    { key: 'Enter/→', label: 'Next' },
    { key: '←/Esc', label: 'Back' },
  ];

  return (
    <Box flexDirection="column" padding={1}>
      <Header title="Pipeline Wizard" subtitle="Step 2: Select Stages" />

      <Box flexDirection="column" marginY={1}>
        <Text color="cyan" bold>Select stages to include in pipeline:</Text>
      </Box>

      <Box flexDirection="column" marginY={1}>
        {STAGES.map((stage, index) => {
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

      <Footer shortcuts={shortcuts} />
    </Box>
  );
}

export default StageSelection;
