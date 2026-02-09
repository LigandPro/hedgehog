import React, { useMemo } from 'react';
import { Box, Text, useInput } from 'ink';
import { Header } from '../../components/Header.js';
import { Footer } from '../../components/Footer.js';
import { useStore } from '../../store/index.js';
import { useTerminalSize } from '../../hooks/useTerminalSize.js';

const STAGE_NAMES: Record<string, string> = {
  descriptors: 'Descriptors',
  struct_filters: 'Struct Filters',
  synthesis: 'Synthesis',
  docking: 'Docking',
  docking_filters: 'Docking Filters',
};

export function StageOrder(): React.ReactElement {
  const setScreen = useStore((state) => state.setScreen);
  const wizard = useStore((state) => state.wizard);
  const setWizardDependency = useStore((state) => state.setWizardDependency);
  const getWizardSelectedStagesInOrder = useStore((state) => state.getWizardSelectedStagesInOrder);

  const { width: terminalWidth } = useTerminalSize();

  const selectedStagesInOrder = useMemo(() => getWizardSelectedStagesInOrder(), [wizard.stageOrder, wizard.selectedStages]);

  // Check if both descriptors and struct_filters are selected (for dependency option)
  const showDependencyOption = wizard.selectedStages.includes('descriptors') && wizard.selectedStages.includes('struct_filters');

  useInput((input, key) => {
    if (showDependencyOption) {
      if (input === ' ') {
        setWizardDependency('runFiltersBeforeDescriptors', !wizard.dependencies.runFiltersBeforeDescriptors);
      }
    }

    if (key.return || key.rightArrow) {
      goNext();
    } else if (key.escape || key.leftArrow || input === 'q') {
      setScreen('wizardStageSelection');
    }
  });

  const goNext = () => {
    // Navigate to first selected stage config
    const firstStage = selectedStagesInOrder[0];
    if (firstStage === 'descriptors') {
      setScreen('wizardConfigDescriptors');
    } else if (firstStage === 'struct_filters') {
      setScreen('wizardConfigFilters');
    } else if (firstStage === 'synthesis') {
      setScreen('wizardConfigSynthesis');
    } else if (firstStage === 'docking') {
      setScreen('wizardConfigDocking');
    } else if (firstStage === 'docking_filters') {
      setScreen('wizardConfigDockingFilters');
    } else {
      setScreen('wizardReview');
    }
  };

  const totalSteps = 2 + selectedStagesInOrder.length + 1; // selection + order + configs + review

  const shortcuts = showDependencyOption
    ? [
        { key: 'Space', label: 'Toggle dependency' },
        { key: '→/Enter', label: 'Next' },
        { key: '←/Esc', label: 'Back' },
      ]
    : [
        { key: '→/Enter', label: 'Next' },
        { key: '←/Esc', label: 'Back' },
      ];

  return (
    <Box flexDirection="column" padding={1}>
      <Header title="Pipeline Wizard" subtitle={`Stage Order (2/${totalSteps})`} />

      <Box flexDirection="column" marginY={1}>
        <Text color="cyan" bold>Selected stages will run in this order:</Text>
      </Box>

      <Box flexDirection="column" marginY={1}>
        {selectedStagesInOrder.map((stage, index) => (
          <Box key={stage}>
            <Text dimColor>{`${index + 1}. `}</Text>
            <Text color="white">{STAGE_NAMES[stage] || stage}</Text>
          </Box>
        ))}
      </Box>

      {showDependencyOption && (
        <>
          <Box marginY={1}>
            <Text color="gray">{'─'.repeat(terminalWidth - 2)}</Text>
          </Box>

          <Box flexDirection="column">
            <Text color="cyan" bold>Dependencies:</Text>
            <Box marginTop={1}>
              <Text color="cyan">{'▸ '}</Text>
              <Text color={wizard.dependencies.runFiltersBeforeDescriptors ? 'green' : 'gray'}>
                [{wizard.dependencies.runFiltersBeforeDescriptors ? '✓' : ' '}]
              </Text>
              <Text color="white">
                {' '}Run Struct Filters before Descriptors
              </Text>
            </Box>
          </Box>
        </>
      )}

      <Footer shortcuts={shortcuts} />
    </Box>
  );
}

export default StageOrder;
