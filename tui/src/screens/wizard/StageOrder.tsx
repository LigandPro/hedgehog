import React, { useMemo } from 'react';
import { Box, Text, useInput } from 'ink';
import { Header } from '../../components/Header.js';
import { Footer } from '../../components/Footer.js';
import { useStore } from '../../store/index.js';

const STAGE_NAMES: Record<string, string> = {
  mol_prep: 'Mol Prep',
  descriptors: 'Descriptors',
  struct_filters: 'Struct Filters',
  synthesis: 'Synthesis',
  docking: 'Docking',
  docking_filters: 'Docking Filters',
};

export function StageOrder(): React.ReactElement {
  const setScreen = useStore((state) => state.setScreen);
  const wizard = useStore((state) => state.wizard);
  const getWizardSelectedStagesInOrder = useStore((state) => state.getWizardSelectedStagesInOrder);

  const selectedStagesInOrder = useMemo(() => getWizardSelectedStagesInOrder(), [wizard.stageOrder, wizard.selectedStages]);

  useInput((input, key) => {
    if (key.return || key.rightArrow) {
      goNext();
    } else if (key.escape || key.leftArrow || input === 'q') {
      setScreen('wizardStageSelection');
    }
  });

  const goNext = () => {
    // Navigate to first selected stage config (Mol Prep is not configured in the wizard UI).
    const firstStage = selectedStagesInOrder.find((s) => s !== 'mol_prep') || selectedStagesInOrder[0];
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

  const shortcuts = [
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

      <Footer shortcuts={shortcuts} />
    </Box>
  );
}

export default StageOrder;
