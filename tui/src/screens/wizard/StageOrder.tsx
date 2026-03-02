import React, { useMemo } from 'react';
import { Box, Text, useInput } from 'ink';
import { Header } from '../../components/Header.js';
import { Footer } from '../../components/Footer.js';
import { AppShell } from '../../components/AppShell.js';
import { useStore } from '../../store/index.js';
import { useTheme } from '../../theme/context.js';

const STAGE_NAMES: Record<string, string> = {
  mol_prep: 'Mol Prep',
  descriptors: 'Descriptors',
  struct_filters: 'Struct Filters',
  synthesis: 'Synthesis',
  docking: 'Docking',
  docking_filters: 'Docking Filters',
};

export function StageOrder(): React.ReactElement {
  const { theme } = useTheme();
  const setScreen = useStore((state) => state.setScreen);
  const goBack = useStore((state) => state.goBack);
  const wizard = useStore((state) => state.wizard);
  const getWizardSelectedStagesInOrder = useStore((state) => state.getWizardSelectedStagesInOrder);

  const selectedStagesInOrder = useMemo(() => getWizardSelectedStagesInOrder(), [wizard.stageOrder, wizard.selectedStages]);

  useInput((input, key) => {
    if (key.return || key.rightArrow) {
      goNext();
    } else if (key.escape || key.leftArrow) {
      goBack();
    }
  });

  const goNext = () => {
    // Navigate to first selected stage config.
    const firstStage = selectedStagesInOrder[0];
    if (firstStage === 'mol_prep') {
      setScreen('wizardConfigMolPrep');
    } else if (firstStage === 'descriptors') {
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

  return (
    <AppShell
      padding={1}
      header={<Header title="Pipeline Wizard" subtitle={`Stage Order (2/${totalSteps})`} />}
      footer={<Footer />}
    >

      <Box flexDirection="column" marginY={1}>
        <Text color={theme.palette.accent} bold>Selected stages will run in this order:</Text>
      </Box>

      <Box flexDirection="column" marginY={1}>
        {selectedStagesInOrder.map((stage, index) => (
          <Box key={stage}>
            <Text color={theme.palette.textMuted}>{`${index + 1}. `}</Text>
            <Text color={theme.palette.text}>{STAGE_NAMES[stage] || stage}</Text>
          </Box>
        ))}
      </Box>
    </AppShell>
  );
}

export default StageOrder;
