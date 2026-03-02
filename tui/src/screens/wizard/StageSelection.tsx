import React, { useMemo, useState } from 'react';
import { Box, Text, useInput } from 'ink';
import { Header } from '../../components/Header.js';
import { Footer } from '../../components/Footer.js';
import { AppShell } from '../../components/AppShell.js';
import { useStore } from '../../store/index.js';
import { useTheme } from '../../theme/context.js';
import { useTerminalSize } from '../../hooks/useTerminalSize.js';
import { STAGE_METADATA, WIZARD_STAGE_ORDER } from './stageMetadata.js';

export function StageSelection(): React.ReactElement {
  const { theme } = useTheme();
  const { width: terminalWidth } = useTerminalSize();
  const setScreen = useStore((state) => state.setScreen);
  const goBack = useStore((state) => state.goBack);
  const wizard = useStore((state) => state.wizard);
  const toggleWizardStage = useStore((state) => state.toggleWizardStage);
  const showToast = useStore((state) => state.showToast);

  const [focusedIndex, setFocusedIndex] = useState(0);

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
  const rowWidth = Math.max(1, terminalWidth - 4);
  const nameWidth = Math.min(16, Math.max(8, rowWidth - 8));
  const descriptionWidth = Math.max(1, rowWidth - nameWidth - 8);

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

  useInput((input, key) => {
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
      openReview();
    } else if (key.escape || key.leftArrow) {
      goBack();
    }
  });

  return (
    <AppShell
      padding={1}
      header={<Header title="Pipeline Wizard" subtitle="Step 2: Select Stages" />}
      footer={<Footer />}
    >
      <Box flexDirection="column" marginY={1}>
        <Text color={theme.palette.accent} bold>Select stages to include in pipeline:</Text>
      </Box>

      <Box flexDirection="column" marginY={1}>
        {stages.map((stage, index) => {
          const isSelected = wizard.selectedStages.includes(stage.key);
          const isFocused = focusedIndex === index;

          return (
            <Box key={stage.key} width={rowWidth}>
              <Text
                color={isFocused ? theme.palette.primary : theme.palette.text}
                backgroundColor={isFocused ? theme.palette.panelStrong : undefined}
              >
                {isFocused ? '▸ ' : '  '}
              </Text>
              <Text
                color={isSelected ? theme.palette.success : theme.palette.textMuted}
                backgroundColor={isFocused ? theme.palette.panelStrong : undefined}
              >
                [{isSelected ? '✓' : ' '}]
              </Text>
              <Text backgroundColor={isFocused ? theme.palette.panelStrong : undefined}> </Text>
              <Box width={nameWidth}>
                <Text
                  color={isFocused ? theme.palette.text : theme.palette.textMuted}
                  backgroundColor={isFocused ? theme.palette.panelStrong : undefined}
                >
                  {stage.name}
                </Text>
              </Box>
              <Box width={descriptionWidth}>
                <Text
                  color={theme.palette.textMuted}
                  backgroundColor={isFocused ? theme.palette.panelStrong : undefined}
                  wrap="truncate-end"
                >
                  {stage.description}
                </Text>
              </Box>
            </Box>
          );
        })}
      </Box>

      <Box marginY={1}>
        <Text color={theme.palette.primary}>{selectedCount}</Text>
        <Text color={theme.palette.textMuted}> stage{selectedCount !== 1 ? 's' : ''} selected</Text>
      </Box>
    </AppShell>
  );
}

export default StageSelection;
