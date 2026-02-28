import React, { useState, useEffect } from 'react';
import { Box, Text, useInput } from 'ink';
import TextInput from 'ink-text-input';
import { Header } from '../components/Header.js';
import { Footer } from '../components/Footer.js';
import { Spinner } from '../components/Spinner.js';
import { useStore } from '../store/index.js';
import { getBridge } from '../services/python-bridge.js';
import { useInputLock } from '../hooks/useInputLock.js';
import { useTheme } from '../theme/context.js';
import { useTerminalSize } from '../hooks/useTerminalSize.js';
import type { SynthesisConfig } from '../types/index.js';

interface FormField {
  key: keyof SynthesisConfig;
  label: string;
  type: 'number' | 'boolean';
  description?: string;
  section?: string;
}

function truncateMiddle(value: string, width: number): string {
  if (width <= 0) return '';
  if (value.length <= width) return value;
  if (width <= 1) return '…';
  if (width <= 3) return `${value.slice(0, width - 1)}…`;
  const left = Math.ceil((width - 1) / 2);
  const right = Math.floor((width - 1) / 2);
  return `${value.slice(0, left)}…${value.slice(value.length - right)}`;
}

const fields: FormField[] = [
  { key: 'run', label: 'Run Stage', type: 'boolean', section: 'General' },
  { key: 'run_retrosynthesis', label: 'Run Retrosynthesis', type: 'boolean', description: 'Enable AiZynthFinder retrosynthesis route search (press [r] to configure)' },
  { key: 'filter_solved_only', label: 'Filter Solved Only', type: 'boolean', description: 'Only keep molecules with found retrosynthesis path' },
  { key: 'sa_score_min', label: 'SA Score Min', type: 'number', section: 'Score Thresholds', description: 'Synthetic Accessibility (1-10, lower = easier)' },
  { key: 'sa_score_max', label: 'SA Score Max', type: 'number' },
  { key: 'syba_score_min', label: 'SYBA Score Min', type: 'number', description: 'Synthetic Bayesian Accessibility' },
  { key: 'syba_score_max', label: 'SYBA Score Max', type: 'number' },
  { key: 'ra_score_min', label: 'RA Score Min', type: 'number', description: 'Retrosynthesis probability (0-1)' },
  { key: 'ra_score_max', label: 'RA Score Max', type: 'number' },
];

export function ConfigSynthesis(): React.ReactElement {
  const { theme } = useTheme();
  const { width: terminalWidth } = useTerminalSize();
  const setScreen = useStore((state) => state.setScreen);
  const goBack = useStore((state) => state.goBack);
  const config = useStore((state) => state.configs.synthesis);
  const setConfig = useStore((state) => state.setConfig);
  const isBackendReady = useStore((state) => state.isBackendReady);
  const showToast = useStore((state) => state.showToast);
  const showConfirm = useStore((state) => state.showConfirm);

  const [loading, setLoading] = useState(!config);
  const [selectedIndex, setSelectedIndex] = useState(0);
  const [editMode, setEditMode] = useState(false);
  const [editValue, setEditValue] = useState('');
  const [values, setValues] = useState<Partial<SynthesisConfig>>(config || {});
  const [error, setError] = useState<string | null>(null);
  const [isDirty, setIsDirty] = useState(false);
  useInputLock(editMode);
  const selectedField = fields[selectedIndex];
  const selectedDescription = selectedField?.description || '';

  useEffect(() => {
    if (!config && isBackendReady) {
      loadConfig();
    } else if (config) {
      setValues(config);
      setLoading(false);
    }
  }, [config, isBackendReady]);

  const loadConfig = async () => {
    try {
      const bridge = getBridge();
      const data = await bridge.loadConfig('synthesis') as unknown as SynthesisConfig;
      setConfig('synthesis', data);
      setValues(data);
    } catch (err) {
      setError(String(err));
    } finally {
      setLoading(false);
    }
  };

  const saveConfig = async () => {
    try {
      const bridge = getBridge();
      await bridge.saveConfig('synthesis', values as Record<string, unknown>);
      setConfig('synthesis', values as unknown as SynthesisConfig);
      setIsDirty(false);
    } catch (err) {
      setError(String(err));
      showToast('error', `Failed to save: ${err}`);
    }
  };

  const handleExit = () => {
    if (isDirty) {
      showConfirm({
        title: 'Unsaved Changes',
        message: 'Discard changes to synthesis configuration?',
        confirmLabel: 'Discard',
        cancelLabel: 'Cancel',
        onConfirm: () => {
          goBack();
        },
      });
    } else {
      goBack();
    }
  };

  useInput((input, key) => {
    if (loading) return;

    if (editMode) {
      if (key.escape) {
        setEditMode(false);
      } else if (key.return) {
        const field = fields[selectedIndex];
        let newValue: number | string = editValue;

        if (editValue !== 'inf' && editValue !== '-inf') {
          newValue = parseFloat(editValue) || 0;
        }

        setValues({ ...values, [field.key]: newValue });
        setIsDirty(true);
        setEditMode(false);
      }
      return;
    }

    if (key.upArrow) {
      setSelectedIndex(Math.max(0, selectedIndex - 1));
    } else if (key.downArrow) {
      setSelectedIndex(Math.min(fields.length - 1, selectedIndex + 1));
    } else if (key.return || input === 'e') {
      const field = fields[selectedIndex];
      if (field.type === 'boolean') {
        setValues({ ...values, [field.key]: !values[field.key] });
        setIsDirty(true);
      } else {
        setEditValue(String(values[field.key] || ''));
        setEditMode(true);
      }
    } else if (input === 's') {
      saveConfig();
    } else if (input === 'r') {
      // Navigate to retrosynthesis config
      if (isDirty) {
        showConfirm({
          title: 'Unsaved Changes',
          message: 'Save changes before configuring retrosynthesis?',
          confirmLabel: 'Save & Continue',
          cancelLabel: 'Discard',
          onConfirm: async () => {
            await saveConfig();
            setScreen('configRetrosynthesis');
          },
          onCancel: () => {
            setScreen('configRetrosynthesis');
          },
        });
      } else {
        setScreen('configRetrosynthesis');
      }
    } else if (key.escape || key.leftArrow) {
      handleExit();
    }
  });

  if (loading) {
    return (
      <Box flexDirection="column" flexGrow={1} padding={1}>
        <Header title="Synthesis Config" />
        <Spinner label="Loading configuration..." />
      </Box>
    );
  }

  return (
    <Box flexDirection="column" flexGrow={1} padding={1}>
      <Header title="Synthesis Config" subtitle={isDirty ? 'config_synthesis.yml *' : 'config_synthesis.yml'} />

      {error && (
        <Box marginY={1}>
          <Text color={theme.palette.error}>Error: {error}</Text>
        </Box>
      )}

      <Box flexDirection="column" marginY={1}>
        {fields.map((field, index) => {
          const isSelected = index === selectedIndex;
          const isEditing = isSelected && editMode;
          const value = values[field.key];
          const showSection = field.section && (index === 0 || fields[index - 1]?.section !== field.section);
          const rowWidth = Math.max(1, terminalWidth - 4);
          const labelWidth = Math.min(20, Math.max(8, rowWidth - 6));
          const valueWidth = Math.max(1, rowWidth - labelWidth - 4);

          return (
            <React.Fragment key={field.key}>
              {showSection && (
                <Box marginTop={index > 0 ? 1 : 0}>
                  <Text color={theme.palette.accent} bold>─ {field.section} ─</Text>
                </Box>
              )}
              <Box flexDirection="column">
                <Box width={rowWidth}>
                  <Text color={isSelected ? theme.palette.primary : theme.palette.text}>
                    {isSelected ? '▶ ' : '  '}
                  </Text>
                  <Box width={labelWidth}>
                    <Text color={theme.palette.textMuted}>{field.label}</Text>
                  </Box>
                  {isEditing ? (
                    <Box width={valueWidth}>
                      <TextInput
                        value={editValue}
                        onChange={setEditValue}
                        focus={true}
                      />
                    </Box>
                    ) : (
                    <Box width={valueWidth}>
                      <Text
                        color={field.type === 'boolean'
                          ? (value ? theme.palette.success : theme.palette.error)
                          : theme.palette.accent}
                      >
                        {truncateMiddle(
                          field.type === 'boolean' ? (value ? 'Yes' : 'No') : String(value ?? ''),
                          valueWidth,
                        ).padEnd(valueWidth, ' ')}
                      </Text>
                    </Box>
                  )}
                </Box>
              </Box>
            </React.Fragment>
          );
        })}
      </Box>

      <Box height={1} paddingLeft={2}>
        <Text color={theme.palette.textMuted} wrap="truncate-end">
          {selectedDescription ? `Info: ${selectedDescription}` : ' '}
        </Text>
      </Box>
      
      <Footer />
    </Box>
  );
}

export default ConfigSynthesis;
