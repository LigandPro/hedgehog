import React, { useState, useEffect } from 'react';
import { Box, Text, useInput } from 'ink';
import TextInput from 'ink-text-input';
import { Header } from '../components/Header.js';
import { Footer } from '../components/Footer.js';
import { Spinner } from '../components/Spinner.js';
import { useStore } from '../store/index.js';
import { getBridge } from '../services/python-bridge.js';
import type { SynthesisConfig } from '../types/index.js';

interface FormField {
  key: keyof SynthesisConfig;
  label: string;
  type: 'number' | 'boolean';
  description?: string;
}

const fields: FormField[] = [
  { key: 'run', label: 'Run Stage', type: 'boolean' },
  { key: 'filter_solved_only', label: 'Filter Solved Only', type: 'boolean', description: 'Only keep molecules with retrosynthesis path' },
  { key: 'sa_score_min', label: 'SA Score Min', type: 'number', description: 'Synthetic Accessibility (1-10, lower = easier)' },
  { key: 'sa_score_max', label: 'SA Score Max', type: 'number' },
  { key: 'syba_score_min', label: 'SYBA Score Min', type: 'number', description: 'Synthetic Bayesian Accessibility' },
  { key: 'syba_score_max', label: 'SYBA Score Max', type: 'number' },
  { key: 'ra_score_min', label: 'RA Score Min', type: 'number', description: 'Retrosynthesis probability (0-1)' },
  { key: 'ra_score_max', label: 'RA Score Max', type: 'number' },
];

export function ConfigSynthesis(): React.ReactElement {
  const setScreen = useStore((state) => state.setScreen);
  const config = useStore((state) => state.configs.synthesis);
  const setConfig = useStore((state) => state.setConfig);
  const isBackendReady = useStore((state) => state.isBackendReady);
  
  const [loading, setLoading] = useState(!config);
  const [selectedIndex, setSelectedIndex] = useState(0);
  const [editMode, setEditMode] = useState(false);
  const [editValue, setEditValue] = useState('');
  const [values, setValues] = useState<Partial<SynthesisConfig>>(config || {});
  const [error, setError] = useState<string | null>(null);
  const [saved, setSaved] = useState(false);

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
      setSaved(true);
      setTimeout(() => setSaved(false), 2000);
    } catch (err) {
      setError(String(err));
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
      } else {
        setEditValue(String(values[field.key] || ''));
        setEditMode(true);
      }
    } else if (input === 's') {
      saveConfig();
    } else if (key.escape || input === 'q') {
      setScreen('welcome');
    }
  });

  const shortcuts = [
    { key: '↑↓', label: 'Navigate' },
    { key: 'e/Enter', label: 'Edit' },
    { key: 's', label: 'Save' },
    { key: 'Esc', label: 'Back' },
  ];

  if (loading) {
    return (
      <Box flexDirection="column" padding={1}>
        <Header title="Synthesis Config" />
        <Spinner label="Loading configuration..." />
      </Box>
    );
  }

  return (
    <Box flexDirection="column" padding={1}>
      <Header title="Synthesis Config" subtitle="config_synthesis.yml" />
      
      {error && (
        <Box marginY={1}>
          <Text color="red">Error: {error}</Text>
        </Box>
      )}
      
      {saved && (
        <Box marginY={1}>
          <Text color="green">✓ Configuration saved</Text>
        </Box>
      )}
      
      <Box flexDirection="column" marginY={1}>
        {fields.map((field, index) => {
          const isSelected = index === selectedIndex;
          const isEditing = isSelected && editMode;
          const value = values[field.key];

          return (
            <Box key={field.key} flexDirection="column">
              <Box>
                <Text color={isSelected ? 'cyan' : 'white'}>
                  {isSelected ? '▶ ' : '  '}
                </Text>
                <Text dimColor>{field.label.padEnd(20)}</Text>
                {isEditing ? (
                  <Box>
                    <TextInput
                      value={editValue}
                      onChange={setEditValue}
                      focus={true}
                    />
                  </Box>
                ) : (
                  <Text color={field.type === 'boolean' ? (value ? 'green' : 'red') : 'yellow'}>
                    {field.type === 'boolean' ? (value ? 'Yes' : 'No') : String(value ?? '')}
                  </Text>
                )}
              </Box>
              {isSelected && field.description && (
                <Box paddingLeft={4}>
                  <Text dimColor italic>{field.description}</Text>
                </Box>
              )}
            </Box>
          );
        })}
      </Box>
      
      <Footer shortcuts={shortcuts} />
    </Box>
  );
}

export default ConfigSynthesis;
