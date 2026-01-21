import React, { useState, useEffect } from 'react';
import { Box, Text, useInput } from 'ink';
import TextInput from 'ink-text-input';
import { Header } from '../components/Header.js';
import { Footer } from '../components/Footer.js';
import { Spinner } from '../components/Spinner.js';
import { useStore } from '../store/index.js';
import { getBridge } from '../services/python-bridge.js';
import type { MainConfig } from '../types/index.js';

interface FormField {
  key: keyof MainConfig;
  label: string;
  type: 'text' | 'number' | 'boolean' | 'path';
}

const fields: FormField[] = [
  { key: 'generated_mols_path', label: 'Input Molecules', type: 'path' },
  { key: 'target_mols_path', label: 'Target Molecules', type: 'path' },
  { key: 'folder_to_save', label: 'Output Folder', type: 'path' },
  { key: 'n_jobs', label: 'Parallel Jobs', type: 'number' },
  { key: 'sample_size', label: 'Sample Size', type: 'number' },
  { key: 'batch_size', label: 'Batch Size', type: 'number' },
  { key: 'save_sampled_mols', label: 'Save Sampled Mols', type: 'boolean' },
];

export function ConfigMain(): React.ReactElement {
  const setScreen = useStore((state) => state.setScreen);
  const config = useStore((state) => state.configs.main);
  const setConfig = useStore((state) => state.setConfig);
  const isBackendReady = useStore((state) => state.isBackendReady);
  
  const [loading, setLoading] = useState(!config);
  const [selectedIndex, setSelectedIndex] = useState(0);
  const [editMode, setEditMode] = useState(false);
  const [editValue, setEditValue] = useState('');
  const [values, setValues] = useState<Partial<MainConfig>>(config || {});
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
      const data = await bridge.loadConfig('main') as unknown as MainConfig;
      setConfig('main', data);
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
      await bridge.saveConfig('main', values as Record<string, unknown>);
      setConfig('main', values as unknown as MainConfig);
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
        let newValue: string | number | boolean = editValue;
        
        if (field.type === 'number') {
          newValue = parseInt(editValue, 10) || 0;
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
        <Header title="Main Config" />
        <Spinner label="Loading configuration..." />
      </Box>
    );
  }

  return (
    <Box flexDirection="column" padding={1}>
      <Header title="Main Config" subtitle="config.yml" />
      
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
            <Box key={field.key}>
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
                  {field.type === 'boolean' ? (value ? 'Yes' : 'No') : String(value || '')}
                </Text>
              )}
            </Box>
          );
        })}
      </Box>
      
      <Footer shortcuts={shortcuts} />
    </Box>
  );
}

export default ConfigMain;
