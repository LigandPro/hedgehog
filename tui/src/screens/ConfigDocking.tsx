import React, { useState, useEffect } from 'react';
import { Box, Text, useInput } from 'ink';
import TextInput from 'ink-text-input';
import { Header } from '../components/Header.js';
import { Footer } from '../components/Footer.js';
import { Spinner } from '../components/Spinner.js';
import { FileBrowser } from '../components/FileBrowser.js';
import { useStore } from '../store/index.js';
import { getBridge } from '../services/python-bridge.js';
import type { DockingConfig } from '../types/index.js';

interface FormField {
  key: string;
  label: string;
  type: 'text' | 'number' | 'boolean' | 'select' | 'path';
  extensions?: string[];
  options?: string[];
  path?: string[];  // nested path in config
}

const fields: FormField[] = [
  { key: 'run', label: 'Run Stage', type: 'boolean' },
  { key: 'tools', label: 'Tools', type: 'select', options: ['both', 'smina', 'gnina'] },
  { key: 'receptor_pdb', label: 'Receptor PDB', type: 'path', extensions: ['.pdb'] },
  { key: 'auto_run', label: 'Auto Run', type: 'boolean' },
  { key: 'run_in_background', label: 'Run in Background', type: 'boolean' },
  { key: 'smina_autobox_ligand', label: 'Smina Autobox Ligand', type: 'text', path: ['smina_config', 'autobox_ligand'] },
  { key: 'smina_autobox_add', label: 'Smina Autobox Add', type: 'number', path: ['smina_config', 'autobox_add'] },
  { key: 'smina_exhaustiveness', label: 'Smina Exhaustiveness', type: 'number', path: ['smina_config', 'exhaustiveness'] },
  { key: 'smina_num_modes', label: 'Smina Num Modes', type: 'number', path: ['smina_config', 'num_modes'] },
  { key: 'smina_cpu', label: 'Smina CPUs', type: 'number', path: ['smina_config', 'cpu'] },
  { key: 'gnina_autobox_ligand', label: 'GNINA Autobox Ligand', type: 'text', path: ['gnina_config', 'autobox_ligand'] },
  { key: 'gnina_autobox_add', label: 'GNINA Autobox Add', type: 'number', path: ['gnina_config', 'autobox_add'] },
];

export function ConfigDocking(): React.ReactElement {
  const setScreen = useStore((state) => state.setScreen);
  const config = useStore((state) => state.configs.docking);
  const setConfig = useStore((state) => state.setConfig);
  const isBackendReady = useStore((state) => state.isBackendReady);
  
  const [loading, setLoading] = useState(!config);
  const [selectedIndex, setSelectedIndex] = useState(0);
  const [editMode, setEditMode] = useState(false);
  const [editValue, setEditValue] = useState('');
  const [values, setValues] = useState<Record<string, any>>({});
  const [error, setError] = useState<string | null>(null);
  const [saved, setSaved] = useState(false);
  const [browsingField, setBrowsingField] = useState<FormField | null>(null);

  const getValue = (field: FormField): any => {
    if (field.path) {
      let val: any = values;
      for (const key of field.path) {
        val = val?.[key];
      }
      return val;
    }
    return values[field.key];
  };

  const setValue = (field: FormField, newValue: any) => {
    const newValues = { ...values };
    if (field.path) {
      let obj = newValues;
      for (let i = 0; i < field.path.length - 1; i++) {
        const key = field.path[i];
        obj[key] = { ...obj[key] };
        obj = obj[key];
      }
      obj[field.path[field.path.length - 1]] = newValue;
    } else {
      newValues[field.key] = newValue;
    }
    setValues(newValues);
  };

  useEffect(() => {
    if (!config && isBackendReady) {
      loadConfig();
    } else if (config) {
      setValues(config as any);
      setLoading(false);
    }
  }, [config, isBackendReady]);

  const loadConfig = async () => {
    try {
      const bridge = getBridge();
      const data = await bridge.loadConfig('docking') as unknown as DockingConfig;
      setConfig('docking', data);
      setValues(data as any);
    } catch (err) {
      setError(String(err));
    } finally {
      setLoading(false);
    }
  };

  const saveConfig = async () => {
    try {
      const bridge = getBridge();
      await bridge.saveConfig('docking', values);
      setConfig('docking', values as unknown as DockingConfig);
      setSaved(true);
      setTimeout(() => setSaved(false), 2000);
    } catch (err) {
      setError(String(err));
    }
  };

  useInput((input, key) => {
    if (loading || browsingField) return;
    
    if (editMode) {
      if (key.escape) {
        setEditMode(false);
      } else if (key.return) {
        const field = fields[selectedIndex];
        let newValue: any = editValue;
        
        if (field.type === 'number') {
          newValue = parseFloat(editValue) || 0;
        }
        
        setValue(field, newValue);
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
      const value = getValue(field);
      
      if (field.type === 'boolean') {
        setValue(field, !value);
      } else if (field.type === 'select' && field.options) {
        const currentIdx = field.options.indexOf(value);
        const nextIdx = (currentIdx + 1) % field.options.length;
        setValue(field, field.options[nextIdx]);
      } else {
        setEditValue(String(value || ''));
        setEditMode(true);
      }
    } else if (input === 's') {
      saveConfig();
    } else if (key.escape || input === 'q') {
      setScreen('welcome');
    }
  });

  const handlePathSelect = (path: string) => {
    if (browsingField) {
      setValue(browsingField, path);
    }
    setBrowsingField(null);
  };

  const handleBrowseCancel = () => {
    setBrowsingField(null);
  };

  const shortcuts = [
    { key: '↑↓', label: 'Navigate' },
    { key: 'e/Enter', label: 'Edit' },
    { key: 's', label: 'Save' },
    { key: 'Esc', label: 'Back' },
  ];

  if (loading) {
    return (
      <Box flexDirection="column" padding={1}>
        <Header title="Docking Config" />
        <Spinner label="Loading configuration..." />
      </Box>
    );
  }

  if (browsingField) {
    const currentValue = String(getValue(browsingField) || process.cwd());
    return (
      <Box flexDirection="column" padding={1}>
        <Header title="Select File" subtitle={browsingField.label} />
        <FileBrowser
          initialPath={currentValue}
          extensions={browsingField.extensions}
          onSelect={handlePathSelect}
          onCancel={handleBrowseCancel}
        />
      </Box>
    );
  }

  return (
    <Box flexDirection="column" padding={1}>
      <Header title="Docking Config" subtitle="config_docking.yml" />
      
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
          const value = getValue(field);

          return (
            <Box key={field.key}>
              <Text color={isSelected ? 'cyan' : 'white'}>
                {isSelected ? '▶ ' : '  '}
              </Text>
              <Text dimColor>{field.label.padEnd(22)}</Text>
              {isEditing ? (
                <Box>
                  <TextInput
                    value={editValue}
                    onChange={setEditValue}
                    focus={true}
                  />
                </Box>
              ) : (
                <Text color={
                  field.type === 'boolean' 
                    ? (value ? 'green' : 'red') 
                    : field.type === 'select' 
                      ? 'magenta' 
                      : 'yellow'
                }>
                  {field.type === 'boolean' 
                    ? (value ? 'Yes' : 'No') 
                    : String(value ?? '')}
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

export default ConfigDocking;
