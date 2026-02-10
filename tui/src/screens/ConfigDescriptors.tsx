import React, { useState, useEffect } from 'react';
import { Box, Text, useInput } from 'ink';
import TextInput from 'ink-text-input';
import { Header } from '../components/Header.js';
import { Footer } from '../components/Footer.js';
import { Spinner } from '../components/Spinner.js';
import { DescriptorTable, type DescriptorRow } from '../components/DescriptorTable.js';
import { useStore } from '../store/index.js';
import { getBridge } from '../services/python-bridge.js';
import { getByPath, setByPath } from '../utils/object-path.js';
import type { DescriptorsConfig } from '../types/index.js';


// Descriptor presets
const PRESETS: Record<string, Record<string, [number, number]>> = {
  'Drug-like': {
    molWt: [200, 500],
    logP: [-0.4, 5.6],
    hbd: [0, 5],
    hba: [0, 10],
    tpsa: [0, 140],
    n_rot_bonds: [0, 10],
  },
  'Lead-like': {
    molWt: [200, 350],
    logP: [-1, 3],
    hbd: [0, 3],
    hba: [0, 6],
    n_rot_bonds: [0, 7],
  },
  'Fragment-like': {
    molWt: [0, 300],
    logP: [-3, 3],
    hbd: [0, 3],
    hba: [0, 3],
    n_heavy_atoms: [0, 17],
  },
};

// List of descriptors with display names
const DESCRIPTOR_NAMES: Record<string, string> = {
  'n_atoms': 'Number of Atoms',
  'n_heavy_atoms': 'Heavy Atoms',
  'n_het_atoms': 'Heteroatoms',
  'n_N_atoms': 'Nitrogen Atoms',
  'fN_atoms': 'Fraction N Atoms',
  'molWt': 'Molecular Weight',
  'logP': 'logP',
  'clogP': 'clogP',
  'sw': 'Water Solubility (Sw)',
  'ring_size': 'Ring Size',
  'n_rings': 'Number of Rings',
  'n_aroma_rings': 'Aromatic Rings',
  'n_fused_aromatic_rings': 'Fused Aromatic Rings',
  'n_rigid_bonds': 'Rigid Bonds',
  'n_rot_bonds': 'Rotatable Bonds',
  'hbd': 'H-Bond Donors',
  'hba': 'H-Bond Acceptors',
  'fsp3': 'Fraction SP3',
  'tpsa': 'TPSA',
  'qed': 'QED',
};

interface SettingField {
  key: string;
  label: string;
  type: 'boolean' | 'number' | 'text';
  path?: string[];
  description?: string;
}

const settingsFields: SettingField[] = [
  { key: 'run', label: 'Run Stage', type: 'boolean', description: 'Enable/disable descriptors calculation stage' },
  { key: 'batch_size', label: 'Batch Size', type: 'number', description: 'Molecules per batch for memory efficiency' },
  { key: 'filter_data', label: 'Filter Data', type: 'boolean', description: 'Apply descriptor borders to filter molecules' },
  { key: 'allowed_chars', label: 'Allowed Atoms', type: 'text', path: ['borders', 'allowed_chars'], description: 'Allowed atom types (e.g., C, N, O, S, F, Cl, Br)' },
  { key: 'charged_mol_allowed', label: 'Charged Allowed', type: 'boolean', path: ['borders', 'charged_mol_allowed'], description: 'Allow charged molecules to pass' },
  { key: 'filter_charged_mol', label: 'Filter Charged', type: 'boolean', path: ['borders', 'filter_charged_mol'], description: 'Remove molecules with formal charges' },
];

type ViewMode = 'settings' | 'borders';

export function ConfigDescriptors(): React.ReactElement {
  const setScreen = useStore((state) => state.setScreen);
  const config = useStore((state) => state.configs.descriptors);
  const setConfig = useStore((state) => state.setConfig);
  const isBackendReady = useStore((state) => state.isBackendReady);
  const showToast = useStore((state) => state.showToast);
  const showConfirm = useStore((state) => state.showConfirm);

  const [loading, setLoading] = useState(!config);
  const [viewMode, setViewMode] = useState<ViewMode>('settings');
  const [descriptors, setDescriptors] = useState<DescriptorRow[]>([]);
  const [rawConfig, setRawConfig] = useState<Record<string, unknown> | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [isDirty, setIsDirty] = useState(false);

  // Settings view state
  const [settingsIndex, setSettingsIndex] = useState(0);
  const [settingsEditMode, setSettingsEditMode] = useState(false);
  const [settingsEditValue, setSettingsEditValue] = useState('');

  useEffect(() => {
    if (!config && isBackendReady) {
      loadConfig();
    } else if (config) {
      parseConfig(config as unknown as Record<string, unknown>);
      setRawConfig(config as unknown as Record<string, unknown>);
      setLoading(false);
    }
  }, [config, isBackendReady]);

  const parseConfig = (cfg: Record<string, unknown>) => {
    const rows: DescriptorRow[] = [];
    const borders = (cfg.borders || {}) as Record<string, unknown>;
    
    for (const [key, displayName] of Object.entries(DESCRIPTOR_NAMES)) {
      const minKey = `${key}_min`;
      const maxKey = `${key}_max`;
      
      rows.push({
        name: key,
        displayName,
        min: (borders[minKey] as number | string) ?? 0,
        max: (borders[maxKey] as number | string) ?? 100,
      });
    }
    
    setDescriptors(rows);
  };

  const getSettingValue = (field: SettingField): unknown => {
    if (!rawConfig) return undefined;
    if (field.path) {
      return getByPath(rawConfig, field.path);
    }
    return (rawConfig as Record<string, unknown>)[field.key];
  };

  const setSettingValue = (field: SettingField, newValue: unknown) => {
    if (!rawConfig) return;
    let newConfig: Record<string, unknown>;
    if (field.path) {
      newConfig = setByPath(rawConfig, field.path, newValue);
    } else {
      newConfig = { ...rawConfig, [field.key]: newValue };
    }
    setRawConfig(newConfig);
    setIsDirty(true);
  };

  const loadConfig = async () => {
    try {
      const bridge = getBridge();
      const data = await bridge.loadConfig('descriptors');
      setConfig('descriptors', data as unknown as DescriptorsConfig);
      setRawConfig(data);
      parseConfig(data);
    } catch (err) {
      setError(String(err));
    } finally {
      setLoading(false);
    }
  };

  const applyPreset = (presetName: string) => {
    const preset = PRESETS[presetName];
    if (!preset) return;

    const newDescriptors = descriptors.map(desc => {
      const presetValues = preset[desc.name];
      if (presetValues) {
        return { ...desc, min: presetValues[0], max: presetValues[1] };
      }
      return desc;
    });

    setDescriptors(newDescriptors);
    setIsDirty(true);
    showToast('info', `Applied preset: ${presetName}`);
  };

  const saveConfig = async () => {
    if (!rawConfig) return;

    try {
      const bridge = getBridge();
      const newBorders: Record<string, unknown> = { ...(rawConfig.borders as Record<string, unknown>) };

      for (const desc of descriptors) {
        newBorders[`${desc.name}_min`] = desc.min;
        newBorders[`${desc.name}_max`] = desc.max;
      }

      const newConfig = { ...rawConfig, borders: newBorders };
      await bridge.saveConfig('descriptors', newConfig);
      setConfig('descriptors', newConfig as unknown as DescriptorsConfig);
      setRawConfig(newConfig);
      setIsDirty(false);
    } catch (err) {
      setError(String(err));
      showToast('error', `Failed to save: ${err}`);
    }
  };

  const handleDescriptorChange = (newDescriptors: DescriptorRow[]) => {
    setDescriptors(newDescriptors);
    setIsDirty(true);
  };

  const handleBack = () => {
    if (viewMode === 'borders') {
      setViewMode('settings');
    } else if (isDirty) {
      showConfirm({
        title: 'Unsaved Changes',
        message: 'Discard changes to descriptors configuration?',
        confirmLabel: 'Discard',
        cancelLabel: 'Cancel',
        onConfirm: () => {
          setScreen('welcome');
        },
      });
    } else {
      setScreen('welcome');
    }
  };

  useInput((input, key) => {
    if (loading) return;
    
    // Settings view mode
    if (viewMode === 'settings') {
      if (settingsEditMode) {
        if (key.escape) {
          setSettingsEditMode(false);
        } else if (key.return) {
          const field = settingsFields[settingsIndex];
          let newValue: string | number | string[] = settingsEditValue;
          if (field.type === 'number') {
            newValue = parseInt(settingsEditValue, 10) || 0;
          } else if (field.key === 'allowed_chars') {
            newValue = settingsEditValue.split(',').map(s => s.trim());
          }
          setSettingValue(field, newValue);
          setSettingsEditMode(false);
        }
        return;
      }

      if (key.upArrow) {
        setSettingsIndex(Math.max(0, settingsIndex - 1));
      } else if (key.downArrow) {
        setSettingsIndex(Math.min(settingsFields.length - 1, settingsIndex + 1));
      } else if (key.return || input === 'e') {
        const field = settingsFields[settingsIndex];
        if (field.type === 'boolean') {
          setSettingValue(field, !getSettingValue(field));
        } else {
          const val = getSettingValue(field);
          setSettingsEditValue(Array.isArray(val) ? val.join(', ') : String(val || ''));
          setSettingsEditMode(true);
        }
      } else if (input === 'b') {
        setViewMode('borders');
      } else if (input === 's') {
        saveConfig();
      } else if (key.escape || key.leftArrow || input === 'q') {
        handleBack();
      }
      return;
    }
    
    // Borders view mode - handled by DescriptorTable
    if (input === 's') {
      saveConfig();
    } else if (input === '1') {
      applyPreset('Drug-like');
    } else if (input === '2') {
      applyPreset('Lead-like');
    } else if (input === '3') {
      applyPreset('Fragment-like');
    }
  });

  const settingsShortcuts = [
    { key: '↑↓', label: 'Navigate' },
    { key: 'e/Enter', label: 'Edit' },
    { key: 'b', label: 'Borders' },
    { key: 's', label: 'Save' },
    { key: '←/Esc', label: 'Back' },
  ];

  const bordersShortcuts = [
    { key: '↑↓', label: 'Navigate' },
    { key: 'm/x', label: 'Edit min/max' },
    { key: '1/2/3', label: 'Presets' },
    { key: 's', label: 'Save' },
    { key: '←/Esc', label: 'Settings' },
  ];

  if (loading) {
    return (
      <Box flexDirection="column" padding={1}>
        <Header title="Descriptors Config" />
        <Spinner label="Loading configuration..." />
      </Box>
    );
  }

  return (
    <Box flexDirection="column" padding={1}>
      <Header title="Descriptors Config" subtitle={isDirty ? 'config_descriptors.yml *' : 'config_descriptors.yml'} />

      {error && (
        <Box marginY={1}>
          <Text color="red">Error: {error}</Text>
        </Box>
      )}

      {viewMode === 'settings' ? (
        <>
          <Box marginY={1}>
            <Text color="cyan" bold>─ Settings ─</Text>
          </Box>
          <Box flexDirection="column">
            {settingsFields.map((field, index) => {
              const isSelected = index === settingsIndex;
              const isEditing = isSelected && settingsEditMode;
              const value = getSettingValue(field);

              return (
                <Box key={field.key} flexDirection="column">
                  <Box>
                    <Text color={isSelected ? 'cyan' : 'white'}>
                      {isSelected ? '▶ ' : '  '}
                    </Text>
                    <Text dimColor>{field.label.padEnd(18)}</Text>
                    {isEditing ? (
                      <Box>
                        <TextInput
                          value={settingsEditValue}
                          onChange={setSettingsEditValue}
                          focus={true}
                        />
                      </Box>
                    ) : (
                      <Text color={
                        field.type === 'boolean'
                          ? (value ? 'green' : 'red')
                          : 'yellow'
                      }>
                        {field.type === 'boolean'
                          ? (value ? 'Yes' : 'No')
                          : Array.isArray(value)
                            ? value.join(', ')
                            : String(value ?? '')}
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
          <Box marginTop={1}>
            <Text dimColor>Press 'b' to edit descriptor borders</Text>
          </Box>
          <Footer shortcuts={settingsShortcuts} />
        </>
      ) : (
        <>
          <DescriptorTable
            descriptors={descriptors}
            onChange={handleDescriptorChange}
            onBack={handleBack}
          />
          <Footer shortcuts={bordersShortcuts} />
        </>
      )}
    </Box>
  );
}

export default ConfigDescriptors;
