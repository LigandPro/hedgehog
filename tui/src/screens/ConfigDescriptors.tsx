import React, { useState, useEffect } from 'react';
import { Box, Text, useInput } from 'ink';
import TextInput from 'ink-text-input';
import { Header } from '../components/Header.js';
import { Footer } from '../components/Footer.js';
import { AppShell } from '../components/AppShell.js';
import { Spinner } from '../components/Spinner.js';
import { DescriptorTable, type DescriptorRow } from '../components/DescriptorTable.js';
import { useStore } from '../store/index.js';
import { getBridge } from '../services/python-bridge.js';
import { useInputLock } from '../hooks/useInputLock.js';
import { getByPath, setByPath } from '../utils/object-path.js';
import { useTheme } from '../theme/context.js';
import { useTerminalSize } from '../hooks/useTerminalSize.js';
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
  'has_spider_side_chains': 'Spider Side Chains',
  'fraction_ring_system': 'Fraction Ring System',
  'tpsa': 'TPSA',
  'qed': 'QED',
};

interface SettingField {
  key: string;
  label: string;
  type: 'boolean' | 'number' | 'text' | 'map';
  path?: string[];
  description?: string;
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

const settingsFields: SettingField[] = [
  { key: 'run', label: 'Run Stage', type: 'boolean', description: 'Enable/disable descriptors calculation stage' },
  { key: 'batch_size', label: 'Batch Size', type: 'number', description: 'Molecules per batch for memory efficiency' },
  { key: 'filter_data', label: 'Filter Data', type: 'boolean', description: 'Apply descriptor borders to filter molecules' },
  {
    key: 'enabled',
    label: 'Enable Structural Constraints',
    type: 'boolean',
    path: ['structural_constraints', 'enabled'],
    description: 'Apply configured type, element, ring, and chain limits',
  },
  {
    key: 'max_n_or_o_atoms',
    label: 'Max N or O Atoms',
    type: 'number',
    path: ['structural_constraints', 'max_n_or_o_atoms'],
    description: 'Upper bound for total N + O atoms in a molecule',
  },
  {
    key: 'max_small_rings_3_4',
    label: 'Max 3/4-Atom Rings',
    type: 'number',
    path: ['structural_constraints', 'max_small_rings_3_4'],
    description: 'Upper bound for the number of 3/4-membered rings',
  },
  {
    key: 'max_acyclic_chain_length',
    label: 'Max Acyclic Chain',
    type: 'number',
    path: ['structural_constraints', 'max_acyclic_chain_length'],
    description: 'Upper bound for longest acyclic chain length',
  },
  {
    key: 'type_limits',
    label: 'Type Limits',
    type: 'map',
    path: ['structural_constraints', 'type_limits'],
    description: 'Type-specific maxima, e.g. .=O=2, Car=6, Nd+=1',
  },
  {
    key: 'element_limits',
    label: 'Element Limits',
    type: 'map',
    path: ['structural_constraints', 'element_limits'],
    description: 'Element-specific maxima, e.g. N=4, O=6, S=2',
  },
];

function formatSettingValue(fieldType: string, value: unknown): string {
  if (fieldType === 'boolean') {
    return value ? 'Yes' : 'No';
  }
  if (fieldType === 'map') {
    if (!value || typeof value !== 'object' || Array.isArray(value)) {
      return '';
    }
    return Object.entries(value as Record<string, unknown>)
      .map(([k, v]) => `${k}=${String(v)}`)
      .join(', ');
  }
  if (Array.isArray(value)) {
    return value.join(', ');
  }
  if (value === undefined || value === null) {
    return '';
  }
  return String(value);
}

function parseMapInput(input: string): Record<string, number> {
  const result: Record<string, number> = {};
  const parts = input
    .split(',')
    .map((part) => part.trim())
    .filter(Boolean);

  for (const part of parts) {
    const sepIndex = part.includes('=') ? part.indexOf('=') : part.indexOf(':');
    if (sepIndex <= 0) {
      continue;
    }
    const key = part.slice(0, sepIndex).trim();
    const rawValue = part.slice(sepIndex + 1).trim();
    const value = Number(rawValue);
    if (!key || !Number.isFinite(value) || value < 0) {
      continue;
    }
    result[key] = Math.floor(value);
  }

  return result;
}

function parseSettingsEditValue(
  field: SettingField,
  value: string,
): string | number | string[] | Record<string, number> {
  if (field.type === 'number') {
    const parsed = Number.parseInt(value, 10);
    return Number.isNaN(parsed) ? 0 : parsed;
  }
  if (field.type === 'map') {
    return parseMapInput(value);
  }
  return value;
}

type ViewMode = 'settings' | 'borders';

export function ConfigDescriptors(): React.ReactElement {
  const { theme } = useTheme();
  const { width: terminalWidth } = useTerminalSize();
  const setScreen = useStore((state) => state.setScreen);
  const goBack = useStore((state) => state.goBack);
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
  useInputLock(settingsEditMode);
  const selectedSettingsField = settingsFields[settingsIndex];
  const selectedDescription = selectedSettingsField?.description || '';

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
          goBack();
        },
      });
    } else {
      goBack();
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
          const newValue = parseSettingsEditValue(field, settingsEditValue);
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
          const formatted = formatSettingValue(field.type, val);
          setSettingsEditValue(formatted);
          setSettingsEditMode(true);
        }
      } else if (input === 'b') {
        setViewMode('borders');
      } else if (input === 's') {
        saveConfig();
      } else if (key.escape || key.leftArrow) {
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
      <AppShell
        padding={1}
        header={<Header title="Descriptors Config" />}
        footer={<Footer />}
      >
        <Spinner label="Loading configuration..." />
      </AppShell>
    );
  }

  return (
    <AppShell
      padding={1}
      header={<Header title="Descriptors Config" subtitle={isDirty ? 'config_descriptors.yml *' : 'config_descriptors.yml'} />}
      footer={<Footer shortcuts={viewMode === 'settings' ? settingsShortcuts : bordersShortcuts} />}
    >
      {error && (
        <Box marginY={1}>
          <Text color={theme.palette.error}>Error: {error}</Text>
        </Box>
      )}

      {viewMode === 'settings' ? (
        <>
          <Box marginY={1}>
            <Text color={theme.palette.accent} bold>─ Settings ─</Text>
          </Box>
          <Box flexDirection="column">
            {settingsFields.map((field, index) => {
              const isSelected = index === settingsIndex;
              const isEditing = isSelected && settingsEditMode;
              const value = getSettingValue(field);
              const rowWidth = Math.max(1, terminalWidth - 4);
              const labelWidth = Math.min(18, Math.max(8, rowWidth - 6));
              const valueWidth = Math.max(1, rowWidth - labelWidth - 4);

              return (
                <Box key={field.key} flexDirection="column">
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
                          value={settingsEditValue}
                          onChange={setSettingsEditValue}
                          focus={true}
                        />
                      </Box>
                    ) : (
                      <Box width={valueWidth}>
                        <Text
                          color={
                            field.type === 'boolean'
                              ? (value ? theme.palette.success : theme.palette.error)
                              : theme.palette.accent
                          }
                        >
                          {truncateMiddle(formatSettingValue(field.type, value), valueWidth).padEnd(valueWidth, ' ')}
                        </Text>
                      </Box>
                    )}
                  </Box>
                </Box>
              );
            })}
          </Box>
          <Box height={1} paddingLeft={2}>
            <Text color={theme.palette.textMuted} wrap="truncate-end">
              {selectedDescription ? `Info: ${selectedDescription}` : ' '}
            </Text>
          </Box>
          <Box marginTop={1}>
            <Text color={theme.palette.textMuted}>Press 'b' to edit descriptor borders</Text>
          </Box>
        </>
      ) : (
        <>
          <DescriptorTable
            descriptors={descriptors}
            onChange={handleDescriptorChange}
            onBack={handleBack}
          />
        </>
      )}
    </AppShell>
  );
}

export default ConfigDescriptors;
