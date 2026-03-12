import React, { useState, useEffect } from 'react';
import { Box, Text, useInput } from 'ink';
import TextInput from 'ink-text-input';
import { Header } from '../components/Header.js';
import { Footer } from '../components/Footer.js';
import { AppShell } from '../components/AppShell.js';
import { Spinner } from '../components/Spinner.js';
import { FileBrowser } from '../components/FileBrowser.js';
import { SearchIndicator } from '../components/SearchIndicator.js';
import { useStore } from '../store/index.js';
import { getBridge } from '../services/python-bridge.js';
import { useSearch } from '../hooks/useSearch.js';
import { useInputLock } from '../hooks/useInputLock.js';
import { getByPath, setByPath } from '../utils/object-path.js';
import { useTheme } from '../theme/context.js';
import { useTerminalSize } from '../hooks/useTerminalSize.js';
import type { DockingConfig } from '../types/index.js';

interface FormField {
  key: string;
  label: string;
  type: 'text' | 'number' | 'boolean' | 'select' | 'path';
  extensions?: string[];
  options?: string[];
  path?: string[];
  section?: string;
  description?: string;
}

const fields: FormField[] = [
  // General section
  { key: 'run', label: 'Run Stage', type: 'boolean', section: 'General', description: 'Enable/disable docking stage' },
  { key: 'tools', label: 'Tools', type: 'select', options: ['all', 'smina', 'gnina', 'matcha'], description: 'Docking engines: all, or select one engine here and use YAML for comma-separated combinations' },
  { key: 'receptor_pdb', label: 'Receptor PDB', type: 'path', extensions: ['.pdb'], description: 'Target protein structure file' },
  { key: 'auto_run', label: 'Auto Run', type: 'boolean', description: 'Start docking automatically when ready' },
  { key: 'run_in_background', label: 'Run in Background', type: 'boolean', description: 'Run docking as background process' },

  // Smina - Basic
  { key: 'smina_bin', label: 'Binary Path', type: 'path', path: ['smina_config', 'bin'], section: 'Smina - Basic', description: 'Path to smina executable' },
  { key: 'smina_cpu', label: 'CPUs', type: 'number', path: ['smina_config', 'cpu'], description: 'Number of CPU cores to use' },
  { key: 'smina_seed', label: 'Random Seed', type: 'number', path: ['smina_config', 'seed'], description: 'Random seed for reproducibility' },
  { key: 'smina_exhaustiveness', label: 'Exhaustiveness', type: 'number', path: ['smina_config', 'exhaustiveness'], description: 'Search thoroughness (default: 8, higher = slower but better)' },
  { key: 'smina_num_modes', label: 'Num Modes', type: 'number', path: ['smina_config', 'num_modes'], description: 'Max number of binding modes to generate' },

  // Smina - Search Space
  { key: 'smina_autobox_ligand', label: 'Autobox Ligand', type: 'path', path: ['smina_config', 'autobox_ligand'], section: 'Smina - Search Space', description: 'Reference ligand to define search box' },
  { key: 'smina_autobox_add', label: 'Autobox Add', type: 'number', path: ['smina_config', 'autobox_add'], description: 'Extra space around autobox (Angstroms)' },
  { key: 'smina_center_x', label: 'Center X', type: 'number', path: ['smina_config', 'center_x'], description: 'Search box center X coordinate' },
  { key: 'smina_center_y', label: 'Center Y', type: 'number', path: ['smina_config', 'center_y'], description: 'Search box center Y coordinate' },
  { key: 'smina_center_z', label: 'Center Z', type: 'number', path: ['smina_config', 'center_z'], description: 'Search box center Z coordinate' },
  { key: 'smina_size_x', label: 'Size X', type: 'number', path: ['smina_config', 'size_x'], description: 'Search box size X (Angstroms)' },
  { key: 'smina_size_y', label: 'Size Y', type: 'number', path: ['smina_config', 'size_y'], description: 'Search box size Y (Angstroms)' },
  { key: 'smina_size_z', label: 'Size Z', type: 'number', path: ['smina_config', 'size_z'], description: 'Search box size Z (Angstroms)' },

  // Smina - Flexible Docking
  { key: 'smina_flex', label: 'Flex PDBQT', type: 'path', path: ['smina_config', 'flex'], section: 'Smina - Flexible', description: 'Flexible receptor side chains file' },
  { key: 'smina_flexres', label: 'Flex Residues', type: 'text', path: ['smina_config', 'flexres'], description: 'Flexible residues (e.g., "A:SER45,A:TYR67")' },
  { key: 'smina_flexdist_ligand', label: 'Flexdist Ligand', type: 'path', path: ['smina_config', 'flexdist_ligand'], description: 'Reference for auto-selecting flex residues' },
  { key: 'smina_flexdist', label: 'Flexdist (A)', type: 'number', path: ['smina_config', 'flexdist'], description: 'Distance from ligand for flex residues' },

  // Smina - Scoring
  { key: 'smina_scoring', label: 'Scoring Function', type: 'select', options: ['default', 'vinardo', 'ad4_scoring'], path: ['smina_config', 'scoring'], section: 'Smina - Scoring', description: 'Scoring function: default, vinardo, or AutoDock4' },
  { key: 'smina_custom_scoring', label: 'Custom Scoring File', type: 'path', path: ['smina_config', 'custom_scoring'], description: 'Custom scoring function weights file' },
  { key: 'smina_custom_atoms', label: 'Custom Atoms File', type: 'path', path: ['smina_config', 'custom_atoms'], description: 'Custom atom type parameters' },

  // Smina - Modes
  { key: 'smina_score_only', label: 'Score Only', type: 'boolean', path: ['smina_config', 'score_only'], section: 'Smina - Modes', description: 'Only score input pose, no docking' },
  { key: 'smina_local_only', label: 'Local Only', type: 'boolean', path: ['smina_config', 'local_only'], description: 'Local optimization only (no global search)' },
  { key: 'smina_minimize', label: 'Minimize', type: 'boolean', path: ['smina_config', 'minimize'], description: 'Energy minimize output poses' },
  { key: 'smina_minimize_iters', label: 'Minimize Iters', type: 'number', path: ['smina_config', 'minimize_iters'], description: 'Max minimization iterations' },
  { key: 'smina_randomize_only', label: 'Randomize Only', type: 'boolean', path: ['smina_config', 'randomize_only'], description: 'Generate random poses (no optimization)' },

  // Smina - Output
  { key: 'smina_energy_range', label: 'Energy Range', type: 'number', path: ['smina_config', 'energy_range'], section: 'Smina - Output', description: 'Max energy diff from best pose (kcal/mol)' },
  { key: 'smina_min_rmsd_filter', label: 'Min RMSD Filter', type: 'number', path: ['smina_config', 'min_rmsd_filter'], description: 'Min RMSD between output poses' },
  { key: 'smina_out_flex', label: 'Output Flex', type: 'path', path: ['smina_config', 'out_flex'], description: 'Output file for flexible residues' },
  { key: 'smina_log', label: 'Log File', type: 'path', path: ['smina_config', 'log'], description: 'Docking log file path' },
  { key: 'smina_quiet', label: 'Quiet', type: 'boolean', path: ['smina_config', 'quiet'], description: 'Suppress output messages' },
  { key: 'smina_addH', label: 'Add Hydrogens', type: 'boolean', path: ['smina_config', 'addH'], description: 'Add hydrogens to ligands' },

  // GNINA - Basic
  { key: 'gnina_bin', label: 'Binary Path', type: 'path', path: ['gnina_config', 'bin'], section: 'GNINA - Basic', description: 'Path to gnina executable' },
  { key: 'gnina_env_path', label: 'Environment Path', type: 'path', path: ['gnina_config', 'env_path'], description: 'Conda/venv environment path' },
  { key: 'gnina_cpu', label: 'CPUs', type: 'number', path: ['gnina_config', 'cpu'], description: 'Number of CPU cores' },
  { key: 'gnina_seed', label: 'Random Seed', type: 'number', path: ['gnina_config', 'seed'], description: 'Random seed for reproducibility' },
  { key: 'gnina_exhaustiveness', label: 'Exhaustiveness', type: 'number', path: ['gnina_config', 'exhaustiveness'], description: 'Search thoroughness' },
  { key: 'gnina_num_modes', label: 'Num Modes', type: 'number', path: ['gnina_config', 'num_modes'], description: 'Max binding modes to generate' },

  // GNINA - Search Space
  { key: 'gnina_autobox_ligand', label: 'Autobox Ligand', type: 'path', path: ['gnina_config', 'autobox_ligand'], section: 'GNINA - Search Space', description: 'Reference ligand for search box' },
  { key: 'gnina_autobox_add', label: 'Autobox Add', type: 'number', path: ['gnina_config', 'autobox_add'], description: 'Extra space around autobox' },
  { key: 'gnina_center_x', label: 'Center X', type: 'number', path: ['gnina_config', 'center_x'], description: 'Search box center X' },
  { key: 'gnina_center_y', label: 'Center Y', type: 'number', path: ['gnina_config', 'center_y'], description: 'Search box center Y' },
  { key: 'gnina_center_z', label: 'Center Z', type: 'number', path: ['gnina_config', 'center_z'], description: 'Search box center Z' },
  { key: 'gnina_size_x', label: 'Size X', type: 'number', path: ['gnina_config', 'size_x'], description: 'Search box size X' },
  { key: 'gnina_size_y', label: 'Size Y', type: 'number', path: ['gnina_config', 'size_y'], description: 'Search box size Y' },
  { key: 'gnina_size_z', label: 'Size Z', type: 'number', path: ['gnina_config', 'size_z'], description: 'Search box size Z' },

  // GNINA - Options
  { key: 'gnina_scoring', label: 'Scoring', type: 'select', options: ['default', 'vinardo'], path: ['gnina_config', 'scoring'], section: 'GNINA - Options', description: 'Scoring: default (CNN) or vinardo' },
  { key: 'gnina_minimize', label: 'Minimize', type: 'boolean', path: ['gnina_config', 'minimize'], description: 'Minimize output poses' },
  { key: 'gnina_log', label: 'Log File', type: 'path', path: ['gnina_config', 'log'], description: 'Log file path' },
  { key: 'gnina_output_dir', label: 'Output Directory', type: 'path', path: ['gnina_config', 'output_dir'], description: 'Output directory for results' },
];

type FieldType = 'text' | 'number' | 'boolean' | 'select' | 'path';

function getFieldColor(
  theme: ReturnType<typeof useTheme>['theme'],
  type: FieldType,
  value: unknown,
): string {
  if (type === 'boolean') {
    return value ? theme.palette.success : theme.palette.error;
  }
  if (type === 'select') {
    return theme.palette.primary;
  }
  if (type === 'path') {
    return theme.palette.info;
  }
  return theme.palette.accent;
}

function formatFieldValue(fieldType: FieldType, value: unknown): string {
  if (fieldType === 'boolean') {
    return value ? 'Yes' : 'No';
  }
  if (value !== undefined && value !== null) {
    return String(value);
  }
  return '(not set)';
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

export function ConfigDocking(): React.ReactElement {
  const { theme } = useTheme();
  const { width: terminalWidth, height: terminalHeight } = useTerminalSize();
  const setScreen = useStore((state) => state.setScreen);
  const goBack = useStore((state) => state.goBack);
  const config = useStore((state) => state.configs.docking);
  const setConfig = useStore((state) => state.setConfig);
  const isBackendReady = useStore((state) => state.isBackendReady);
  const showToast = useStore((state) => state.showToast);
  const showConfirm = useStore((state) => state.showConfirm);
  const setSearchActive = useStore((state) => state.setSearchActive);
  const setSearchQuery = useStore((state) => state.setSearchQuery);

  const [loading, setLoading] = useState(!config);
  const [selectedIndex, setSelectedIndex] = useState(0);
  const [editMode, setEditMode] = useState(false);
  const [editValue, setEditValue] = useState('');
  const [values, setValues] = useState<Record<string, unknown>>({});
  const [error, setError] = useState<string | null>(null);
  const [browsingField, setBrowsingField] = useState<FormField | null>(null);
  const [scrollOffset, setScrollOffset] = useState(0);
  const [isDirty, setIsDirty] = useState(false);
  useInputLock(editMode);

  const reservedRows = 12;
  const visibleRows = Math.max(2, terminalHeight - reservedRows);

  // Search functionality
  const {
    searchQuery,
    searchActive,
    filteredItems: filteredFields,
    handleSearchInput,
  } = useSearch({
    items: fields,
    searchFields: (field) => [field.label, field.key, field.section || ''],
  });

  const getValue = (field: FormField): unknown => {
    if (field.path) {
      return getByPath(values, field.path);
    }
    return values[field.key];
  };

  const setValue = (field: FormField, newValue: unknown) => {
    let newValues: Record<string, unknown>;
    if (field.path) {
      newValues = setByPath(values, field.path, newValue);
    } else {
      newValues = { ...values, [field.key]: newValue };
    }
    setValues(newValues);
    setIsDirty(true);
  };

  useEffect(() => {
    if (!config && isBackendReady) {
      loadConfig();
    } else if (config) {
      setValues(config as unknown as Record<string, unknown>);
      setLoading(false);
    }
  }, [config, isBackendReady]);

  const loadConfig = async () => {
    try {
      const bridge = getBridge();
      const data = await bridge.loadConfig('docking') as unknown as DockingConfig;
      setConfig('docking', data);
      setValues(data as unknown as Record<string, unknown>);
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
        message: 'Discard changes to docking configuration?',
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

  // Sync search state with global store for Footer display
  useEffect(() => {
    setSearchActive(searchActive);
    setSearchQuery(searchQuery);
    return () => {
      setSearchActive(false);
      setSearchQuery('');
    };
  }, [searchActive, searchQuery, setSearchActive, setSearchQuery]);

  // Reset selection when search changes
  useEffect(() => {
    setSelectedIndex(0);
    setScrollOffset(0);
  }, [searchQuery]);

  // Items to display (filtered or all)
  const displayFields = filteredFields;
  const totalFields = displayFields.length;

  const handleEditModeInput = (key: { escape?: boolean; return?: boolean }) => {
    if (key.escape) {
      setEditMode(false);
      return;
    }
    if (key.return) {
      const field = displayFields[selectedIndex];
      if (field) {
        let newValue: string | number | undefined = editValue;
        if (field.type === 'number') {
          newValue = editValue === '' ? undefined : parseFloat(editValue);
        }
        setValue(field, newValue);
      }
      setEditMode(false);
    }
  };

  const handleFieldAction = () => {
    const field = displayFields[selectedIndex];
    if (!field) return;
    const value = getValue(field);

    if (field.type === 'boolean') {
      setValue(field, !value);
    } else if (field.type === 'select' && field.options) {
      const currentIdx = field.options.indexOf(value as string);
      const nextIdx = (currentIdx + 1) % field.options.length;
      setValue(field, field.options[nextIdx]);
    } else if (field.type === 'path') {
      setBrowsingField(field);
    } else {
      setEditValue(String(value ?? ''));
      setEditMode(true);
    }
  };

  useInput((input, key) => {
    if (loading || browsingField) return;

    if (handleSearchInput(input, key)) {
      return;
    }

    if (editMode) {
      handleEditModeInput(key);
      return;
    }

    if (key.upArrow) {
      const newIndex = Math.max(0, selectedIndex - 1);
      setSelectedIndex(newIndex);
      if (newIndex < scrollOffset) {
        setScrollOffset(newIndex);
      }
    } else if (key.downArrow) {
      const newIndex = Math.min(totalFields - 1, selectedIndex + 1);
      setSelectedIndex(newIndex);
      if (newIndex >= scrollOffset + visibleRows) {
        setScrollOffset(newIndex - visibleRows + 1);
      }
    } else if (key.pageUp) {
      const newIndex = Math.max(0, selectedIndex - visibleRows);
      setSelectedIndex(newIndex);
      setScrollOffset(Math.max(0, scrollOffset - visibleRows));
    } else if (key.pageDown) {
      const newIndex = Math.min(totalFields - 1, selectedIndex + visibleRows);
      setSelectedIndex(newIndex);
      setScrollOffset(Math.min(Math.max(0, totalFields - visibleRows), scrollOffset + visibleRows));
    } else if (key.return || input === 'e') {
      handleFieldAction();
    } else if (input === 's') {
      saveConfig();
    } else if (key.escape || key.leftArrow) {
      handleExit();
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

  if (loading) {
    return (
      <AppShell
        padding={1}
        header={<Header title="Docking Config" />}
        footer={<Footer />}
      >
        <Spinner label="Loading configuration..." />
      </AppShell>
    );
  }

  if (browsingField) {
    const currentValue = String(getValue(browsingField) || process.cwd());
    const browserShortcuts = [
      { key: '↑↓', label: 'Navigate' },
      { key: 'Enter', label: 'Open/Select' },
      { key: '→/e', label: 'Edit path' },
      { key: 'Space', label: 'Search' },
      { key: 'Ctrl+F', label: 'Search' },
      { key: '←', label: 'Back' },
    ];
    return (
      <AppShell
        padding={1}
        header={<Header title="Select File" subtitle={browsingField.label} />}
        footer={<Footer shortcuts={browserShortcuts} />}
      >
        <FileBrowser
          initialPath={currentValue}
          extensions={browsingField.extensions}
          onSelect={handlePathSelect}
          onCancel={handleBrowseCancel}
        />
      </AppShell>
    );
  }

  const visibleFields = displayFields.slice(scrollOffset, scrollOffset + visibleRows);
  let lastSection = scrollOffset > 0 ? displayFields[scrollOffset - 1]?.section : undefined;
  const rowWidth = Math.max(1, terminalWidth - 4);
  const labelWidth = Math.min(20, Math.max(8, rowWidth - 6));
  const valueWidth = Math.max(1, rowWidth - labelWidth - 4);
  const selectedField = displayFields[selectedIndex];
  const selectedDescription = selectedField?.description || '';

  return (
    <AppShell
      padding={1}
      header={<Header title="Docking Config" subtitle={isDirty ? 'config_docking.yml *' : 'config_docking.yml'} />}
      footer={<Footer />}
    >
      <SearchIndicator active={searchActive} query={searchQuery} />

      {error && (
        <Box marginY={1}>
          <Text color={theme.palette.error}>Error: {error}</Text>
        </Box>
      )}

      {totalFields === 0 ? (
        <Box marginY={2}>
          <Text color={theme.palette.textMuted}>No parameters match "{searchQuery}"</Text>
        </Box>
      ) : (
        <Box flexDirection="column" marginY={1}>
          {visibleFields.map((field, index) => {
            const actualIndex = scrollOffset + index;
            const isSelected = actualIndex === selectedIndex;
            const isEditing = isSelected && editMode;
            const value = getValue(field);
            const showSection = field.section && field.section !== lastSection;
            if (field.section) lastSection = field.section;

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
                      <Text color={theme.palette.textMuted} wrap="truncate-end">
                        {field.label.padEnd(labelWidth, ' ')}
                      </Text>
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
                        <Text color={getFieldColor(theme, field.type, value)}>
                          {truncateMiddle(formatFieldValue(field.type, value), valueWidth).padEnd(valueWidth, ' ')}
                        </Text>
                      </Box>
                    )}
                  </Box>
                </Box>
              </React.Fragment>
            );
          })}
        </Box>
      )}

      <Box height={1} paddingLeft={2}>
        <Text color={theme.palette.textMuted} wrap="truncate-end">
          {selectedDescription ? `Info: ${selectedDescription}` : ' '}
        </Text>
      </Box>

      <Box>
        <Text color={theme.palette.textMuted}>
          {totalFields > 0 ? `${scrollOffset + 1}-${Math.min(scrollOffset + visibleRows, totalFields)} of ${totalFields}` : '0'} params
          {searchQuery && ` (filtered from ${fields.length})`}
        </Text>
      </Box>
    </AppShell>
  );
}

export default ConfigDocking;
