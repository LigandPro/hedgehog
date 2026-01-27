import React, { useState, useEffect } from 'react';
import { Box, Text, useInput } from 'ink';
import TextInput from 'ink-text-input';
import { Header } from '../components/Header.js';
import { Footer } from '../components/Footer.js';
import { Spinner } from '../components/Spinner.js';
import { FileBrowser } from '../components/FileBrowser.js';
import { SearchIndicator } from '../components/SearchIndicator.js';
import { useStore } from '../store/index.js';
import { getBridge } from '../services/python-bridge.js';
import { useSearch } from '../hooks/useSearch.js';
import { getByPath, setByPath } from '../utils/object-path.js';
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
  { key: 'tools', label: 'Tools', type: 'select', options: ['both', 'smina', 'gnina'], description: 'Docking engine: smina (CPU), gnina (GPU+CNN), or both' },
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

function getFieldColor(type: FieldType, value: unknown): 'green' | 'red' | 'cyan' | 'blue' | 'yellow' {
  if (type === 'boolean') {
    return value ? 'green' : 'red';
  }
  if (type === 'select') {
    return 'cyan';
  }
  if (type === 'path') {
    return 'blue';
  }
  return 'yellow';
}

export function ConfigDocking(): React.ReactElement {
  const setScreen = useStore((state) => state.setScreen);
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
  const [values, setValues] = useState<Record<string, any>>({});
  const [error, setError] = useState<string | null>(null);
  const [browsingField, setBrowsingField] = useState<FormField | null>(null);
  const [scrollOffset, setScrollOffset] = useState(0);
  const [isDirty, setIsDirty] = useState(false);

  const visibleRows = 16;

  // Search functionality
  const {
    searchQuery,
    searchActive,
    filteredItems: filteredFields,
    handleSearchInput,
    highlightMatch,
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
          setScreen('welcome');
        },
      });
    } else {
      setScreen('welcome');
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

  useInput((input, key) => {
    if (loading || browsingField) return;

    // Handle search input first
    if (handleSearchInput(input, key)) {
      return;
    }

    if (editMode) {
      if (key.escape) {
        setEditMode(false);
      } else if (key.return) {
        const field = displayFields[selectedIndex];
        if (field) {
          let newValue: any = editValue;

          if (field.type === 'number') {
            newValue = editValue === '' ? undefined : parseFloat(editValue);
          }

          setValue(field, newValue);
        }
        setEditMode(false);
      }
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
    } else if (input === 's') {
      saveConfig();
    } else if (key.escape || key.leftArrow || input === 'q') {
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

  const shortcuts = [
    { key: '↑↓', label: 'Navigate' },
    { key: 'PgUp/Dn', label: 'Page' },
    { key: 'e/Enter', label: 'Edit' },
    { key: 's', label: 'Save' },
    { key: '←/Esc', label: 'Back' },
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

  const visibleFields = displayFields.slice(scrollOffset, scrollOffset + visibleRows);
  let lastSection = scrollOffset > 0 ? displayFields[scrollOffset - 1]?.section : undefined;

  return (
    <Box flexDirection="column" padding={1}>
      <Header title="Docking Config" subtitle={isDirty ? 'config_docking.yml *' : 'config_docking.yml'} />

      <SearchIndicator active={searchActive} query={searchQuery} />

      {error && (
        <Box marginY={1}>
          <Text color="red">Error: {error}</Text>
        </Box>
      )}

      {totalFields === 0 ? (
        <Box marginY={2}>
          <Text dimColor>No parameters match "{searchQuery}"</Text>
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

            // Highlight matching text
            const labelParts = searchQuery ? highlightMatch(field.label) : [{ text: field.label, highlighted: false }];

            return (
              <React.Fragment key={field.key}>
                {showSection && (
                  <Box marginTop={index > 0 ? 1 : 0}>
                    <Text color="cyan" bold>─ {field.section} ─</Text>
                  </Box>
                )}
                <Box flexDirection="column">
                  <Box>
                    <Text color={isSelected ? 'cyan' : 'white'}>
                      {isSelected ? '▶ ' : '  '}
                    </Text>
                    <Box width={20}>
                      {labelParts.map((part, pi) => (
                        <Text key={pi} dimColor={!part.highlighted} color={part.highlighted ? 'yellow' : undefined} bold={part.highlighted}>
                          {part.text}
                        </Text>
                      ))}
                    </Box>
                    {isEditing ? (
                      <Box>
                        <TextInput
                          value={editValue}
                          onChange={setEditValue}
                          focus={true}
                        />
                      </Box>
                    ) : (
                      <Text color={getFieldColor(field.type, value)}>
                        {field.type === 'boolean'
                          ? (value ? 'Yes' : 'No')
                          : value !== undefined && value !== null
                            ? String(value)
                            : '(not set)'}
                      </Text>
                    )}
                  </Box>
                  {isSelected && field.description && (
                    <Box paddingLeft={4}>
                      <Text dimColor italic>{field.description}</Text>
                    </Box>
                  )}
                </Box>
              </React.Fragment>
            );
          })}
        </Box>
      )}

      <Box>
        <Text dimColor>
          {totalFields > 0 ? `${scrollOffset + 1}-${Math.min(scrollOffset + visibleRows, totalFields)} of ${totalFields}` : '0'} params
          {searchQuery && ` (filtered from ${fields.length})`}
        </Text>
      </Box>

      <Footer shortcuts={shortcuts} />
    </Box>
  );
}

export default ConfigDocking;
