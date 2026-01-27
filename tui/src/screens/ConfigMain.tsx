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
import type { MainConfig } from '../types/index.js';

interface FormField {
  key: keyof MainConfig;
  label: string;
  type: 'text' | 'number' | 'boolean' | 'path';
  section?: string;
  description?: string;
}

const fields: FormField[] = [
  // Input/Output section
  { key: 'generated_mols_path', label: 'Input Molecules', type: 'path', section: 'Input/Output', description: 'Path to input molecules file (CSV, SMI, or SDF format)' },
  { key: 'target_mols_path', label: 'Target Molecules', type: 'path', description: 'Reference molecules for comparison (optional)' },
  { key: 'folder_to_save', label: 'Output Folder', type: 'path', description: 'Directory where all results will be saved' },
  // Processing section
  { key: 'n_jobs', label: 'Parallel Jobs', type: 'number', section: 'Processing', description: 'Number of parallel workers (-1 = all CPUs)' },
  { key: 'sample_size', label: 'Sample Size', type: 'number', description: 'Max molecules to process (use for testing)' },
  { key: 'batch_size', label: 'Batch Size', type: 'number', description: 'Molecules per batch for memory efficiency' },
  { key: 'save_sampled_mols', label: 'Save Sampled Mols', type: 'boolean', description: 'Save sampled molecules to separate file' },
  // Filter files section
  { key: 'pains_file_path', label: 'PAINS Filter File', type: 'path', section: 'Filter Files', description: 'SMARTS patterns for PAINS (pan-assay interference)' },
  { key: 'mcf_file_path', label: 'MCF File', type: 'path', description: 'Medicinal chemistry filter patterns' },
  // External tools section
  { key: 'ligand_preparation_tool', label: 'Ligand Prep Tool', type: 'path', section: 'External Tools', description: 'Path to ligand preparation binary (e.g., obabel)' },
  { key: 'protein_preparation_tool', label: 'Protein Prep Tool', type: 'path', description: 'Path to protein preparation binary' },
];

type ConfigFieldType = 'text' | 'number' | 'boolean' | 'path';

function getFieldColor(type: ConfigFieldType, value: unknown): 'green' | 'red' | 'blue' | 'yellow' {
  if (type === 'boolean') {
    return value ? 'green' : 'red';
  }
  if (type === 'path') {
    return 'blue';
  }
  return 'yellow';
}

export function ConfigMain(): React.ReactElement {
  const setScreen = useStore((state) => state.setScreen);
  const config = useStore((state) => state.configs.main);
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
  const [values, setValues] = useState<Partial<MainConfig>>(config || {});
  const [error, setError] = useState<string | null>(null);
  const [browsingField, setBrowsingField] = useState<keyof MainConfig | null>(null);
  const [scrollOffset, setScrollOffset] = useState(0);
  const [isDirty, setIsDirty] = useState(false);

  const visibleRows = 14;

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

  // Sync search state with global store
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
        message: 'Discard changes to main configuration?',
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
          let newValue: string | number | boolean = editValue;

          if (field.type === 'number') {
            newValue = parseInt(editValue, 10) || 0;
          }

          setValues({ ...values, [field.key]: newValue });
          setIsDirty(true);
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
    } else if (key.return || input === 'e') {
      const field = displayFields[selectedIndex];
      if (!field) return;
      if (field.type === 'boolean') {
        setValues({ ...values, [field.key]: !values[field.key] });
        setIsDirty(true);
      } else if (field.type === 'path') {
        setBrowsingField(field.key);
      } else {
        setEditValue(String(values[field.key] || ''));
        setEditMode(true);
      }
    } else if (input === 'b') {
      const field = displayFields[selectedIndex];
      if (field?.type === 'path') {
        setBrowsingField(field.key);
      }
    } else if (input === 's') {
      saveConfig();
    } else if (key.escape || key.leftArrow || input === 'q') {
      handleExit();
    }
  });

  const handlePathSelect = (path: string) => {
    if (browsingField) {
      setValues({ ...values, [browsingField]: path });
      setIsDirty(true);
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
    { key: '←/Esc', label: 'Back' },
  ];

  if (loading) {
    return (
      <Box flexDirection="column" padding={1}>
        <Header title="Main Config" />
        <Spinner label="Loading configuration..." />
      </Box>
    );
  }

  if (browsingField) {
    const currentValue = String(values[browsingField] || process.cwd());
    return (
      <Box flexDirection="column" padding={1}>
        <Header title="Select Path" subtitle={browsingField} />
        <FileBrowser
          initialPath={currentValue}
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
      <Header title="Main Config" subtitle={isDirty ? 'config.yml *' : 'config.yml'} />

      <SearchIndicator active={searchActive} query={searchQuery} />

      {error && (
        <Box marginY={1}>
          <Text color="red">Error: {error}</Text>
        </Box>
      )}

      {totalFields === 0 ? (
        <Box marginY={2}>
          <Text dimColor>No fields match "{searchQuery}"</Text>
        </Box>
      ) : (
        <Box flexDirection="column" marginY={1}>
          {visibleFields.map((field, index) => {
            const actualIndex = scrollOffset + index;
            const isSelected = actualIndex === selectedIndex;
            const isEditing = isSelected && editMode;
            const value = values[field.key];
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
                        {field.type === 'boolean' ? (value ? 'Yes' : 'No') : String(value || '(not set)')}
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
          {totalFields > 0 ? `${scrollOffset + 1}-${Math.min(scrollOffset + visibleRows, totalFields)} of ${totalFields}` : '0'} fields
          {searchQuery && ` (filtered from ${fields.length})`}
        </Text>
      </Box>

      <Footer shortcuts={shortcuts} />
    </Box>
  );
}

export default ConfigMain;
