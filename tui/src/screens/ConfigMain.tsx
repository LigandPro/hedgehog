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
import { useTheme } from '../theme/context.js';
import { useTerminalSize } from '../hooks/useTerminalSize.js';
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

function getFieldColor(
  theme: ReturnType<typeof useTheme>['theme'],
  type: ConfigFieldType,
  value: unknown,
): string {
  if (type === 'boolean') {
    return value ? theme.palette.success : theme.palette.error;
  }
  if (type === 'path') {
    return theme.palette.info;
  }
  return theme.palette.accent;
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

export function ConfigMain(): React.ReactElement {
  const { theme } = useTheme();
  const { width: terminalWidth, height: terminalHeight } = useTerminalSize();
  const setScreen = useStore((state) => state.setScreen);
  const goBack = useStore((state) => state.goBack);
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
          goBack();
        },
      });
    } else {
      goBack();
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
    } else if (key.escape || key.leftArrow) {
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

  if (loading) {
    return (
      <AppShell
        padding={1}
        header={<Header title="Main Config" />}
        footer={<Footer />}
      >
        <Spinner label="Loading configuration..." />
      </AppShell>
    );
  }

  if (browsingField) {
    const currentValue = String(values[browsingField] || process.cwd());
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
        header={<Header title="Select Path" subtitle={browsingField} />}
        footer={<Footer shortcuts={browserShortcuts} />}
      >
        <FileBrowser
          initialPath={currentValue}
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
      header={<Header title="Main Config" subtitle={isDirty ? 'config.yml *' : 'config.yml'} />}
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
          <Text color={theme.palette.textMuted}>No fields match "{searchQuery}"</Text>
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
                          {truncateMiddle(
                            field.type === 'boolean' ? (value ? 'Yes' : 'No') : String(value || '(not set)'),
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
      )}

      <Box height={1} paddingLeft={2}>
        <Text color={theme.palette.textMuted} wrap="truncate-end">
          {selectedDescription ? `Info: ${selectedDescription}` : ' '}
        </Text>
      </Box>

      <Box>
        <Text color={theme.palette.textMuted}>
          {totalFields > 0 ? `${scrollOffset + 1}-${Math.min(scrollOffset + visibleRows, totalFields)} of ${totalFields}` : '0'} fields
          {searchQuery && ` (filtered from ${fields.length})`}
        </Text>
      </Box>
    </AppShell>
  );
}

export default ConfigMain;
