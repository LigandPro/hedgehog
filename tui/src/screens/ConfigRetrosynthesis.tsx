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
import type { RetrosynthesisConfig, Screen } from '../types/index.js';

interface FormField {
  key: keyof RetrosynthesisConfig;
  label: string;
  type: 'path';
  extensions?: string[];
  section?: string;
  description?: string;
}

const fields: FormField[] = [
  // Expansion Models section
  {
    key: 'expansion_uspto_model',
    label: 'USPTO Model',
    type: 'path',
    extensions: ['.onnx'],
    section: 'Expansion Models',
    description: 'USPTO expansion model (.onnx)',
  },
  {
    key: 'expansion_uspto_templates',
    label: 'USPTO Templates',
    type: 'path',
    extensions: ['.csv.gz', '.csv'],
    description: 'USPTO reaction templates (.csv.gz)',
  },
  {
    key: 'expansion_ringbreaker_model',
    label: 'RingBreaker Model',
    type: 'path',
    extensions: ['.onnx'],
    description: 'RingBreaker expansion model (.onnx)',
  },
  {
    key: 'expansion_ringbreaker_templates',
    label: 'RingBreaker Templates',
    type: 'path',
    extensions: ['.csv.gz', '.csv'],
    description: 'RingBreaker reaction templates (.csv.gz)',
  },

  // Filter Models section
  {
    key: 'filter_uspto',
    label: 'USPTO Filter',
    type: 'path',
    extensions: ['.onnx'],
    section: 'Filter Models',
    description: 'USPTO filter model (.onnx)',
  },

  // Stock Databases section
  {
    key: 'stock_zinc',
    label: 'ZINC Stock',
    type: 'path',
    extensions: ['.hdf5', '.h5'],
    section: 'Stock Databases',
    description: 'ZINC stock database (.hdf5)',
  },
];

export function ConfigRetrosynthesis(): React.ReactElement {
  const setScreen = useStore((state) => state.setScreen);
  const previousScreen = useStore((state) => state.previousScreen);
  const config = useStore((state) => state.configs.retrosynthesis);
  const setConfig = useStore((state) => state.setConfig);
  const isBackendReady = useStore((state) => state.isBackendReady);
  const showToast = useStore((state) => state.showToast);
  const showConfirm = useStore((state) => state.showConfirm);
  const setSearchActive = useStore((state) => state.setSearchActive);
  const setSearchQuery = useStore((state) => state.setSearchQuery);

  // Determine where to go back to (Wizard or ConfigSynthesis)
  const backScreen: Screen = previousScreen === 'wizardConfigSynthesis' ? 'wizardConfigSynthesis' : 'configSynthesis';

  const [loading, setLoading] = useState(!config);
  const [selectedIndex, setSelectedIndex] = useState(0);
  const [editMode, setEditMode] = useState(false);
  const [editValue, setEditValue] = useState('');
  const [values, setValues] = useState<Partial<RetrosynthesisConfig>>(config || {});
  const [error, setError] = useState<string | null>(null);
  const [browsingField, setBrowsingField] = useState<FormField | null>(null);
  const [scrollOffset, setScrollOffset] = useState(0);
  const [isDirty, setIsDirty] = useState(false);

  const visibleRows = 12;

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
      const data = await bridge.loadConfig('retrosynthesis') as unknown as RetrosynthesisConfig;
      setConfig('retrosynthesis', data);
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
      await bridge.saveConfig('retrosynthesis', values as Record<string, unknown>);
      setConfig('retrosynthesis', values as RetrosynthesisConfig);
      setIsDirty(false);
      showToast('success', 'Configuration saved');
    } catch (err) {
      setError(String(err));
      showToast('error', `Failed to save: ${err}`);
    }
  };

  const handleExit = () => {
    if (isDirty) {
      showConfirm({
        title: 'Unsaved Changes',
        message: 'Discard changes to retrosynthesis configuration?',
        confirmLabel: 'Discard',
        cancelLabel: 'Cancel',
        onConfirm: () => {
          setScreen(backScreen);
        },
      });
    } else {
      setScreen(backScreen);
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
          setValues({ ...values, [field.key]: editValue });
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
    } else if (key.pageUp) {
      const newIndex = Math.max(0, selectedIndex - visibleRows);
      setSelectedIndex(newIndex);
      setScrollOffset(Math.max(0, scrollOffset - visibleRows));
    } else if (key.pageDown) {
      const newIndex = Math.min(totalFields - 1, selectedIndex + visibleRows);
      setSelectedIndex(newIndex);
      setScrollOffset(Math.min(Math.max(0, totalFields - visibleRows), scrollOffset + visibleRows));
    } else if (key.return || input === 'e' || input === 'b') {
      const field = displayFields[selectedIndex];
      if (!field) return;

      // All fields are paths, open file browser
      setBrowsingField(field);
    } else if (input === 's') {
      saveConfig();
    } else if (key.escape || key.leftArrow || input === 'q') {
      handleExit();
    }
  });

  const handlePathSelect = (path: string) => {
    if (browsingField) {
      setValues({ ...values, [browsingField.key]: path });
      setIsDirty(true);
    }
    setBrowsingField(null);
  };

  const handleBrowseCancel = () => {
    setBrowsingField(null);
  };

  const shortcuts = [
    { key: '↑↓', label: 'Navigate' },
    { key: 'b/Enter', label: 'Browse' },
    { key: '/', label: 'Search' },
    { key: 's', label: 'Save' },
    { key: '←/Esc', label: 'Back' },
  ];

  if (loading) {
    return (
      <Box flexDirection="column" padding={1}>
        <Header title="Retrosynthesis Config" />
        <Spinner label="Loading configuration..." />
      </Box>
    );
  }

  if (browsingField) {
    const currentValue = String(values[browsingField.key] || process.cwd());
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
      <Header title="Retrosynthesis Config" subtitle={isDirty ? 'aizynthfinder/config.yml *' : 'aizynthfinder/config.yml'} />

      <SearchIndicator active={searchActive} query={searchQuery} />

      {error && (
        <Box marginY={1}>
          <Text color="red">Error: {error}</Text>
        </Box>
      )}

      <Box marginBottom={1}>
        <Text dimColor>AiZynthFinder model and data file paths for retrosynthesis route search.</Text>
      </Box>

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
            const value = values[field.key];
            const showSection = field.section && field.section !== lastSection;
            if (field.section) lastSection = field.section;

            // Highlight matching text
            const labelParts = searchQuery ? highlightMatch(field.label) : [{ text: field.label, highlighted: false }];

            // Truncate long paths for display
            const displayValue = value ? (String(value).length > 50 ? '...' + String(value).slice(-47) : String(value)) : '(not set)';

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
                    <Box width={22}>
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
                      <Text color={value ? 'blue' : 'gray'}>
                        {displayValue}
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

export default ConfigRetrosynthesis;
