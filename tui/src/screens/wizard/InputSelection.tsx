import React, { useState, useEffect } from 'react';
import { Box, Text, useInput } from 'ink';
import TextInput from 'ink-text-input';
import { Header } from '../../components/Header.js';
import { Footer } from '../../components/Footer.js';
import { FileBrowser } from '../../components/FileBrowser.js';
import { useStore } from '../../store/index.js';
import { getBridge } from '../../services/python-bridge.js';
import { useTerminalSize } from '../../hooks/useTerminalSize.js';

/**
 * Truncate a path to fit within maxLen characters.
 * Keeps the end of the path visible with '...' prefix.
 */
function truncatePath(p: string, maxLen: number): string {
  if (!p || p.length <= maxLen) return p;
  return '...' + p.slice(-(maxLen - 3));
}

interface ConfigField {
  key: string;
  label: string;
  type: 'path' | 'number' | 'boolean';
  isDirectory?: boolean;
  extensions?: string[];
}

const CONFIG_FIELDS: ConfigField[] = [
  {
    key: 'generated_mols_path',
    label: 'Input File',
    type: 'path',
    isDirectory: false,
    extensions: ['.csv', '.smi', '.smiles', '.sdf', '.mol2'],
  },
  {
    key: 'folder_to_save',
    label: 'Output Folder',
    type: 'path',
    isDirectory: true,
  },
];

export function InputSelection(): React.ReactElement {
  const setScreen = useStore((state) => state.setScreen);
  const mainConfig = useStore((state) => state.configs.main);
  const setConfig = useStore((state) => state.setConfig);
  const isBackendReady = useStore((state) => state.isBackendReady);
  const showToast = useStore((state) => state.showToast);

  const [focusedIndex, setFocusedIndex] = useState(0);
  const [browsingField, setBrowsingField] = useState<ConfigField | null>(null);
  const [editMode, setEditMode] = useState(false);
  const [editValue, setEditValue] = useState('');
  const [values, setValues] = useState<Record<string, unknown>>({});
  const [loading, setLoading] = useState(!mainConfig);
  const [moleculeCount, setMoleculeCount] = useState<number | null>(null);
  const [countingMolecules, setCountingMolecules] = useState(false);

  // Get terminal size with resize support
  const { width: terminalWidth } = useTerminalSize();

  // Reserve: indicator (2) + label (20) + padding (2) + margin (2)
  const maxPathLen = Math.max(20, terminalWidth - 26);

  useEffect(() => {
    if (!mainConfig && isBackendReady) {
      loadConfig();
    } else if (mainConfig) {
      setValues(mainConfig as unknown as Record<string, unknown>);
      setLoading(false);
    }
  }, [mainConfig, isBackendReady]);

  // Count molecules when input file changes
  useEffect(() => {
    const inputPath = values.generated_mols_path as string;
    if (!inputPath || !isBackendReady) {
      setMoleculeCount(null);
      return;
    }

    const countMolecules = async () => {
      setCountingMolecules(true);
      try {
        const bridge = getBridge();
        const result = await bridge.call<{ count: number }>('count_molecules', { path: inputPath });
        setMoleculeCount(result.count);
      } catch {
        setMoleculeCount(null);
      } finally {
        setCountingMolecules(false);
      }
    };

    countMolecules();
  }, [values.generated_mols_path, isBackendReady]);

  const loadConfig = async () => {
    try {
      const bridge = getBridge();
      const data = await bridge.loadConfig('main');
      setConfig('main', data as any);
      setValues(data as Record<string, unknown>);
    } catch (err) {
      showToast('error', `Failed to load config: ${err}`);
    } finally {
      setLoading(false);
    }
  };

  const saveAndContinue = async () => {
    if (!values.generated_mols_path) {
      showToast('warning', 'Please select an input file');
      return;
    }

    try {
      const bridge = getBridge();
      await bridge.saveConfig('main', values);
      setConfig('main', values as any);
      setScreen('wizardStageSelection');
    } catch (err) {
      showToast('error', `Failed to save: ${err}`);
    }
  };

  useInput((input, key) => {
    if (loading || browsingField) return;

    // Edit mode for number fields
    if (editMode) {
      if (key.escape) {
        setEditMode(false);
      } else if (key.return) {
        const field = CONFIG_FIELDS[focusedIndex];
        const numVal = parseInt(editValue, 10);
        if (!isNaN(numVal)) {
          setValues({ ...values, [field.key]: numVal });
        }
        setEditMode(false);
      }
      return;
    }

    if (key.upArrow) {
      setFocusedIndex(Math.max(0, focusedIndex - 1));
    } else if (key.downArrow) {
      setFocusedIndex(Math.min(CONFIG_FIELDS.length - 1, focusedIndex + 1));
    } else if (key.tab && !key.shift) {
      setFocusedIndex((focusedIndex + 1) % CONFIG_FIELDS.length);
    } else if (key.tab && key.shift) {
      setFocusedIndex((focusedIndex - 1 + CONFIG_FIELDS.length) % CONFIG_FIELDS.length);
    } else if (input === ' ') {
      // Space - edit/browse current field
      const field = CONFIG_FIELDS[focusedIndex];
      if (field.type === 'path') {
        setBrowsingField(field);
      } else if (field.type === 'boolean') {
        setValues({ ...values, [field.key]: !values[field.key] });
      } else if (field.type === 'number') {
        setEditValue(String(values[field.key] || ''));
        setEditMode(true);
      }
    } else if (key.return) {
      // Enter - continue to next step
      saveAndContinue();
    } else if (key.escape || key.leftArrow || input === 'q') {
      setScreen('welcome');
    }
  });

  const handlePathSelect = (path: string) => {
    if (browsingField) {
      setValues({ ...values, [browsingField.key]: path });
    }
    setBrowsingField(null);
  };

  const handleBrowseCancel = () => {
    setBrowsingField(null);
  };

  if (loading) {
    return (
      <Box flexDirection="column" padding={1}>
        <Header title="Pipeline Wizard" subtitle="Loading..." />
        <Text dimColor>Loading configuration...</Text>
      </Box>
    );
  }

  if (browsingField) {
    const currentValue = String(values[browsingField.key] || '') || process.cwd();
    const browserShortcuts = [
      { key: '↑↓', label: 'Navigate' },
      { key: 'Enter', label: browsingField.isDirectory ? 'Open' : 'Select' },
      ...(browsingField.isDirectory ? [{ key: 'Space', label: 'Select folder' }] : []),
      { key: '/', label: 'Search' },
      { key: '←', label: 'Back' },
    ];
    return (
      <Box flexDirection="column" padding={1}>
        <Header title="Pipeline Wizard" subtitle={`Select ${browsingField.label}`} />
        <FileBrowser
          initialPath={currentValue}
          extensions={browsingField.extensions}
          selectDirectory={browsingField.isDirectory}
          onSelect={handlePathSelect}
          onCancel={handleBrowseCancel}
        />
        <Footer shortcuts={browserShortcuts} />
      </Box>
    );
  }

  const shortcuts = [
    { key: '↑↓', label: 'Navigate' },
    { key: 'Space', label: 'Edit' },
    { key: 'Enter', label: 'Next' },
    { key: '←/Esc', label: 'Back' },
  ];

  return (
    <Box flexDirection="column" padding={1}>
      <Header title="Pipeline Wizard" subtitle="Step 1: Input & Output" />

      <Box flexDirection="column" marginY={1}>
        {CONFIG_FIELDS.map((field, index) => {
          const isFocused = focusedIndex === index;
          const value = values[field.key];
          const isEditing = isFocused && editMode && field.type === 'number';

          let displayValue = '';
          if (field.type === 'path') {
            displayValue = value ? truncatePath(String(value), maxPathLen) : '(not set)';
          } else if (field.type === 'boolean') {
            displayValue = value ? 'Yes' : 'No';
          } else if (field.type === 'number') {
            displayValue = String(value ?? '');
          }

          const valueColor = field.type === 'path'
            ? (value ? 'cyan' : 'red')
            : field.type === 'boolean'
              ? (value ? 'green' : 'red')
              : 'yellow';

          return (
            <Box key={field.key}>
              <Text color={isFocused ? 'cyan' : 'gray'}>
                {isFocused ? '▸ ' : '  '}
              </Text>
              <Text color={isFocused ? 'white' : 'gray'}>{field.label.padEnd(20)}</Text>
              {isEditing ? (
                <Box>
                  <TextInput
                    value={editValue}
                    onChange={setEditValue}
                    focus={true}
                  />
                </Box>
              ) : (
                <Text color={valueColor}>
                  {displayValue}
                </Text>
              )}
            </Box>
          );
        })}
      </Box>

      {/* Molecule count */}
      {(values.generated_mols_path as string) && (
        <Box marginY={1}>
          <Text dimColor>Molecules: </Text>
          {countingMolecules ? (
            <Text color="yellow">counting...</Text>
          ) : moleculeCount !== null ? (
            <Text color="green" bold>{moleculeCount.toLocaleString()}</Text>
          ) : (
            <Text color="red">unable to count</Text>
          )}
        </Box>
      )}

      <Footer shortcuts={shortcuts} />
    </Box>
  );
}

export default InputSelection;
