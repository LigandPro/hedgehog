import React, { useState, useEffect } from 'react';
import { Box, Text, useInput } from 'ink';
import TextInput from 'ink-text-input';
import { Header } from '../../components/Header.js';
import { Footer } from '../../components/Footer.js';
import { FileBrowser } from '../../components/FileBrowser.js';
import { useStore } from '../../store/index.js';
import { getBridge } from '../../services/python-bridge.js';
import { useInputLock } from '../../hooks/useInputLock.js';
import { useTerminalSize } from '../../hooks/useTerminalSize.js';
import { useTheme } from '../../theme/context.js';

/**
 * Truncate a path to fit within maxLen characters.
 * Keeps the end of the path visible with '...' prefix.
 */
function truncatePath(p: string, maxLen: number): string {
  if (!p || p.length <= maxLen) return p;
  return '...' + p.slice(-(maxLen - 3));
}

function getMoleculeCountErrorHint(error: unknown): string {
  const message = String(error ?? '').toLowerCase();
  if (message.includes('not found')) {
    return 'check that the file exists';
  }
  if (message.includes('permission')) {
    return 'check read permissions';
  }
  if (message.includes('not a file')) {
    return 'select a file, not a directory';
  }
  return 'check path and supported extension';
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
  const { theme } = useTheme();
  const setScreen = useStore((state) => state.setScreen);
  const goBack = useStore((state) => state.goBack);
  const mainConfig = useStore((state) => state.configs.main);
  const setConfig = useStore((state) => state.setConfig);
  const isBackendReady = useStore((state) => state.isBackendReady);
  const showToast = useStore((state) => state.showToast);

  const [focusedIndex, setFocusedIndex] = useState(0);
  const [browsingField, setBrowsingField] = useState<ConfigField | null>(null);
  const [startBrowserInPathEdit, setStartBrowserInPathEdit] = useState(false);
  const [editMode, setEditMode] = useState(false);
  const [editValue, setEditValue] = useState('');
  const [values, setValues] = useState<Record<string, unknown>>({});
  const [loading, setLoading] = useState(!mainConfig);
  const [moleculeCount, setMoleculeCount] = useState<number | null>(null);
  const [countingMolecules, setCountingMolecules] = useState(false);
  const [moleculeCountErrorHint, setMoleculeCountErrorHint] = useState<string | null>(null);
  useInputLock(editMode);

  // Get terminal size with resize support
  const { width: terminalWidth } = useTerminalSize();

  const panelWidth = Math.max(48, Math.min(96, terminalWidth - 4));
  const labelWidth = Math.max(14, Math.min(20, panelWidth - 30));
  const maxPathLen = Math.max(20, panelWidth - labelWidth - 8);

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
    const rawInputPath = values.generated_mols_path;
    const inputPath = typeof rawInputPath === 'string' ? rawInputPath.trim() : '';
    if (!inputPath || !isBackendReady) {
      setMoleculeCount(null);
      setMoleculeCountErrorHint(null);
      return;
    }

    const countMolecules = async () => {
      setCountingMolecules(true);
      setMoleculeCountErrorHint(null);
      try {
        const bridge = getBridge();
        const result = await bridge.call<{ count: number; error?: string }>('count_molecules', { path: inputPath });
        if (result.error) {
          setMoleculeCount(null);
          setMoleculeCountErrorHint(getMoleculeCountErrorHint(result.error));
          return;
        }
        setMoleculeCount(result.count);
      } catch (error) {
        setMoleculeCount(null);
        setMoleculeCountErrorHint(getMoleculeCountErrorHint(error));
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
      const validation = await bridge.validateConfig('main', values);
      if (!validation.valid) {
        showToast('error', `Invalid config: ${validation.errors.join('; ')}`);
        return;
      }
      setConfig('main', values as any);
      showToast('info', 'Input/output settings saved for this run only');
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
        setStartBrowserInPathEdit(false);
        setBrowsingField(field);
      } else if (field.type === 'boolean') {
        setValues({ ...values, [field.key]: !values[field.key] });
      } else if (field.type === 'number') {
        setEditValue(String(values[field.key] || ''));
        setEditMode(true);
      }
    } else if (key.rightArrow || input === 'e') {
      const field = CONFIG_FIELDS[focusedIndex];
      if (field.type === 'path') {
        setStartBrowserInPathEdit(true);
        setBrowsingField(field);
      }
    } else if (key.return) {
      // Enter - continue to next step
      saveAndContinue();
    } else if (key.escape || key.leftArrow) {
      goBack();
    }
  });

  const handlePathSelect = (path: string) => {
    if (browsingField) {
      setValues({ ...values, [browsingField.key]: path });
    }
    setBrowsingField(null);
    setStartBrowserInPathEdit(false);
  };

  const handleBrowseCancel = () => {
    setBrowsingField(null);
    setStartBrowserInPathEdit(false);
  };

  const trimmedInputPath =
    typeof values.generated_mols_path === 'string'
      ? values.generated_mols_path.trim()
      : '';

  if (loading) {
    return (
      <Box flexDirection="column" flexGrow={1} padding={1}>
        <Header title="Pipeline Wizard" subtitle="Loading..." />
        <Box flexGrow={1} justifyContent="center">
          <Text color={theme.palette.textMuted}>Loading configuration...</Text>
        </Box>
        <Footer />
      </Box>
    );
  }

  if (browsingField) {
    const currentValue = String(values[browsingField.key] || '') || process.cwd();
    const browserShortcuts = [
      { key: '↑↓', label: 'Navigate' },
      { key: 'Enter', label: 'Open/Select' },
      { key: '→/e', label: 'Edit path' },
      { key: 'Space', label: browsingField.isDirectory ? 'Select folder' : 'Search' },
      { key: 'Ctrl+F', label: 'Search' },
      { key: '←', label: 'Back' },
    ];
    return (
      <Box flexDirection="column" flexGrow={1} padding={1}>
        <Header title="Pipeline Wizard" subtitle={`Select ${browsingField.label}`} />
        <FileBrowser
          initialPath={currentValue}
          extensions={browsingField.extensions}
          selectDirectory={browsingField.isDirectory}
          onSelect={handlePathSelect}
          onCancel={handleBrowseCancel}
          startInPathEdit={startBrowserInPathEdit}
        />
        <Footer shortcuts={browserShortcuts} />
      </Box>
    );
  }

  return (
    <Box flexDirection="column" flexGrow={1} padding={1}>
      <Header title="Pipeline Wizard" subtitle="Step 1: Input & Output" />

      <Box flexGrow={1} justifyContent="center">
        <Box flexDirection="column" width={panelWidth}>
          <Text color={theme.palette.textMuted}>Select source molecules and where to store outputs for this run.</Text>
          <Text color={theme.palette.textMuted}>Space opens browser, Right/e edits the path directly.</Text>

          <Box flexDirection="column" marginTop={1}>
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
                ? (value ? theme.palette.primary : theme.palette.error)
                : field.type === 'boolean'
                  ? (value ? theme.palette.success : theme.palette.error)
                  : theme.palette.accent;

              return (
                <Box key={field.key}>
                  <Text color={isFocused ? theme.palette.primary : theme.palette.textMuted}>
                    {isFocused ? '▸ ' : '  '}
                  </Text>
                  <Text color={isFocused ? theme.palette.text : theme.palette.textMuted}>
                    {field.label.padEnd(labelWidth)}
                  </Text>
                  {isEditing ? (
                    <Box>
                      <TextInput
                        value={editValue}
                        onChange={setEditValue}
                        focus={true}
                      />
                    </Box>
                  ) : (
                    <Text color={valueColor}>{displayValue}</Text>
                  )}
                </Box>
              );
            })}
          </Box>

          {trimmedInputPath && (
            <Box marginTop={1}>
              <Text color={theme.palette.textMuted}>Molecules: </Text>
              {countingMolecules ? (
                <Text color={theme.palette.warning}>counting...</Text>
              ) : moleculeCount !== null ? (
                <Text color={theme.palette.success} bold>{moleculeCount.toLocaleString()}</Text>
              ) : (
                <Text color={theme.palette.error}>
                  unable to count ({moleculeCountErrorHint ?? 'check file path'})
                </Text>
              )}
            </Box>
          )}
        </Box>
      </Box>

      <Footer />
    </Box>
  );
}

export default InputSelection;
