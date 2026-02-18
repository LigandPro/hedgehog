import React, { useState, useEffect } from 'react';
import { Box, Text, useInput } from 'ink';
import TextInput from 'ink-text-input';
import { Header } from '../../components/Header.js';
import { Footer } from '../../components/Footer.js';
import { Spinner } from '../../components/Spinner.js';
import { useStore } from '../../store/index.js';
import { getBridge } from '../../services/python-bridge.js';
import { useTerminalSize } from '../../hooks/useTerminalSize.js';

interface ParamDef {
  key: string;
  label: string;
  type: 'boolean' | 'number' | 'range' | 'select' | 'text';
  minKey?: string;
  maxKey?: string;
  options?: string[];
  configPath?: string[]; // Path to nested config value
  description?: string;
}

interface StageConfigDef {
  stageName: string;
  displayName: string;
  configType: 'mol_prep' | 'descriptors' | 'filters' | 'synthesis' | 'docking' | 'docking_filters';
  params: ParamDef[];
}

const MOL_PREP_PARAMS: ParamDef[] = [
  { key: 'run', label: 'Run Stage', type: 'boolean', description: 'Enable/disable Datamol molecule preparation' },
  { key: 'n_jobs', label: 'Parallel Jobs', type: 'number', description: '-1 uses all available CPU cores' },
  { key: 'enabled', label: 'Fix Molecule', type: 'boolean', configPath: ['steps', 'fix_mol'], description: 'Apply Datamol fix_mol cleanup' },
  { key: 'enabled', label: 'Remove Salts', type: 'boolean', configPath: ['steps', 'remove_salts_solvents'], description: 'Remove salts and solvents from structures' },
  { key: 'keep_largest_fragment', label: 'Largest Fragment', type: 'boolean', configPath: ['steps'], description: 'Keep only the largest molecular fragment' },
  { key: 'enabled', label: 'Standardize Mol', type: 'boolean', configPath: ['steps', 'standardize_mol'], description: 'Apply normalization/reionization/uncharging pipeline' },
  { key: 'uncharge', label: 'Uncharge', type: 'boolean', configPath: ['steps', 'standardize_mol'], description: 'Convert charged species to neutral when possible' },
  { key: 'stereo', label: 'Standardize Stereo', type: 'boolean', configPath: ['steps', 'standardize_mol'], description: 'Canonicalize stereochemistry during standardization' },
  { key: 'remove_stereochemistry', label: 'Drop Stereo', type: 'boolean', configPath: ['steps'], description: 'Remove stereochemistry in final standardized molecules' },
  { key: 'allowed_atoms', label: 'Allowed Atoms', type: 'text', configPath: ['filters'], description: 'Comma-separated whitelist of allowed atom symbols' },
  { key: 'require_neutral', label: 'Require Neutral', type: 'boolean', configPath: ['filters'], description: 'Reject molecules with formal charge after prep' },
  { key: 'reject_radicals', label: 'Reject Radicals', type: 'boolean', configPath: ['filters'], description: 'Reject molecules containing radicals' },
  { key: 'reject_isotopes', label: 'Reject Isotopes', type: 'boolean', configPath: ['filters'], description: 'Reject isotopically labelled molecules' },
  { key: 'require_single_fragment', label: 'Single Fragment', type: 'boolean', configPath: ['filters'], description: 'Reject molecules with multiple fragments' },
];

// Descriptors config params - all from config_descriptors.yml
const DESCRIPTORS_PARAMS: ParamDef[] = [
  { key: 'run', label: 'Run Stage', type: 'boolean', description: 'Enable/disable descriptors calculation' },
  { key: 'batch_size', label: 'Batch Size', type: 'number', description: 'Molecules per batch for memory efficiency' },
  { key: 'filter_data', label: 'Filter Data', type: 'boolean', description: 'Filter molecules outside border ranges' },
  { key: 'allowed_chars', label: 'Allowed Atoms', type: 'text', configPath: ['borders'], description: 'Allowed atom types (C, N, O, S, F, Cl, Br)' },
  { key: 'charged_mol_allowed', label: 'Charged Allowed', type: 'boolean', configPath: ['borders'], description: 'Allow charged molecules' },
  // Descriptor ranges
  { key: 'molWt', label: 'Molecular Weight', type: 'range', minKey: 'molWt_min', maxKey: 'molWt_max', configPath: ['borders'], description: 'Molecular weight range (Da)' },
  { key: 'logP', label: 'logP', type: 'range', minKey: 'logP_min', maxKey: 'logP_max', configPath: ['borders'], description: 'Lipophilicity (partition coefficient)' },
  { key: 'hbd', label: 'H-Bond Donors', type: 'range', minKey: 'hbd_min', maxKey: 'hbd_max', configPath: ['borders'], description: 'Hydrogen bond donors count' },
  { key: 'hba', label: 'H-Bond Acceptors', type: 'range', minKey: 'hba_min', maxKey: 'hba_max', configPath: ['borders'], description: 'Hydrogen bond acceptors count' },
  { key: 'tpsa', label: 'TPSA', type: 'range', minKey: 'tpsa_min', maxKey: 'tpsa_max', configPath: ['borders'], description: 'Topological polar surface area (Å²)' },
  { key: 'n_rot_bonds', label: 'Rotatable Bonds', type: 'range', minKey: 'n_rot_bonds_min', maxKey: 'n_rot_bonds_max', configPath: ['borders'], description: 'Flexibility indicator' },
  { key: 'n_rings', label: 'Number of Rings', type: 'range', minKey: 'n_rings_min', maxKey: 'n_rings_max', configPath: ['borders'], description: 'Ring systems count' },
  { key: 'fsp3', label: 'Fraction SP3', type: 'range', minKey: 'fsp3_min', maxKey: 'fsp3_max', configPath: ['borders'], description: 'Saturation: sp3 carbons / total C' },
  { key: 'qed', label: 'QED', type: 'range', minKey: 'qed_min', maxKey: 'qed_max', configPath: ['borders'], description: 'Quantitative drug-likeness (0-1)' },
];

const FILTERS_PARAMS: ParamDef[] = [
  { key: 'run', label: 'Run Stage', type: 'boolean', description: 'Enable/disable structural filters' },
  { key: 'calculate_common_alerts', label: 'Common Alerts', type: 'boolean', description: 'PAINS, Dundee, Glaxo alerts' },
  { key: 'calculate_NIBR', label: 'NIBR Filters', type: 'boolean', description: 'Novartis structural filters' },
  { key: 'calculate_lilly', label: 'Lilly Filters', type: 'boolean', description: 'Eli Lilly medchem demerits' },
  { key: 'filter_data', label: 'Filter Data', type: 'boolean', description: 'Remove flagged molecules' },
];

const SYNTHESIS_PARAMS: ParamDef[] = [
  { key: 'run', label: 'Run Stage', type: 'boolean', description: 'Enable/disable synthesis scoring' },
  { key: 'run_retrosynthesis', label: 'Run Retrosynthesis', type: 'boolean', description: 'Enable AiZynthFinder route search' },
  { key: 'filter_solved_only', label: 'Filter Solved Only', type: 'boolean', description: 'Keep only molecules with found routes' },
  { key: 'sa_score', label: 'SA Score', type: 'range', minKey: 'sa_score_min', maxKey: 'sa_score_max', description: 'Synthetic Accessibility (1-10, lower=easier)' },
  { key: 'syba_score', label: 'SYBA Score', type: 'range', minKey: 'syba_score_min', maxKey: 'syba_score_max', description: 'Synthetic Bayesian Accessibility' },
  { key: 'ra_score', label: 'RA Score', type: 'range', minKey: 'ra_score_min', maxKey: 'ra_score_max', description: 'Retrosynthesis probability (0-1)' },
];

const DOCKING_PARAMS: ParamDef[] = [
  { key: 'run', label: 'Run Stage', type: 'boolean', description: 'Enable/disable docking' },
  { key: 'tools', label: 'Tool', type: 'select', options: ['smina', 'gnina', 'both'], description: 'smina (CPU), gnina (GPU+CNN)' },
  { key: 'exhaustiveness', label: 'Exhaustiveness', type: 'number', configPath: ['smina_config'], description: 'Search depth (8=fast, 32=thorough)' },
  { key: 'num_modes', label: 'Num Modes', type: 'number', configPath: ['smina_config'], description: 'Max binding poses to generate' },
  { key: 'energy_range', label: 'Energy Range', type: 'number', configPath: ['smina_config'], description: 'Max energy diff from best (kcal/mol)' },
];

const DOCKING_FILTERS_PARAMS: ParamDef[] = [
  { key: 'run', label: 'Run Stage', type: 'boolean', description: 'Enable/disable docking filters' },
  { key: 'enabled', label: 'Pose Quality', type: 'boolean', configPath: ['pose_quality'], description: 'PoseCheck: clashes and strain filters' },
  { key: 'max_clashes', label: 'Max Clashes', type: 'number', configPath: ['pose_quality'], description: 'Maximum allowed steric clashes' },
  { key: 'enabled', label: 'Interactions', type: 'boolean', configPath: ['interactions'], description: 'ProLIF: interaction-based filtering' },
  { key: 'min_hbonds', label: 'Min H-Bonds', type: 'number', configPath: ['interactions'], description: 'Minimum required hydrogen bonds' },
  { key: 'enabled', label: 'Conformer Deviation', type: 'boolean', configPath: ['conformer_deviation'], description: 'Reject poses far from plausible conformers' },
  { key: 'max_rmsd_to_conformer', label: 'Max RMSD', type: 'number', configPath: ['conformer_deviation'], description: 'Max RMSD to nearest conformer (Å)' },
  { key: 'mode', label: 'Aggregation Mode', type: 'select', options: ['all', 'any'], configPath: ['aggregation'], description: 'all=pass all enabled filters, any=pass at least one' },
  { key: 'save_failed', label: 'Save Failed', type: 'boolean', configPath: ['aggregation'], description: 'Save failed molecules CSV' },
];

const STAGE_CONFIGS: Record<string, StageConfigDef> = {
  mol_prep: {
    stageName: 'mol_prep',
    displayName: 'Mol Prep',
    configType: 'mol_prep',
    params: MOL_PREP_PARAMS,
  },
  descriptors: {
    stageName: 'descriptors',
    displayName: 'Descriptors',
    configType: 'descriptors',
    params: DESCRIPTORS_PARAMS,
  },
  struct_filters: {
    stageName: 'struct_filters',
    displayName: 'Struct Filters',
    configType: 'filters',
    params: FILTERS_PARAMS,
  },
  synthesis: {
    stageName: 'synthesis',
    displayName: 'Synthesis',
    configType: 'synthesis',
    params: SYNTHESIS_PARAMS,
  },
  docking: {
    stageName: 'docking',
    displayName: 'Docking',
    configType: 'docking',
    params: DOCKING_PARAMS,
  },
  docking_filters: {
    stageName: 'docking_filters',
    displayName: 'Docking Filters',
    configType: 'docking_filters',
    params: DOCKING_FILTERS_PARAMS,
  },
};

interface QuickConfigProps {
  stageName: string;
}

interface RangeValue {
  min: unknown;
  max: unknown;
}

type ParamType = 'boolean' | 'number' | 'range' | 'select' | 'text';

function getParamColor(type: ParamType, value: unknown): 'green' | 'red' | 'cyan' | 'yellow' {
  if (type === 'boolean') {
    return value ? 'green' : 'red';
  }
  if (type === 'select') {
    return 'cyan';
  }
  return 'yellow';
}

export function QuickConfig({ stageName }: QuickConfigProps): React.ReactElement {
  const setScreen = useStore((state) => state.setScreen);
  const showToast = useStore((state) => state.showToast);
  const setConfig = useStore((state) => state.setConfig);
  const configs = useStore((state) => state.configs);

  const [focusedIndex, setFocusedIndex] = useState(0);
  const [editMode, setEditMode] = useState<'none' | 'value' | 'min' | 'max'>('none');
  const [editValue, setEditValue] = useState('');
  const [loading, setLoading] = useState(true);
  const [rawConfig, setRawConfig] = useState<Record<string, unknown> | null>(null);
  const [isDirty, setIsDirty] = useState(false);

  const { width: terminalWidth, height: terminalHeight } = useTerminalSize();
  const listHeight = Math.max(5, terminalHeight - 12);
  const descriptionWidth = Math.min(40, Math.max(24, Math.floor(terminalWidth * 0.35)));
  const showSideDescription = terminalWidth >= 80;
  const contentWidth = Math.max(20, terminalWidth - 2);

  const stageConfig = STAGE_CONFIGS[stageName];

  // Load config on mount
  useEffect(() => {
    loadConfig();
  }, [stageName]);

  const loadConfig = async () => {
    try {
      const bridge = getBridge();
      const configType = stageConfig.configType;
      let cfg = configs[configType] as Record<string, unknown> | null;
      if (!cfg) {
        cfg = await bridge.loadConfig(configType) as Record<string, unknown>;
        setConfig(configType, cfg as any);
      }
      setRawConfig(cfg);
    } catch (err) {
      showToast('error', `Failed to load config: ${err}`);
    } finally {
      setLoading(false);
    }
  };

  const getValue = (param: ParamDef): unknown => {
    if (!rawConfig) return undefined;

    if (param.type === 'range' && param.minKey && param.maxKey) {
      if (param.configPath) {
        let target = rawConfig as Record<string, unknown>;
        for (const key of param.configPath) {
          target = (target[key] || {}) as Record<string, unknown>;
        }
        return { min: target[param.minKey], max: target[param.maxKey] };
      }
      return { min: rawConfig[param.minKey], max: rawConfig[param.maxKey] };
    }

    if (param.configPath) {
      let target = rawConfig as Record<string, unknown>;
      for (const key of param.configPath) {
        target = (target[key] || {}) as Record<string, unknown>;
      }
      return target[param.key];
    }

    return rawConfig[param.key];
  };

  const setValue = (param: ParamDef, value: unknown, isMin?: boolean) => {
    if (!rawConfig) return;

    const newConfig = JSON.parse(JSON.stringify(rawConfig));

    if (param.type === 'range' && param.minKey && param.maxKey) {
      const key = isMin ? param.minKey : param.maxKey;
      if (param.configPath) {
        let target = newConfig;
        for (const k of param.configPath) {
          if (!target[k]) target[k] = {};
          target = target[k];
        }
        target[key] = value;
      } else {
        newConfig[key] = value;
      }
    } else if (param.configPath) {
      let target = newConfig;
      for (const k of param.configPath) {
        if (!target[k]) target[k] = {};
        target = target[k];
      }
      target[param.key] = value;
    } else {
      newConfig[param.key] = value;
    }

    setRawConfig(newConfig);
    setIsDirty(true);
  };

  const formatValue = (param: ParamDef): string => {
    const value = getValue(param);
    if (value === undefined) return '-';

    if (param.type === 'boolean') {
      return value ? 'Yes' : 'No';
    }
    if (param.type === 'range') {
      const { min, max } = value as { min: unknown; max: unknown };
      return `${min ?? '-'} - ${max ?? '-'}`;
    }
    if (param.type === 'text' && Array.isArray(value)) {
      return value.join(', ');
    }
    return String(value ?? '');
  };

  const saveConfig = async () => {
    if (!rawConfig) return;

    try {
      const bridge = getBridge();
      await bridge.saveConfig(stageConfig.configType, rawConfig);
      setConfig(stageConfig.configType, rawConfig as any);
      setIsDirty(false);
      // Return to stage selection after save
      setScreen('wizardStageSelection');
    } catch (err) {
      showToast('error', `Failed to save: ${err}`);
    }
  };

  const goBack = () => {
    setScreen('wizardStageSelection');
  };

  useInput((input, key) => {
    if (loading) return;

    if (editMode !== 'none') {
      if (key.escape) {
        setEditMode('none');
      } else if (key.return) {
        const param = stageConfig.params[focusedIndex];
        if (param.type === 'range') {
          const numValue = parseFloat(editValue);
          if (!isNaN(numValue) || editValue === 'inf') {
            setValue(param, editValue === 'inf' ? 'inf' : numValue, editMode === 'min');
          }
        } else if (param.type === 'number') {
          const numValue = parseFloat(editValue);
          if (!isNaN(numValue)) {
            setValue(param, numValue);
          }
        } else if (param.type === 'text') {
          // For array values like allowed_chars
          const value = getValue(param);
          if (Array.isArray(value)) {
            setValue(param, editValue.split(',').map(s => s.trim()));
          } else {
            setValue(param, editValue);
          }
        }
        setEditMode('none');
      } else if (key.tab && editMode === 'min') {
        // Move to max editing
        const param = stageConfig.params[focusedIndex];
        if (param.type === 'range') {
          const numValue = parseFloat(editValue);
          if (!isNaN(numValue)) {
            setValue(param, numValue, true);
          }
          const val = getValue(param) as { max: unknown };
          setEditValue(String(val?.max ?? ''));
          setEditMode('max');
        }
      }
      return;
    }

    if (key.upArrow) {
      setFocusedIndex(Math.max(0, focusedIndex - 1));
    } else if (key.downArrow) {
      setFocusedIndex(Math.min(stageConfig.params.length - 1, focusedIndex + 1));
    } else if (input === ' ') {
      // Space - edit current field
      const param = stageConfig.params[focusedIndex];
      if (param.type === 'boolean') {
        setValue(param, !getValue(param));
      } else if (param.type === 'select' && param.options) {
        const current = getValue(param);
        const idx = param.options.indexOf(current as string);
        const nextIdx = (idx + 1) % param.options.length;
        setValue(param, param.options[nextIdx]);
      } else if (param.type === 'range') {
        const val = getValue(param) as { min: unknown };
        setEditValue(String(val?.min ?? ''));
        setEditMode('min');
      } else if (param.type === 'number') {
        setEditValue(String(getValue(param) ?? ''));
        setEditMode('value');
      } else if (param.type === 'text') {
        const val = getValue(param);
        setEditValue(Array.isArray(val) ? val.join(', ') : String(val ?? ''));
        setEditMode('value');
      }
    } else if (key.return || input === 's') {
      // Enter or 's' - save and return to stage selection
      saveConfig();
    } else if (input === 'r' && stageName === 'synthesis') {
      // 'r' - go to detailed retrosynthesis config (only for synthesis stage)
      setScreen('configRetrosynthesis');
    } else if (key.leftArrow || key.escape || input === 'q') {
      goBack();
    }
  });

  // Calculate visible range
  const startIdx = Math.max(0, Math.min(focusedIndex - Math.floor(listHeight / 2), stageConfig.params.length - listHeight));
  const endIdx = Math.min(stageConfig.params.length, startIdx + listHeight);
  const visibleParams = stageConfig.params.slice(startIdx, endIdx);

  const shortcuts = [
    { key: '↑↓', label: 'Navigate' },
    { key: 'Space', label: 'Edit' },
    ...(stageName === 'synthesis' ? [{ key: 'r', label: 'Retrosynthesis' }] : []),
    { key: 'Enter', label: 'Save & Back' },
    { key: '←/Esc', label: 'Back' },
  ];

  if (loading) {
    return (
      <Box flexDirection="column" padding={1}>
        <Header title="Pipeline Wizard" subtitle={`Configure ${stageConfig.displayName}`} />
        <Spinner label="Loading configuration..." />
      </Box>
    );
  }

  const focusedParam = stageConfig.params[focusedIndex];
  const focusedDescription = focusedParam?.description || 'No description available.';

  return (
    <Box flexDirection="column" padding={1}>
      <Header
        title="Pipeline Wizard"
        subtitle={`Configure ${stageConfig.displayName}${isDirty ? ' *' : ''}`}
      />

      {/* Stage-specific hints */}
      {stageName === 'synthesis' && (
        <Box marginBottom={1}>
          <Text dimColor>Synthesis scoring (SA, SYBA, RA) and retrosynthesis route search. </Text>
          <Text color="cyan">[r]</Text>
          <Text dimColor> for AiZynthFinder model paths.</Text>
        </Box>
      )}
      {stageName === 'mol_prep' && (
        <Box marginBottom={1}>
          <Text dimColor>Datamol cleanup and strict molecule normalization before downstream stages.</Text>
        </Box>
      )}
      {stageName === 'descriptors' && (
        <Box marginBottom={1}>
          <Text dimColor>Molecular property calculation and filtering by descriptor ranges.</Text>
        </Box>
      )}
      {stageName === 'struct_filters' && (
        <Box marginBottom={1}>
          <Text dimColor>Structural alerts: PAINS, NIBR, Lilly medchem filters.</Text>
        </Box>
      )}
      {stageName === 'docking' && (
        <Box marginBottom={1}>
          <Text dimColor>Molecular docking with smina (CPU) or gnina (GPU+CNN scoring).</Text>
        </Box>
      )}

      <Box flexDirection={showSideDescription ? 'row' : 'column'} marginY={1} width={contentWidth}>
        <Box flexDirection="column" flexGrow={1} flexShrink={1} paddingRight={showSideDescription ? 2 : 0}>
          {visibleParams.map((param, index) => {
            const actualIndex = startIdx + index;
            const isFocused = focusedIndex === actualIndex;
            const isEditing = isFocused && editMode !== 'none';
            const value = getValue(param);

            return (
              <Box key={`${param.key}-${actualIndex}`} flexDirection="column">
                <Box>
                  <Text color={isFocused ? 'cyan' : 'gray'}>
                    {isFocused ? '▸ ' : '  '}
                  </Text>
                  <Text color={isFocused ? 'white' : 'gray'}>{param.label.padEnd(20)}</Text>

                  {isEditing ? (
                    <Box>
                      {param.type === 'range' ? (
                        <>
                          {editMode === 'min' ? (
                            <>
                              <TextInput value={editValue} onChange={setEditValue} focus={true} />
                              <Text dimColor> - {String((value as RangeValue)?.max ?? '')}</Text>
                            </>
                          ) : (
                            <>
                              <Text dimColor>{String((value as RangeValue)?.min ?? '')} - </Text>
                              <TextInput value={editValue} onChange={setEditValue} focus={true} />
                            </>
                          )}
                        </>
                      ) : (
                        <TextInput value={editValue} onChange={setEditValue} focus={true} />
                      )}
                    </Box>
                  ) : (
                    <Text color={getParamColor(param.type, value)}>
                      {formatValue(param)}
                    </Text>
                  )}
                </Box>
                {!showSideDescription && isFocused && param.description && (
                  <Box paddingLeft={4}>
                    <Text dimColor italic>{param.description}</Text>
                  </Box>
                )}
              </Box>
            );
          })}
        </Box>

        {showSideDescription && (
          <Box flexDirection="column" width={descriptionWidth}>
            <Text color="cyan" bold>Description</Text>
            <Text dimColor wrap="wrap">{focusedDescription}</Text>
          </Box>
        )}
      </Box>

      {stageConfig.params.length > listHeight && (
        <Box>
          <Text dimColor>({focusedIndex + 1}/{stageConfig.params.length})</Text>
        </Box>
      )}

      <Footer shortcuts={shortcuts} />
    </Box>
  );
}

// Individual screen wrappers
export function WizardConfigDescriptors(): React.ReactElement {
  return <QuickConfig stageName="descriptors" />;
}

export function WizardConfigMolPrep(): React.ReactElement {
  return <QuickConfig stageName="mol_prep" />;
}

export function WizardConfigFilters(): React.ReactElement {
  return <QuickConfig stageName="struct_filters" />;
}

export function WizardConfigSynthesis(): React.ReactElement {
  return <QuickConfig stageName="synthesis" />;
}

export function WizardConfigDocking(): React.ReactElement {
  return <QuickConfig stageName="docking" />;
}

export function WizardConfigDockingFilters(): React.ReactElement {
  return <QuickConfig stageName="docking_filters" />;
}

export default QuickConfig;
