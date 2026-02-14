import React, { useState, useEffect } from 'react';
import { Box, Text, useInput } from 'ink';
import { Header } from '../components/Header.js';
import { Footer } from '../components/Footer.js';
import { Spinner } from '../components/Spinner.js';
import { FileBrowser } from '../components/FileBrowser.js';
import { SearchIndicator } from '../components/SearchIndicator.js';
import { useStore } from '../store/index.js';
import { getBridge } from '../services/python-bridge.js';
import { useSearch } from '../hooks/useSearch.js';
import type { FiltersConfig } from '../types/index.js';

const ALL_RULESETS = [
  'Dundee', 'BMS', 'Inpharmatica', 'LD50-Oral', 'Glaxo', 'PAINS',
  'AlphaScreen-Hitters', 'Frequent-Hitter', 'Chelator', 'Toxicophore',
  'Skin', 'Reactive-Unstable-Toxic', 'SureChEMBL', 'Electrophilic',
  'Genotoxic-Carcinogenicity', 'MLSMR', 'LINT', 'GST-Hitters',
  'Non-Genotoxic-Carcinogenicity', 'HIS-Hitters', 'LuciferaseInhibitor',
];

interface SettingItem {
  key: string;
  label: string;
  type: 'boolean' | 'path' | 'select';
  options?: string[];
  description?: string;
}

const settingsItems: SettingItem[] = [
  { key: 'run', label: 'Run Stage', type: 'boolean', description: 'Enable/disable structural filters stage' },
  { key: 'filter_data', label: 'Filter Data', type: 'boolean', description: 'Apply filters to remove flagged molecules' },
  { key: 'alerts_data_path', label: 'Alerts Data Path', type: 'path', description: 'Path to structural alerts SMARTS database' },
  { key: 'calculate_common_alerts', label: 'Common Alerts', type: 'boolean', description: 'Check for PAINS, Dundee, Glaxo, etc. alerts' },
  { key: 'calculate_molgraph_stats', label: 'MolGraph Stats', type: 'boolean', description: 'Calculate molecular graph statistics' },
  { key: 'calculate_molcomplexity', label: 'Mol Complexity', type: 'boolean', description: 'Compute molecular complexity metrics' },
  { key: 'calculate_NIBR', label: 'NIBR Filters', type: 'boolean', description: 'Apply Novartis NIBR structural filters' },
  { key: 'nibr_scheduler', label: 'NIBR Scheduler', type: 'select', options: ['threads', 'processes'], description: 'Parallelization: threads (faster) or processes (safer)' },
  { key: 'calculate_bredt', label: 'Bredt Filters', type: 'boolean', description: 'Check for Bredt rule violations (strained rings)' },
  { key: 'calculate_lilly', label: 'Lilly Filters', type: 'boolean', description: 'Apply Eli Lilly medchem demerits rules' },
  { key: 'lilly_scheduler', label: 'Lilly Scheduler', type: 'select', options: ['threads', 'processes'], description: 'Parallelization: threads (faster) or processes (safer)' },
];

type ViewMode = 'settings' | 'rulesets';

type SettingType = 'boolean' | 'path' | 'select';

function getSettingColor(type: SettingType, value: unknown): 'green' | 'red' | 'cyan' | 'blue' {
  if (type === 'boolean') {
    return value ? 'green' : 'red';
  }
  if (type === 'select') {
    return 'cyan';
  }
  return 'blue';
}

export function ConfigFilters(): React.ReactElement {
  const setScreen = useStore((state) => state.setScreen);
  const config = useStore((state) => state.configs.filters);
  const setConfig = useStore((state) => state.setConfig);
  const isBackendReady = useStore((state) => state.isBackendReady);
  const showToast = useStore((state) => state.showToast);
  const showConfirm = useStore((state) => state.showConfirm);
  const setSearchActive = useStore((state) => state.setSearchActive);
  const setSearchQuery = useStore((state) => state.setSearchQuery);

  const [loading, setLoading] = useState(!config);
  const [viewMode, setViewMode] = useState<ViewMode>('settings');
  const [rawConfig, setRawConfig] = useState<FiltersConfig | null>(config);
  const [rulesets, setRulesets] = useState<Record<string, boolean>>({});
  const [selectedIndex, setSelectedIndex] = useState(0);
  const [error, setError] = useState<string | null>(null);
  const [scrollOffset, setScrollOffset] = useState(0);
  const [browsingField, setBrowsingField] = useState<string | null>(null);
  const [isDirty, setIsDirty] = useState(false);

  const visibleRows = 16;

  // Search for rulesets
  const {
    searchQuery,
    searchActive,
    filteredItems: filteredRulesets,
    handleSearchInput,
    highlightMatch,
  } = useSearch({
    items: ALL_RULESETS,
    searchFields: (ruleset) => [ruleset],
    enabled: viewMode === 'rulesets',
  });

  // Sync search state with global store
  useEffect(() => {
    if (viewMode === 'rulesets') {
      setSearchActive(searchActive);
      setSearchQuery(searchQuery);
    }
    return () => {
      setSearchActive(false);
      setSearchQuery('');
    };
  }, [searchActive, searchQuery, viewMode, setSearchActive, setSearchQuery]);

  // Reset selection when search changes
  useEffect(() => {
    setSelectedIndex(0);
    setScrollOffset(0);
  }, [searchQuery]);

  useEffect(() => {
    if (!config && isBackendReady) {
      loadConfig();
    } else if (config) {
      parseConfig(config);
      setRawConfig(config);
      setLoading(false);
    }
  }, [config, isBackendReady]);

  const parseConfig = (cfg: FiltersConfig) => {
    const rs: Record<string, boolean> = {};
    for (const ruleset of ALL_RULESETS) {
      rs[ruleset] = cfg.include_rulesets?.includes(ruleset) ?? false;
    }
    setRulesets(rs);
  };

  const loadConfig = async () => {
    try {
      const bridge = getBridge();
      const data = await bridge.loadConfig('filters') as unknown as FiltersConfig;
      setConfig('filters', data);
      setRawConfig(data);
      parseConfig(data);
    } catch (err) {
      setError(String(err));
    } finally {
      setLoading(false);
    }
  };

  const saveConfig = async () => {
    if (!rawConfig) return;

    try {
      const bridge = getBridge();

      // Build include_rulesets from rulesets state
      const includeRulesets = Object.entries(rulesets)
        .filter(([_, enabled]) => enabled)
        .map(([name]) => name);

      const newConfig: FiltersConfig = {
        ...rawConfig,
        include_rulesets: includeRulesets,
      };

      await bridge.saveConfig('filters', newConfig as unknown as Record<string, unknown>);
      setConfig('filters', newConfig);
      setRawConfig(newConfig);
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
        message: 'Discard changes to filters configuration?',
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

  const getSettingValue = (key: string): unknown => {
    if (!rawConfig) return undefined;
    return rawConfig[key as keyof FiltersConfig];
  };

  const setSettingValue = (key: string, value: unknown): void => {
    if (!rawConfig) return;
    setRawConfig({ ...rawConfig, [key]: value } as FiltersConfig);
    setIsDirty(true);
  };

  const toggleRuleset = (ruleset: string) => {
    setRulesets({ ...rulesets, [ruleset]: !rulesets[ruleset] });
    setIsDirty(true);
  };

  useInput((input, key) => {
    if (loading || browsingField) return;

    // Handle search input in rulesets mode
    if (viewMode === 'rulesets' && handleSearchInput(input, key)) {
      return;
    }

    const currentItems = viewMode === 'settings' ? settingsItems : filteredRulesets;
    const maxIndex = currentItems.length - 1;

    if (key.upArrow) {
      const newIndex = Math.max(0, selectedIndex - 1);
      setSelectedIndex(newIndex);
      if (newIndex < scrollOffset) {
        setScrollOffset(newIndex);
      }
    } else if (key.downArrow) {
      const newIndex = Math.min(maxIndex, selectedIndex + 1);
      setSelectedIndex(newIndex);
      if (newIndex >= scrollOffset + visibleRows) {
        setScrollOffset(newIndex - visibleRows + 1);
      }
    } else if (key.return || input === ' ' || input === 'e') {
      if (viewMode === 'settings') {
        const item = settingsItems[selectedIndex];
        if (item.type === 'boolean') {
          setSettingValue(item.key, !getSettingValue(item.key));
        } else if (item.type === 'select' && item.options) {
          const current = getSettingValue(item.key) as string;
          const idx = item.options.indexOf(current);
          const nextIdx = (idx + 1) % item.options.length;
          setSettingValue(item.key, item.options[nextIdx]);
        } else if (item.type === 'path') {
          setBrowsingField(item.key);
        }
      } else {
        // Rulesets view - use filtered list
        const ruleset = filteredRulesets[selectedIndex];
        if (ruleset) {
          toggleRuleset(ruleset);
        }
      }
    } else if (input === 'r') {
      setViewMode(viewMode === 'settings' ? 'rulesets' : 'settings');
      setSelectedIndex(0);
      setScrollOffset(0);
    } else if (input === 'a' && viewMode === 'rulesets' && !searchActive) {
      // Select all rulesets (only visible/filtered ones when searching)
      const toEnable = searchQuery ? filteredRulesets : ALL_RULESETS;
      const updated = { ...rulesets };
      for (const rs of toEnable) {
        updated[rs] = true;
      }
      setRulesets(updated);
      setIsDirty(true);
    } else if (input === 'n' && viewMode === 'rulesets' && !searchActive) {
      // Deselect all rulesets
      const toDisable = searchQuery ? filteredRulesets : ALL_RULESETS;
      const updated = { ...rulesets };
      for (const rs of toDisable) {
        updated[rs] = false;
      }
      setRulesets(updated);
      setIsDirty(true);
    } else if (input === 's') {
      saveConfig();
    } else if (key.escape || key.leftArrow || input === 'q') {
      if (viewMode === 'rulesets') {
        setViewMode('settings');
        setSelectedIndex(0);
        setScrollOffset(0);
      } else {
        handleExit();
      }
    }
  });

  const handlePathSelect = (path: string) => {
    if (browsingField) {
      setSettingValue(browsingField, path);
    }
    setBrowsingField(null);
  };

  const handleBrowseCancel = () => {
    setBrowsingField(null);
  };

  const settingsShortcuts = [
    { key: '↑↓', label: 'Navigate' },
    { key: 'Space', label: 'Toggle' },
    { key: 'r', label: 'Rulesets' },
    { key: 's', label: 'Save' },
    { key: '←/Esc', label: 'Back' },
  ];

  const rulesetsShortcuts = [
    { key: '↑↓', label: 'Navigate' },
    { key: 'Space', label: 'Toggle' },
    { key: 'a/n', label: 'All/None' },
    { key: 's', label: 'Save' },
    { key: '←/Esc', label: 'Settings' },
  ];

  if (loading) {
    return (
      <Box flexDirection="column" padding={1}>
        <Header title="Struct Filters Config" />
        <Spinner label="Loading configuration..." />
      </Box>
    );
  }

  if (browsingField) {
    const currentValue = String(getSettingValue(browsingField) || process.cwd());
    return (
      <Box flexDirection="column" padding={1}>
        <Header title="Select File" subtitle={browsingField} />
        <FileBrowser
          initialPath={currentValue}
          onSelect={handlePathSelect}
          onCancel={handleBrowseCancel}
        />
      </Box>
    );
  }

  if (viewMode === 'rulesets') {
    const displayRulesets = filteredRulesets;
    const visibleRulesets = displayRulesets.slice(scrollOffset, scrollOffset + visibleRows);
    const enabledCount = Object.values(rulesets).filter(Boolean).length;

    return (
      <Box flexDirection="column" padding={1}>
        <Header title="Struct Filters Config" subtitle={isDirty ? 'Rulesets *' : 'Rulesets'} />

        <SearchIndicator active={searchActive} query={searchQuery} />

        {error && (
          <Box marginY={1}>
            <Text color="red">Error: {error}</Text>
          </Box>
        )}

        <Box marginY={1}>
          <Text color="cyan" bold>─ Rulesets ({enabledCount}/{ALL_RULESETS.length} enabled) ─</Text>
        </Box>

        {displayRulesets.length === 0 ? (
          <Box marginY={1}>
            <Text dimColor>No rulesets match "{searchQuery}"</Text>
          </Box>
        ) : (
          <Box flexDirection="column">
            {visibleRulesets.map((ruleset, index) => {
              const actualIndex = scrollOffset + index;
              const isSelected = actualIndex === selectedIndex;
              const enabled = rulesets[ruleset];
              const nameParts = searchQuery ? highlightMatch(ruleset) : [{ text: ruleset, highlighted: false }];

              return (
                <Box key={ruleset}>
                  <Text color={isSelected ? 'cyan' : 'white'}>
                    {isSelected ? '▶ ' : '  '}
                  </Text>
                  <Text color={enabled ? 'green' : 'red'}>
                    {enabled ? '☑' : '☐'}
                  </Text>
                  <Text> </Text>
                  {nameParts.map((part, pi) => (
                    <Text
                      key={pi}
                      color={isSelected ? (part.highlighted ? 'yellow' : 'white') : (part.highlighted ? 'yellow' : 'gray')}
                      bold={part.highlighted}
                    >
                      {part.text}
                    </Text>
                  ))}
                </Box>
              );
            })}
          </Box>
        )}

        <Box marginTop={1}>
          <Text dimColor>
            {displayRulesets.length > 0
              ? `${scrollOffset + 1}-${Math.min(scrollOffset + visibleRows, displayRulesets.length)} of ${displayRulesets.length}`
              : '0'}
            {searchQuery && ` (filtered from ${ALL_RULESETS.length})`}
          </Text>
        </Box>

        <Footer shortcuts={rulesetsShortcuts} />
      </Box>
    );
  }

  // Settings view
  const visibleSettings = settingsItems.slice(scrollOffset, scrollOffset + visibleRows);

  return (
    <Box flexDirection="column" padding={1}>
      <Header title="Struct Filters Config" subtitle={isDirty ? 'config_structFilters.yml *' : 'config_structFilters.yml'} />

      {error && (
        <Box marginY={1}>
          <Text color="red">Error: {error}</Text>
        </Box>
      )}

      <Box marginY={1}>
        <Text color="cyan" bold>─ Settings ─</Text>
      </Box>

      <Box flexDirection="column">
        {visibleSettings.map((item, index) => {
          const actualIndex = scrollOffset + index;
          const isSelected = actualIndex === selectedIndex;
          const value = getSettingValue(item.key);

          return (
            <Box key={item.key} flexDirection="column">
              <Box>
                <Text color={isSelected ? 'cyan' : 'white'}>
                  {isSelected ? '▶ ' : '  '}
                </Text>
                <Text dimColor>{item.label.padEnd(22)}</Text>
                <Text color={getSettingColor(item.type, value)}>
                  {item.type === 'boolean'
                    ? (value ? 'Yes' : 'No')
                    : String(value || '(not set)')}
                </Text>
              </Box>
              {isSelected && item.description && (
                <Box paddingLeft={4}>
                  <Text dimColor italic>{item.description}</Text>
                </Box>
              )}
            </Box>
          );
        })}
      </Box>

      <Box marginTop={1}>
        <Text dimColor>Press 'r' to edit rulesets ({Object.values(rulesets).filter(Boolean).length} enabled)</Text>
      </Box>

      <Footer shortcuts={settingsShortcuts} />
    </Box>
  );
}

export default ConfigFilters;
