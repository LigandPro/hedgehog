import React, { useState, useEffect } from 'react';
import { Box, Text, useInput } from 'ink';
import { Header } from '../components/Header.js';
import { Footer } from '../components/Footer.js';
import { Spinner } from '../components/Spinner.js';
import { useStore } from '../store/index.js';
import { getBridge } from '../services/python-bridge.js';
import type { FiltersConfig } from '../types/index.js';

const ALL_RULESETS = [
  'Dundee', 'BMS', 'Inpharmatica', 'LD50-Oral', 'Glaxo', 'PAINS',
  'AlphaScreen-Hitters', 'Frequent-Hitter', 'Chelator', 'Toxicophore',
  'Skin', 'Reactive-Unstable-Toxic', 'SureChEMBL', 'Electrophilic',
  'Genotoxic-Carcinogenicity', 'MLSMR', 'LINT', 'GST-Hitters',
  'Non-Genotoxic-Carcinogenicity', 'HIS-Hitters', 'LuciferaseInhibitor',
];

interface ToggleItem {
  key: string;
  label: string;
  enabled: boolean;
  isRuleset?: boolean;
}

export function ConfigFilters(): React.ReactElement {
  const setScreen = useStore((state) => state.setScreen);
  const config = useStore((state) => state.configs.filters);
  const setConfig = useStore((state) => state.setConfig);
  const isBackendReady = useStore((state) => state.isBackendReady);
  
  const [loading, setLoading] = useState(!config);
  const [items, setItems] = useState<ToggleItem[]>([]);
  const [selectedIndex, setSelectedIndex] = useState(0);
  const [error, setError] = useState<string | null>(null);
  const [saved, setSaved] = useState(false);
  const [scrollOffset, setScrollOffset] = useState(0);

  const visibleRows = 18;

  useEffect(() => {
    if (!config && isBackendReady) {
      loadConfig();
    } else if (config) {
      parseConfig(config);
      setLoading(false);
    }
  }, [config, isBackendReady]);

  const parseConfig = (cfg: FiltersConfig) => {
    const toggleItems: ToggleItem[] = [
      { key: 'run', label: 'Run Stage', enabled: cfg.run },
      { key: 'filter_data', label: 'Filter Data', enabled: cfg.filter_data },
      { key: 'calculate_common_alerts', label: 'Common Alerts', enabled: cfg.calculate_common_alerts },
      { key: 'calculate_molgraph_stats', label: 'MolGraph Stats', enabled: cfg.calculate_molgraph_stats },
      { key: 'calculate_molcomplexity', label: 'Mol Complexity', enabled: cfg.calculate_molcomplexity },
      { key: 'calculate_NIBR', label: 'NIBR Filters', enabled: cfg.calculate_NIBR },
      { key: 'calculate_bredt', label: 'Bredt Filters', enabled: cfg.calculate_bredt },
      { key: 'calculate_lilly', label: 'Lilly Filters', enabled: cfg.calculate_lilly },
    ];

    // Add rulesets
    for (const ruleset of ALL_RULESETS) {
      toggleItems.push({
        key: `ruleset_${ruleset}`,
        label: ruleset,
        enabled: cfg.include_rulesets?.includes(ruleset) ?? false,
        isRuleset: true,
      });
    }

    setItems(toggleItems);
  };

  const loadConfig = async () => {
    try {
      const bridge = getBridge();
      const data = await bridge.loadConfig('filters') as unknown as FiltersConfig;
      setConfig('filters', data);
      parseConfig(data);
    } catch (err) {
      setError(String(err));
    } finally {
      setLoading(false);
    }
  };

  const saveConfig = async () => {
    if (!config) return;
    
    try {
      const bridge = getBridge();
      const newConfig: Partial<FiltersConfig> = { ...config };
      
      // Extract settings
      for (const item of items) {
        if (item.isRuleset) continue;
        (newConfig as any)[item.key] = item.enabled;
      }
      
      // Extract rulesets
      newConfig.include_rulesets = items
        .filter(item => item.isRuleset && item.enabled)
        .map(item => item.key.replace('ruleset_', ''));
      
      await bridge.saveConfig('filters', newConfig as unknown as Record<string, unknown>);
      setConfig('filters', newConfig as unknown as FiltersConfig);
      setSaved(true);
      setTimeout(() => setSaved(false), 2000);
    } catch (err) {
      setError(String(err));
    }
  };

  useInput((input, key) => {
    if (loading) return;

    if (key.upArrow) {
      const newIndex = Math.max(0, selectedIndex - 1);
      setSelectedIndex(newIndex);
      if (newIndex < scrollOffset) {
        setScrollOffset(newIndex);
      }
    } else if (key.downArrow) {
      const newIndex = Math.min(items.length - 1, selectedIndex + 1);
      setSelectedIndex(newIndex);
      if (newIndex >= scrollOffset + visibleRows) {
        setScrollOffset(newIndex - visibleRows + 1);
      }
    } else if (key.return || input === ' ') {
      const newItems = [...items];
      newItems[selectedIndex].enabled = !newItems[selectedIndex].enabled;
      setItems(newItems);
    } else if (input === 's') {
      saveConfig();
    } else if (key.escape || input === 'q') {
      setScreen('welcome');
    }
  });

  const shortcuts = [
    { key: '↑↓', label: 'Navigate' },
    { key: 'Space', label: 'Toggle' },
    { key: 's', label: 'Save' },
    { key: 'Esc', label: 'Back' },
  ];

  if (loading) {
    return (
      <Box flexDirection="column" padding={1}>
        <Header title="Struct Filters Config" />
        <Spinner label="Loading configuration..." />
      </Box>
    );
  }

  const visibleItems = items.slice(scrollOffset, scrollOffset + visibleRows);

  return (
    <Box flexDirection="column" padding={1}>
      <Header title="Struct Filters Config" subtitle="config_structFilters.yml" />
      
      {error && (
        <Box marginY={1}>
          <Text color="red">Error: {error}</Text>
        </Box>
      )}
      
      {saved && (
        <Box marginY={1}>
          <Text color="green">✓ Configuration saved</Text>
        </Box>
      )}
      
      <Box flexDirection="column" marginY={1}>
        {visibleItems.map((item, index) => {
          const actualIndex = scrollOffset + index;
          const isSelected = actualIndex === selectedIndex;
          const isHeader = actualIndex === 8; // First ruleset
          
          return (
            <React.Fragment key={item.key}>
              {isHeader && actualIndex === scrollOffset && (
                <Box marginTop={1} marginBottom={0}>
                  <Text color="cyan" bold>Rulesets:</Text>
                </Box>
              )}
              <Box>
                <Text color={isSelected ? 'cyan' : 'white'}>
                  {isSelected ? '▶ ' : '  '}
                </Text>
                <Text color={item.enabled ? 'green' : 'red'}>
                  {item.enabled ? '☑' : '☐'}
                </Text>
                <Text color={isSelected ? 'white' : 'gray'}> {item.label}</Text>
              </Box>
            </React.Fragment>
          );
        })}
      </Box>
      
      {items.length > visibleRows && (
        <Box>
          <Text dimColor>
            Showing {scrollOffset + 1}-{Math.min(scrollOffset + visibleRows, items.length)} of {items.length}
          </Text>
        </Box>
      )}
      
      <Footer shortcuts={shortcuts} />
    </Box>
  );
}

export default ConfigFilters;
