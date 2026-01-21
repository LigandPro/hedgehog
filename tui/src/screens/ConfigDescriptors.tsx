import React, { useState, useEffect } from 'react';
import { Box, Text, useInput } from 'ink';
import { Header } from '../components/Header.js';
import { Footer } from '../components/Footer.js';
import { Spinner } from '../components/Spinner.js';
import { DescriptorTable, type DescriptorRow } from '../components/DescriptorTable.js';
import { useStore } from '../store/index.js';
import { getBridge } from '../services/python-bridge.js';
import type { DescriptorsConfig } from '../types/index.js';


// Descriptor presets
const PRESETS: Record<string, Record<string, [number, number]>> = {
  'Drug-like': {
    molWt: [200, 500],
    logP: [-0.4, 5.6],
    hbd: [0, 5],
    hba: [0, 10],
    tpsa: [0, 140],
    n_rot_bonds: [0, 10],
  },
  'Lead-like': {
    molWt: [200, 350],
    logP: [-1, 3],
    hbd: [0, 3],
    hba: [0, 6],
    n_rot_bonds: [0, 7],
  },
  'Fragment-like': {
    molWt: [0, 300],
    logP: [-3, 3],
    hbd: [0, 3],
    hba: [0, 3],
    n_heavy_atoms: [0, 17],
  },
};

// List of descriptors with display names
const DESCRIPTOR_NAMES: Record<string, string> = {
  'n_atoms': 'Number of Atoms',
  'n_heavy_atoms': 'Heavy Atoms',
  'n_het_atoms': 'Heteroatoms',
  'n_N_atoms': 'Nitrogen Atoms',
  'fN_atoms': 'Fraction N Atoms',
  'molWt': 'Molecular Weight',
  'logP': 'logP',
  'clogP': 'clogP',
  'sw': 'Sw',
  'ring_size': 'Ring Size',
  'n_rings': 'Number of Rings',
  'n_aroma_rings': 'Aromatic Rings',
  'n_fused_aromatic_rings': 'Fused Aromatic Rings',
  'n_rigid_bonds': 'Rigid Bonds',
  'n_rot_bonds': 'Rotatable Bonds',
  'hbd': 'H-Bond Donors',
  'hba': 'H-Bond Acceptors',
  'fsp3': 'Fraction SP3',
  'tpsa': 'TPSA',
  'qed': 'QED',
};

export function ConfigDescriptors(): React.ReactElement {
  const setScreen = useStore((state) => state.setScreen);
  const config = useStore((state) => state.configs.descriptors);
  const setConfig = useStore((state) => state.setConfig);
  const isBackendReady = useStore((state) => state.isBackendReady);
  
  const [loading, setLoading] = useState(!config);
  const [descriptors, setDescriptors] = useState<DescriptorRow[]>([]);
  const [rawConfig, setRawConfig] = useState<Record<string, unknown> | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [saved, setSaved] = useState(false);
  const [presetApplied, setPresetApplied] = useState<string | null>(null);

  useEffect(() => {
    if (!config && isBackendReady) {
      loadConfig();
    } else if (config) {
      parseConfig(config as unknown as Record<string, unknown>);
      setRawConfig(config as unknown as Record<string, unknown>);
      setLoading(false);
    }
  }, [config, isBackendReady]);

  const parseConfig = (cfg: Record<string, unknown>) => {
    const rows: DescriptorRow[] = [];
    const borders = (cfg.borders || {}) as Record<string, unknown>;
    
    for (const [key, displayName] of Object.entries(DESCRIPTOR_NAMES)) {
      const minKey = `${key}_min`;
      const maxKey = `${key}_max`;
      
      rows.push({
        name: key,
        displayName,
        min: (borders[minKey] as number | string) ?? 0,
        max: (borders[maxKey] as number | string) ?? 100,
      });
    }
    
    setDescriptors(rows);
  };

  const loadConfig = async () => {
    try {
      const bridge = getBridge();
      const data = await bridge.loadConfig('descriptors');
      setConfig('descriptors', data as unknown as DescriptorsConfig);
      setRawConfig(data);
      parseConfig(data);
    } catch (err) {
      setError(String(err));
    } finally {
      setLoading(false);
    }
  };

  const applyPreset = (presetName: string) => {
    const preset = PRESETS[presetName];
    if (!preset) return;

    const newDescriptors = descriptors.map(desc => {
      const presetValues = preset[desc.name];
      if (presetValues) {
        return { ...desc, min: presetValues[0], max: presetValues[1] };
      }
      return desc;
    });

    setDescriptors(newDescriptors);
    setPresetApplied(presetName);
    setTimeout(() => setPresetApplied(null), 2000);
  };

  const saveConfig = async () => {
    if (!rawConfig) return;
    
    try {
      const bridge = getBridge();
      const newBorders: Record<string, number | string> = { ...(rawConfig.borders as Record<string, unknown>) } as Record<string, number | string>;
      
      for (const desc of descriptors) {
        newBorders[`${desc.name}_min`] = desc.min;
        newBorders[`${desc.name}_max`] = desc.max;
      }
      
      const newConfig = { ...rawConfig, borders: newBorders };
      await bridge.saveConfig('descriptors', newConfig);
      setConfig('descriptors', newConfig as unknown as DescriptorsConfig);
      setRawConfig(newConfig);
      setSaved(true);
      setTimeout(() => setSaved(false), 2000);
    } catch (err) {
      setError(String(err));
    }
  };

  const handleDescriptorChange = (newDescriptors: DescriptorRow[]) => {
    setDescriptors(newDescriptors);
  };

  const handleBack = () => {
    setScreen('welcome');
  };

  useInput((input) => {
    if (input === 's') {
      saveConfig();
    } else if (input === '1') {
      applyPreset('Drug-like');
    } else if (input === '2') {
      applyPreset('Lead-like');
    } else if (input === '3') {
      applyPreset('Fragment-like');
    }
  });

  const shortcuts = [
    { key: '↑↓', label: 'Navigate' },
    { key: 'm/x', label: 'Edit min/max' },
    { key: '1/2/3', label: 'Presets' },
    { key: 's', label: 'Save' },
    { key: 'Esc', label: 'Back' },
  ];

  if (loading) {
    return (
      <Box flexDirection="column" padding={1}>
        <Header title="Descriptors Config" />
        <Spinner label="Loading configuration..." />
      </Box>
    );
  }

  return (
    <Box flexDirection="column" padding={1}>
      <Header title="Descriptors Config" subtitle="config_descriptors.yml" />
      
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
      
      <DescriptorTable
        descriptors={descriptors}
        onChange={handleDescriptorChange}
        onBack={handleBack}
      />
      
      <Footer shortcuts={shortcuts} />
    </Box>
  );
}

export default ConfigDescriptors;
