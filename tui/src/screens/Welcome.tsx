import React, { useState } from 'react';
import { Box, Text, useInput, useApp } from 'ink';
import SelectInput from 'ink-select-input';
import { Header } from '../components/Header.js';
import { Footer } from '../components/Footer.js';
import { useStore } from '../store/index.js';
import type { Screen } from '../types/index.js';

interface MenuItem {
  label: string;
  value: Screen | 'exit';
}

const menuItems: MenuItem[] = [
  { label: 'New Pipeline Run', value: 'pipelineRunner' },
  { label: 'Edit Main Config', value: 'configMain' },
  { label: 'Edit Descriptors', value: 'configDescriptors' },
  { label: 'Edit Struct Filters', value: 'configFilters' },
  { label: 'Edit Synthesis', value: 'configSynthesis' },
  { label: 'Edit Docking', value: 'configDocking' },
  { label: 'Exit', value: 'exit' },
];

export function Welcome(): React.ReactElement {
  const { exit } = useApp();
  const setScreen = useStore((state) => state.setScreen);
  const isBackendReady = useStore((state) => state.isBackendReady);
  const [selectedIndex, setSelectedIndex] = useState(0);

  useInput((input, key) => {
    if (input === 'q' || (key.escape && selectedIndex === menuItems.length - 1)) {
      exit();
    }
  });

  const handleSelect = (item: MenuItem) => {
    if (item.value === 'exit') {
      exit();
    } else {
      setScreen(item.value);
    }
  };

  const shortcuts = [
    { key: '↑↓', label: 'Navigate' },
    { key: 'Enter', label: 'Select' },
    { key: 'q', label: 'Quit' },
  ];

  return (
    <Box flexDirection="column" padding={1}>
      <Header 
        showLogo={true}
        subtitle="Molecular Design Pipeline"
      />
      
      {!isBackendReady && (
        <Box marginY={1}>
          <Text color="yellow">⚠ Backend not connected. Some features may not work.</Text>
        </Box>
      )}
      
      <Box flexDirection="column" marginY={1}>
        {menuItems.map((item, index) => {
          const isSelected = index === selectedIndex;
          return (
            <Box key={item.value}>
              <Text color={isSelected ? 'cyan' : 'white'}>
                {isSelected ? '▶ ' : '  '}
                {item.label}
              </Text>
            </Box>
          );
        })}
      </Box>
      
      {/* Hidden SelectInput for keyboard handling */}
      <Box display="none">
        <SelectInput
          items={menuItems}
          onSelect={handleSelect}
          onHighlight={(item) => setSelectedIndex(menuItems.indexOf(item))}
        />
      </Box>
      
      <Footer shortcuts={shortcuts} />
    </Box>
  );
}

export default Welcome;
