import React, { useState } from 'react';
import { Box, Text, useInput } from 'ink';
import TextInput from 'ink-text-input';
import { useInputLock } from '../hooks/useInputLock.js';
import { useTheme } from '../theme/context.js';

export interface DescriptorRow {
  name: string;
  displayName: string;
  min: number | string;
  max: number | string;
}

interface DescriptorTableProps {
  descriptors: DescriptorRow[];
  onChange: (descriptors: DescriptorRow[]) => void;
  onBack: () => void;
}

export function DescriptorTable({ 
  descriptors, 
  onChange, 
  onBack,
}: DescriptorTableProps): React.ReactElement {
  const { theme } = useTheme();
  const [selectedRow, setSelectedRow] = useState(0);
  const [selectedCol, setSelectedCol] = useState<'min' | 'max' | null>(null);
  const [editValue, setEditValue] = useState('');
  const [scrollOffset, setScrollOffset] = useState(0);
  useInputLock(selectedCol !== null);

  const visibleRows = 15;
  const visibleDescriptors = descriptors.slice(scrollOffset, scrollOffset + visibleRows);

  useInput((input, key) => {
    if (selectedCol !== null) {
      // Edit mode
      if (key.escape) {
        setSelectedCol(null);
      } else if (key.return) {
        const newDescriptors = [...descriptors];
        const descriptor = newDescriptors[selectedRow];
        const value = editValue === '' ? 0 : (isNaN(parseFloat(editValue)) ? editValue : parseFloat(editValue));
        if (selectedCol === 'min') {
          descriptor.min = value;
        } else {
          descriptor.max = value;
        }
        onChange(newDescriptors);
        setSelectedCol(null);
      }
      return;
    }

    if (key.upArrow) {
      const newRow = Math.max(0, selectedRow - 1);
      setSelectedRow(newRow);
      if (newRow < scrollOffset) {
        setScrollOffset(newRow);
      }
    } else if (key.downArrow) {
      const newRow = Math.min(descriptors.length - 1, selectedRow + 1);
      setSelectedRow(newRow);
      if (newRow >= scrollOffset + visibleRows) {
        setScrollOffset(newRow - visibleRows + 1);
      }
    } else if (key.leftArrow || key.rightArrow) {
      // Could implement column navigation here
    } else if (input === 'm') {
      setEditValue(String(descriptors[selectedRow].min));
      setSelectedCol('min');
    } else if (input === 'x') {
      setEditValue(String(descriptors[selectedRow].max));
      setSelectedCol('max');
    } else if (key.escape) {
      onBack();
    }
  });

  return (
    <Box flexDirection="column">
      {/* Header */}
      <Box marginBottom={1}>
        <Text color={theme.palette.textMuted}>  </Text>
        <Text bold color={theme.palette.primary}>{'Descriptor'.padEnd(25)}</Text>
        <Text bold color={theme.palette.primary}>{'Min'.padStart(12)}</Text>
        <Text bold color={theme.palette.primary}>{'Max'.padStart(12)}</Text>
      </Box>
      
      {/* Rows */}
      {visibleDescriptors.map((desc, index) => {
        const actualIndex = scrollOffset + index;
        const isSelected = actualIndex === selectedRow;
        
        return (
          <Box key={desc.name}>
            <Text color={isSelected ? theme.palette.primary : theme.palette.text}>
              {isSelected ? '▶ ' : '  '}
            </Text>
            <Text color={isSelected ? theme.palette.text : theme.palette.textMuted}>
              {desc.displayName.padEnd(25)}
            </Text>
            {selectedCol === 'min' && isSelected ? (
              <Box width={12}>
                <TextInput
                  value={editValue}
                  onChange={setEditValue}
                  focus={true}
                />
              </Box>
            ) : (
              <Text color={theme.palette.accent}>{String(desc.min).padStart(12)}</Text>
            )}
            {selectedCol === 'max' && isSelected ? (
              <Box width={12}>
                <TextInput
                  value={editValue}
                  onChange={setEditValue}
                  focus={true}
                />
              </Box>
            ) : (
              <Text color={theme.palette.accent}>{String(desc.max).padStart(12)}</Text>
            )}
          </Box>
        );
      })}
      
      {/* Scroll indicator */}
      {descriptors.length > visibleRows && (
        <Box marginTop={1}>
          <Text color={theme.palette.textMuted}>
            Showing {scrollOffset + 1}-{Math.min(scrollOffset + visibleRows, descriptors.length)} of {descriptors.length}
          </Text>
        </Box>
      )}
      
      <Box marginTop={1} gap={2}>
        <Text color={theme.palette.text}><Text color={theme.palette.primary}>[m]</Text> Edit min</Text>
        <Text color={theme.palette.text}><Text color={theme.palette.primary}>[x]</Text> Edit max</Text>
        <Text color={theme.palette.text}><Text color={theme.palette.primary}>[Esc]</Text> Back</Text>
      </Box>
    </Box>
  );
}

export default DescriptorTable;
