import React, { useState } from 'react';
import { Box, Text, useInput } from 'ink';
import TextInput from 'ink-text-input';

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
  const [selectedRow, setSelectedRow] = useState(0);
  const [selectedCol, setSelectedCol] = useState<'min' | 'max' | null>(null);
  const [editValue, setEditValue] = useState('');
  const [scrollOffset, setScrollOffset] = useState(0);

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
    } else if (key.escape || input === 'q') {
      onBack();
    }
  });

  return (
    <Box flexDirection="column">
      {/* Header */}
      <Box marginBottom={1}>
        <Text color="gray">  </Text>
        <Text bold color="cyan">{'Descriptor'.padEnd(25)}</Text>
        <Text bold color="cyan">{'Min'.padStart(12)}</Text>
        <Text bold color="cyan">{'Max'.padStart(12)}</Text>
      </Box>
      
      {/* Rows */}
      {visibleDescriptors.map((desc, index) => {
        const actualIndex = scrollOffset + index;
        const isSelected = actualIndex === selectedRow;
        
        return (
          <Box key={desc.name}>
            <Text color={isSelected ? 'cyan' : 'white'}>
              {isSelected ? '▶ ' : '  '}
            </Text>
            <Text color={isSelected ? 'white' : 'gray'}>
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
              <Text color="yellow">{String(desc.min).padStart(12)}</Text>
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
              <Text color="yellow">{String(desc.max).padStart(12)}</Text>
            )}
          </Box>
        );
      })}
      
      {/* Scroll indicator */}
      {descriptors.length > visibleRows && (
        <Box marginTop={1}>
          <Text dimColor>
            Showing {scrollOffset + 1}-{Math.min(scrollOffset + visibleRows, descriptors.length)} of {descriptors.length}
          </Text>
        </Box>
      )}
      
      <Box marginTop={1} gap={2}>
        <Text><Text color="cyan">[m]</Text> Edit min</Text>
        <Text><Text color="cyan">[x]</Text> Edit max</Text>
        <Text><Text color="cyan">[Esc]</Text> Back</Text>
      </Box>
    </Box>
  );
}

export default DescriptorTable;
