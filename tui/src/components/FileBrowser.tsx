import React, { useState, useEffect } from 'react';
import { Box, Text, useInput } from 'ink';
import TextInput from 'ink-text-input';
import { getBridge } from '../services/python-bridge.js';

interface FileBrowserProps {
  initialPath: string;
  extensions?: string[];
  onSelect: (path: string) => void;
  onCancel: () => void;
}

interface FileEntry {
  name: string;
  path: string;
  isDirectory: boolean;
}

export function FileBrowser({ 
  initialPath, 
  extensions,
  onSelect, 
  onCancel,
}: FileBrowserProps): React.ReactElement {
  const [currentPath, setCurrentPath] = useState(initialPath);
  const [files, setFiles] = useState<FileEntry[]>([]);
  const [selectedIndex, setSelectedIndex] = useState(0);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [manualInput, setManualInput] = useState(false);
  const [inputValue, setInputValue] = useState(initialPath);

  useEffect(() => {
    loadFiles();
  }, [currentPath]);

  const loadFiles = async () => {
    setLoading(true);
    setError(null);
    try {
      const bridge = getBridge();
      const result = await bridge.call<FileEntry[]>('list_directory', { 
        path: currentPath,
        extensions,
      });
      setFiles([{ name: '..', path: currentPath + '/..', isDirectory: true }, ...result]);
      setSelectedIndex(0);
    } catch (err) {
      setError(String(err));
    } finally {
      setLoading(false);
    }
  };

  useInput((input, key) => {
    if (manualInput) {
      if (key.escape) {
        setManualInput(false);
      } else if (key.return) {
        onSelect(inputValue);
      }
      return;
    }

    if (key.upArrow) {
      setSelectedIndex(Math.max(0, selectedIndex - 1));
    } else if (key.downArrow) {
      setSelectedIndex(Math.min(files.length - 1, selectedIndex + 1));
    } else if (key.return) {
      const selected = files[selectedIndex];
      if (selected) {
        if (selected.isDirectory) {
          setCurrentPath(selected.path);
        } else {
          onSelect(selected.path);
        }
      }
    } else if (input === 'p') {
      setInputValue(currentPath);
      setManualInput(true);
    } else if (key.escape || input === 'q') {
      onCancel();
    }
  });

  if (loading) {
    return <Text>Loading...</Text>;
  }

  if (manualInput) {
    return (
      <Box flexDirection="column">
        <Text>Enter path:</Text>
        <TextInput
          value={inputValue}
          onChange={setInputValue}
          focus={true}
        />
        <Text dimColor>[Enter] Confirm  [Esc] Cancel</Text>
      </Box>
    );
  }

  return (
    <Box flexDirection="column">
      <Box marginBottom={1}>
        <Text dimColor>Path: </Text>
        <Text color="cyan">{currentPath}</Text>
      </Box>
      
      {error && (
        <Text color="red">Error: {error}</Text>
      )}
      
      <Box flexDirection="column" height={15} overflowY="hidden">
        {files.slice(Math.max(0, selectedIndex - 7), selectedIndex + 8).map((file, index) => {
          const actualIndex = Math.max(0, selectedIndex - 7) + index;
          const isSelected = actualIndex === selectedIndex;
          return (
            <Box key={file.path}>
              <Text color={isSelected ? 'cyan' : 'white'}>
                {isSelected ? '▶ ' : '  '}
              </Text>
              <Text color={file.isDirectory ? 'blue' : 'white'}>
                {file.isDirectory ? '📁 ' : '📄 '}
                {file.name}
              </Text>
            </Box>
          );
        })}
      </Box>
      
      <Box marginTop={1} gap={2}>
        <Text><Text color="cyan">[Enter]</Text> Select</Text>
        <Text><Text color="cyan">[p]</Text> Manual path</Text>
        <Text><Text color="cyan">[Esc]</Text> Cancel</Text>
      </Box>
    </Box>
  );
}

export default FileBrowser;
