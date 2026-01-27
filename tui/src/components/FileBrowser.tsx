import React, { useState, useEffect, useMemo } from 'react';
import { Box, Text, useInput } from 'ink';
import TextInput from 'ink-text-input';
import path from 'path';
import { getBridge } from '../services/python-bridge.js';
import { useTerminalSize } from '../hooks/useTerminalSize.js';

interface FileBrowserProps {
  initialPath: string;
  extensions?: string[];
  selectDirectory?: boolean;
  onSelect: (path: string) => void;
  onCancel: () => void;
  title?: string;
}

interface FileEntry {
  name: string;
  path: string;
  isDirectory: boolean;
}

/**
 * Get a valid starting directory from a path.
 * If path looks like a file (has extension), return parent directory.
 * If path is empty or invalid, return cwd.
 */
function getStartingDirectory(inputPath: string): string {
  if (!inputPath || inputPath.trim() === '') {
    return process.cwd();
  }

  // Check if path looks like a file (has extension in the last segment)
  const basename = path.basename(inputPath);
  const hasExtension = basename.includes('.') && !basename.startsWith('.');

  if (hasExtension) {
    // It's likely a file, use parent directory
    return path.dirname(inputPath);
  }

  return inputPath;
}

/**
 * Truncate path for display, keeping the end visible.
 */
function truncatePath(p: string, maxLen: number): string {
  if (!p || p.length <= maxLen) return p;
  return '...' + p.slice(-(maxLen - 3));
}

/**
 * Get text color for file entry based on selection state and type.
 */
function getFileTextColor(isSelected: boolean, isDirectory: boolean, matchesExtension: boolean): string {
  if (isSelected) return 'white';
  if (isDirectory) return 'blue';
  if (matchesExtension) return 'green';
  return 'gray';
}

/**
 * Get icon color for file entry based on selection state and type.
 */
function getFileIconColor(isSelected: boolean, isDirectory: boolean, matchesExtension: boolean): string {
  if (isSelected) return 'cyan';
  if (isDirectory) return 'blue';
  if (matchesExtension) return 'green';
  return 'gray';
}

export function FileBrowser({
  initialPath,
  extensions,
  selectDirectory = false,
  onSelect,
  onCancel,
  title,
}: FileBrowserProps): React.ReactElement {
  const startDir = getStartingDirectory(initialPath);
  const [currentPath, setCurrentPath] = useState(startDir);
  const [files, setFiles] = useState<FileEntry[]>([]);
  const [selectedIndex, setSelectedIndex] = useState(0);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [manualInput, setManualInput] = useState(false);
  const [inputValue, setInputValue] = useState(initialPath || process.cwd());

  // Search mode
  const [searchMode, setSearchMode] = useState(false);
  const [searchQuery, setSearchQuery] = useState('');

  // Get terminal size with resize support
  const { width: terminalWidth, height: terminalHeight } = useTerminalSize();

  useEffect(() => {
    loadFiles();
  }, [currentPath]);

  // Reset search when directory changes
  useEffect(() => {
    setSearchQuery('');
    setSearchMode(false);
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

  // Filter files based on search query
  const filteredFiles = useMemo(() => {
    if (!searchQuery) return files;
    return files.filter(f =>
      f.name.toLowerCase().includes(searchQuery.toLowerCase())
    );
  }, [files, searchQuery]);

  // Reset selection if it exceeds filtered list
  useEffect(() => {
    if (selectedIndex >= filteredFiles.length) {
      setSelectedIndex(Math.max(0, filteredFiles.length - 1));
    }
  }, [filteredFiles.length, selectedIndex]);

  useInput((input, key) => {
    // Manual path input mode
    if (manualInput) {
      if (key.escape) {
        setManualInput(false);
      } else if (key.return) {
        onSelect(inputValue);
      }
      return;
    }

    // Search mode input handling
    if (searchMode) {
      if (key.escape) {
        setSearchMode(false);
        setSearchQuery('');
      } else if (key.return) {
        setSearchMode(false);
      } else if (key.backspace || key.delete) {
        setSearchQuery(searchQuery.slice(0, -1));
      } else if (input && input.length === 1 && !key.ctrl && !key.meta) {
        // Allow any character including space in search
        setSearchQuery(searchQuery + input);
      }
      return;
    }

    // Normal navigation mode
    if (key.upArrow) {
      setSelectedIndex(Math.max(0, selectedIndex - 1));
    } else if (key.downArrow) {
      setSelectedIndex(Math.min(filteredFiles.length - 1, selectedIndex + 1));
    } else if (key.return) {
      const selected = filteredFiles[selectedIndex];
      if (selected) {
        if (selected.isDirectory) {
          setCurrentPath(selected.path);
        } else {
          onSelect(selected.path);
        }
      }
    } else if (input === ' ' && selectDirectory) {
      // Space to select current directory
      onSelect(currentPath);
    } else if (input === '/') {
      // Activate search mode
      setSearchMode(true);
      setSearchQuery('');
    } else if (input === 'p') {
      setInputValue(currentPath);
      setManualInput(true);
    } else if (key.escape || key.leftArrow || input === 'q') {
      // Exit on Esc, left arrow, or q
      onCancel();
    }
  });

  // Calculate available height for file list
  const listHeight = Math.max(5, terminalHeight - 8);
  const halfHeight = Math.floor(listHeight / 2);

  // Calculate visible range centered on selected item
  const startIdx = Math.max(0, Math.min(selectedIndex - halfHeight, filteredFiles.length - listHeight));
  const endIdx = Math.min(filteredFiles.length, startIdx + listHeight);
  const visibleFiles = filteredFiles.slice(startIdx, endIdx);

  // Calculate scrollbar position
  const showScrollbar = filteredFiles.length > listHeight;
  const scrollbarHeight = listHeight;
  const scrollThumbSize = Math.max(1, Math.floor((listHeight / filteredFiles.length) * scrollbarHeight));
  const scrollThumbPos = Math.floor((startIdx / Math.max(1, filteredFiles.length - listHeight)) * (scrollbarHeight - scrollThumbSize));

  const displayTitle = title || (selectDirectory ? 'Select Folder' : 'Select File');
  const displayPath = truncatePath(currentPath, terminalWidth - 20);

  if (loading) {
    return (
      <Box flexDirection="column" padding={1}>
        <Text dimColor>Loading...</Text>
      </Box>
    );
  }

  if (manualInput) {
    return (
      <Box flexDirection="column" padding={1}>
        <Text bold color="cyan">Enter path manually:</Text>
        <Box marginTop={1}>
          <TextInput
            value={inputValue}
            onChange={setInputValue}
            focus={true}
          />
        </Box>
        <Box marginTop={1}>
          <Text dimColor>[Enter] Confirm  [Esc] Cancel</Text>
        </Box>
      </Box>
    );
  }

  return (
    <Box flexDirection="column" flexGrow={1}>
      {/* Header: title and path */}
      <Box paddingX={1} marginBottom={1}>
        <Text bold color="cyan">{displayTitle}</Text>
        <Text dimColor> - </Text>
        <Text color="white">{displayPath}</Text>
        {filteredFiles.length > listHeight && (
          <Text dimColor> ({selectedIndex + 1}/{filteredFiles.length})</Text>
        )}
      </Box>

      {/* Search indicator */}
      {searchMode && (
        <Box paddingX={1}>
          <Text color="yellow">/ </Text>
          <Text color="white">{searchQuery}</Text>
          <Text color="yellow">_</Text>
        </Box>
      )}
      {!searchMode && searchQuery && (
        <Box paddingX={1}>
          <Text dimColor>Filter: </Text>
          <Text color="yellow">{searchQuery}</Text>
          <Text dimColor> ({filteredFiles.length} matches)</Text>
        </Box>
      )}

      {/* Error display */}
      {error && (
        <Box paddingX={1}>
          <Text color="red">Error: {error}</Text>
        </Box>
      )}

      {/* File list with scrollbar */}
      <Box flexDirection="row" flexGrow={1}>
        <Box flexDirection="column" flexGrow={1} paddingX={1}>
          {visibleFiles.map((file, index) => {
            const actualIndex = startIdx + index;
            const isSelected = actualIndex === selectedIndex;

            // Single indicator: ▸ for selected, ▶ for folders, ○ for files
            let icon = '  ';
            if (isSelected) {
              icon = '▸ ';
            } else if (file.isDirectory) {
              icon = '▶ ';
            } else {
              icon = '○ ';
            }

            // Check if file matches allowed extensions
            const fileMatchesExtension = (() => {
              if (file.isDirectory || file.name === '..') return true;
              if (!extensions || extensions.length === 0) return true;
              const ext = path.extname(file.name).toLowerCase().replace('.', '');
              return extensions.some(e => e.toLowerCase().replace('.', '') === ext);
            })();

            const textColor = getFileTextColor(isSelected, file.isDirectory, fileMatchesExtension);
            const iconColor = getFileIconColor(isSelected, file.isDirectory, fileMatchesExtension);

            // Truncate filename if too long (reserve space for icons and scrollbar)
            const maxNameLen = terminalWidth - 8;
            const displayName = truncatePath(file.name, maxNameLen);

            return (
              <Box key={file.path}>
                <Text color={iconColor}>
                  {icon}
                </Text>
                <Text color={textColor}>
                  {displayName}
                </Text>
              </Box>
            );
          })}
          {filteredFiles.length === 0 && (
            <Text dimColor>No matching files</Text>
          )}
        </Box>

        {/* Scrollbar */}
        {showScrollbar && (
          <Box flexDirection="column">
            {Array.from({ length: scrollbarHeight }).map((_, i) => {
              const isThumb = i >= scrollThumbPos && i < scrollThumbPos + scrollThumbSize;
              return (
                <Text key={i} color="gray">
                  {isThumb ? '▓' : '░'}
                </Text>
              );
            })}
          </Box>
        )}
      </Box>

    </Box>
  );
}

export default FileBrowser;
