import React, { useState, useEffect, useMemo } from 'react';
import { Box, Text, useInput } from 'ink';
import TextInput from 'ink-text-input';
import path from 'path';
import fs from 'node:fs';
import { getBridge } from '../services/python-bridge.js';
import { useTerminalSize } from '../hooks/useTerminalSize.js';
import { useInputLock } from '../hooks/useInputLock.js';
import { useTheme } from '../theme/context.js';
import { LineFill } from './LineFill.js';
import type { ThemeDefinition } from '../theme/themes.js';

interface FileBrowserProps {
  initialPath: string;
  extensions?: string[];
  selectDirectory?: boolean;
  onSelect: (path: string) => void;
  onCancel: () => void;
  title?: string;
  startInPathEdit?: boolean;
}

interface FileEntry {
  name: string;
  path: string;
  isDirectory: boolean;
}

type InteractionMode = 'browse' | 'pathEdit' | 'search';

/**
 * Get a valid starting directory from a path.
 * If path looks like a file (has extension), return parent directory.
 * If path is empty or invalid, return cwd.
 */
function getStartingDirectory(inputPath: string): string {
  if (!inputPath || inputPath.trim() === '') {
    return process.cwd();
  }

  const expanded = expandHomePath(inputPath.trim());
  const normalizedExpanded = path.normalize(expanded);
  if (
    isExistingDirectory(normalizedExpanded)
    || isExistingDirectory(path.resolve(process.cwd(), normalizedExpanded))
  ) {
    return normalizedExpanded;
  }

  // Check if path looks like a file (has extension in the last segment)
  const basename = path.basename(normalizedExpanded);
  const hasExtension = basename.includes('.') && !basename.startsWith('.');

  if (hasExtension) {
    // It's likely a file, use parent directory
    return path.dirname(normalizedExpanded);
  }

  return normalizedExpanded;
}

function expandHomePath(inputPath: string): string {
  if (!inputPath.startsWith('~')) {
    return inputPath;
  }
  const home = process.env.HOME || '';
  if (!home) {
    return inputPath;
  }
  if (inputPath === '~') {
    return home;
  }
  if (inputPath.startsWith('~/')) {
    return path.join(home, inputPath.slice(2));
  }
  return inputPath;
}

function resolveInputPath(inputPath: string, basePath: string): string {
  const trimmed = inputPath.trim();
  if (!trimmed) {
    return basePath;
  }
  const expanded = expandHomePath(trimmed);
  if (path.isAbsolute(expanded)) {
    return path.normalize(expanded);
  }

  const normalizedExpanded = path.normalize(expanded);
  const normalizedBase = path.normalize(basePath);

  const hasBasePrefix = (value: string, prefix: string): boolean => (
    value === prefix || value.startsWith(`${prefix}${path.sep}`)
  );

  if (
    normalizedBase &&
    normalizedBase !== '.' &&
    hasBasePrefix(normalizedExpanded, normalizedBase)
  ) {
    return normalizedExpanded;
  }

  if (path.isAbsolute(normalizedBase)) {
    const relativeBase = path.normalize(path.relative(process.cwd(), normalizedBase));
    if (
      relativeBase &&
      relativeBase !== '.' &&
      !relativeBase.startsWith('..') &&
      !path.isAbsolute(relativeBase) &&
      hasBasePrefix(normalizedExpanded, relativeBase)
    ) {
      return normalizedExpanded;
    }
  }

  return path.resolve(basePath, normalizedExpanded);
}

function isExistingDirectory(candidatePath: string): boolean {
  try {
    return fs.existsSync(candidatePath) && fs.statSync(candidatePath).isDirectory();
  } catch {
    return false;
  }
}

function findNearestExistingDirectory(candidatePath: string): string | null {
  let current = path.normalize(candidatePath);
  const root = path.parse(current).root || '/';

  while (true) {
    if (isExistingDirectory(current)) {
      return current;
    }
    if (current === root) {
      return null;
    }
    const parent = path.dirname(current);
    if (!parent || parent === current) {
      return null;
    }
    current = parent;
  }
}

function getPathTailFilter(inputValue: string, resolvedPath: string): string {
  const trimmed = inputValue.trim();
  if (!trimmed) {
    return '';
  }
  if (/[\\/]$/.test(trimmed)) {
    return '';
  }

  const basename = path.basename(path.normalize(resolvedPath));
  if (!basename || basename === '.' || basename === '..') {
    return '';
  }
  if (isExistingDirectory(resolvedPath)) {
    return '';
  }
  return basename;
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
function getFileTextColor(
  theme: ThemeDefinition,
  isSelected: boolean,
  isDirectory: boolean,
  matchesExtension: boolean,
): string {
  if (isSelected) return theme.palette.text;
  if (isDirectory) return theme.palette.info;
  if (matchesExtension) return theme.palette.success;
  return theme.palette.textMuted;
}

/**
 * Get icon color for file entry based on selection state and type.
 */
function getFileIconColor(
  theme: ThemeDefinition,
  isSelected: boolean,
  isDirectory: boolean,
  matchesExtension: boolean,
): string {
  if (isSelected) return theme.palette.primary;
  if (isDirectory) return theme.palette.info;
  if (matchesExtension) return theme.palette.success;
  return theme.palette.textMuted;
}

export function FileBrowser({
  initialPath,
  extensions,
  selectDirectory = false,
  onSelect,
  onCancel,
  title,
  startInPathEdit = false,
}: FileBrowserProps): React.ReactElement {
  const startDir = getStartingDirectory(initialPath);
  const [currentPath, setCurrentPath] = useState(startDir);
  const [files, setFiles] = useState<FileEntry[]>([]);
  const [selectedIndex, setSelectedIndex] = useState(0);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [interactionMode, setInteractionMode] = useState<InteractionMode>(
    startInPathEdit ? 'pathEdit' : 'browse'
  );
  const [inputValue, setInputValue] = useState(initialPath || process.cwd());
  const [inputBasePath, setInputBasePath] = useState(startDir);
  const [searchQuery, setSearchQuery] = useState('');
  const { theme } = useTheme();

  // Get terminal size with resize support
  const { width: terminalWidth, height: terminalHeight } = useTerminalSize();
  const contentWidth = Math.max(8, terminalWidth - 2);
  const isPathEditMode = interactionMode === 'pathEdit';
  const isSearchMode = interactionMode === 'search';
  useInputLock(isPathEditMode || isSearchMode);

  const resolvedInputPath = useMemo(
    () => resolveInputPath(inputValue, inputBasePath),
    [inputValue, inputBasePath]
  );

  const previewDir = useMemo(
    () => findNearestExistingDirectory(getStartingDirectory(resolvedInputPath)),
    [resolvedInputPath]
  );

  const pathTailFilter = useMemo(
    () => (isPathEditMode ? getPathTailFilter(inputValue, resolvedInputPath) : ''),
    [inputValue, isPathEditMode, resolvedInputPath]
  );

  const activeFilterQuery = isPathEditMode ? pathTailFilter : searchQuery;

  useEffect(() => {
    loadFiles();
  }, [currentPath]);

  useEffect(() => {
    if (!isPathEditMode) {
      return;
    }
    if (previewDir && previewDir !== currentPath) {
      setCurrentPath(previewDir);
    }
  }, [isPathEditMode, previewDir, currentPath]);

  // Reset search when directory changes
  useEffect(() => {
    setSearchQuery('');
    setInteractionMode((mode) => (mode === 'search' ? 'browse' : mode));
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
      setFiles([{ name: '..', path: path.join(currentPath, '..'), isDirectory: true }, ...result]);
      setSelectedIndex(0);
    } catch (err) {
      setError(String(err));
    } finally {
      setLoading(false);
    }
  };

  // Filter files based on search query or manual path tail
  const filteredFiles = useMemo(() => {
    if (!activeFilterQuery) return files;
    return files.filter(f =>
      f.name.toLowerCase().includes(activeFilterQuery.toLowerCase())
    );
  }, [files, activeFilterQuery]);

  // Reset selection if it exceeds filtered list
  useEffect(() => {
    if (selectedIndex >= filteredFiles.length) {
      setSelectedIndex(Math.max(0, filteredFiles.length - 1));
    }
  }, [filteredFiles.length, selectedIndex]);

  const enterPathEditMode = (seed: string): void => {
    setInputBasePath(currentPath);
    setInputValue(seed);
    setSearchQuery('');
    setInteractionMode('pathEdit');
  };

  const enterSearchMode = (): void => {
    setSearchQuery('');
    setInteractionMode('search');
  };

  useInput((input, key) => {
    // Manual path input mode
    if (isPathEditMode) {
      if (key.escape) {
        setInteractionMode('browse');
      } else if (key.return) {
        if (!selectDirectory && isExistingDirectory(resolvedInputPath)) {
          setCurrentPath(resolvedInputPath);
          setInteractionMode('browse');
          return;
        }
        if (selectDirectory) {
          // Keep user-entered directory path as-is to allow creating new folders.
          onSelect(resolvedInputPath);
          return;
        }
        onSelect(resolvedInputPath);
      }
      return;
    }

    // Search mode input handling
    if (isSearchMode) {
      if (key.escape) {
        setInteractionMode('browse');
        setSearchQuery('');
      } else if (key.return) {
        setInteractionMode('browse');
      } else if (key.backspace || key.delete) {
        setSearchQuery((prev) => prev.slice(0, -1));
      } else if (input && input.length === 1 && !key.ctrl && !key.meta) {
        // Allow any character including space in search
        setSearchQuery((prev) => prev + input);
      }
      return;
    }

    // Normal navigation mode
    if (key.upArrow) {
      if (selectedIndex === 0) {
        enterPathEditMode(currentPath);
        return;
      }
      setSelectedIndex(Math.max(0, selectedIndex - 1));
    } else if (key.downArrow) {
      if (filteredFiles.length === 0) {
        return;
      }
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
    } else if (input === ' ') {
      if (selectDirectory) {
        // Space to select current directory
        onSelect(currentPath);
      } else {
        // Space opens quick search in file mode
        enterSearchMode();
      }
    } else if (key.ctrl && input.toLowerCase() === 'f') {
      // Activate search mode
      enterSearchMode();
    } else if (key.rightArrow || input === 'e') {
      enterPathEditMode(currentPath);
    } else if (key.escape || key.leftArrow) {
      // Exit on Esc or left arrow
      onCancel();
    }
  });

  // Reserve rows for stable layout across header/filter/error/footer wrappers.
  const reservedRows = 12;
  const maxVisibleRows = Math.max(5, terminalHeight - reservedRows);
  const listHeight = Math.max(
    5,
    Math.min(filteredFiles.length || 1, maxVisibleRows),
  );
  const halfHeight = Math.floor(listHeight / 2);

  // Calculate visible range centered on selected item
  const startIdx = Math.max(0, Math.min(selectedIndex - halfHeight, filteredFiles.length - listHeight));
  const endIdx = Math.min(filteredFiles.length, startIdx + listHeight);
  const visibleFiles = filteredFiles.slice(startIdx, endIdx);

  // Calculate scrollbar position
  const showScrollbar = filteredFiles.length > listHeight;
  const scrollbarHeight = listHeight;
  const scrollThumbSize = Math.max(1, Math.floor((listHeight / Math.max(1, filteredFiles.length)) * scrollbarHeight));
  const scrollThumbPos = showScrollbar
    ? Math.floor((startIdx / Math.max(1, filteredFiles.length - listHeight)) * (scrollbarHeight - scrollThumbSize))
    : 0;

  const displayTitle = title || (selectDirectory ? 'Select Folder' : 'Select File');
  const displayPath = truncatePath(currentPath, terminalWidth - 20);

  if (loading) {
    return (
      <Box flexDirection="column" padding={1}>
        <Text color={theme.palette.textMuted}>Loading...</Text>
      </Box>
    );
  }

  if (isPathEditMode) {
    const previewPathForDisplay = previewDir || currentPath;
    const manualDisplayPath = truncatePath(previewPathForDisplay, terminalWidth - 25);
    return (
      <Box flexDirection="column" padding={1}>
        <Text bold color={theme.palette.primary}>Edit path:</Text>
        <Box marginTop={1}>
          <TextInput
            value={inputValue}
            onChange={setInputValue}
            focus={true}
          />
        </Box>
        <Box marginTop={1}>
          <Text color={theme.palette.textMuted}>Preview folder: </Text>
          <Text color={theme.palette.text}>{manualDisplayPath}</Text>
        </Box>
        {pathTailFilter && (
          <Box marginTop={1}>
            <Text color={theme.palette.textMuted}>Name filter: </Text>
            <Text color={theme.palette.accent}>{pathTailFilter}</Text>
          </Box>
        )}
        <Box marginTop={1}>
          <Text color={theme.palette.textMuted}>[Enter] Apply  [Esc] Cancel</Text>
        </Box>
        <Box marginTop={1} flexDirection="column">
          <Text color={theme.palette.textMuted}>Contents:</Text>
          {filteredFiles.slice(0, 8).map((file) => (
            <Box key={file.path} width={contentWidth}>
              <Text color={file.isDirectory ? theme.palette.info : theme.palette.textMuted} wrap="truncate-end">
                {file.isDirectory ? '▶ ' : '○ '}
                {file.name}
              </Text>
              <LineFill width={contentWidth} />
            </Box>
          ))}
          {filteredFiles.length > 8 && (
            <Text color={theme.palette.textMuted}>... and {filteredFiles.length - 8} more</Text>
          )}
        </Box>
      </Box>
    );
  }

  return (
    <Box flexDirection="column" flexGrow={1}>
      {/* Header: title and path */}
      <Box paddingX={1} marginBottom={1} width={terminalWidth}>
        <Text bold color={theme.palette.primary}>{displayTitle}</Text>
        <Text color={theme.palette.textMuted}> - </Text>
        <Text color={theme.palette.text}>{displayPath}</Text>
        {filteredFiles.length > listHeight && (
          <Text color={theme.palette.textMuted}> ({selectedIndex + 1}/{filteredFiles.length})</Text>
        )}
        <LineFill width={terminalWidth} />
      </Box>

      {/* Search indicator */}
      {isSearchMode && (
        <Box paddingX={1} width={terminalWidth}>
          <Text color={theme.palette.accent}>Ctrl+F </Text>
          <Text color={theme.palette.text}>{searchQuery}</Text>
          <Text color={theme.palette.accent}>_</Text>
          <LineFill width={terminalWidth} />
        </Box>
      )}
      {!isSearchMode && activeFilterQuery && (
        <Box paddingX={1} width={terminalWidth}>
          <Text color={theme.palette.textMuted}>Filter: </Text>
          <Text color={theme.palette.accent}>{activeFilterQuery}</Text>
          <Text color={theme.palette.textMuted}> ({filteredFiles.length} matches)</Text>
          <LineFill width={terminalWidth} />
        </Box>
      )}

      {/* Error display */}
      {error && (
        <Box paddingX={1} width={terminalWidth}>
          <Text color={theme.palette.error}>Error: {error}</Text>
          <LineFill width={terminalWidth} />
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

            const textColor = getFileTextColor(theme, isSelected, file.isDirectory, fileMatchesExtension);
            const iconColor = getFileIconColor(theme, isSelected, file.isDirectory, fileMatchesExtension);

            // Truncate filename if too long (reserve space for icons and scrollbar)
            const maxNameLen = terminalWidth - 8;
            const displayName = truncatePath(file.name, maxNameLen);

            return (
              <Box key={file.path} width={contentWidth}>
                <Text color={iconColor}>
                  {icon}
                </Text>
                <Text color={textColor}>
                  {displayName}
                </Text>
                <LineFill width={contentWidth} />
              </Box>
            );
          })}
          {filteredFiles.length === 0 && (
            <Text color={theme.palette.textMuted}>No matching files</Text>
          )}
        </Box>

        {/* Scrollbar */}
        {showScrollbar && (
          <Box flexDirection="column">
            {Array.from({ length: scrollbarHeight }).map((_, i) => {
              const isThumb = i >= scrollThumbPos && i < scrollThumbPos + scrollThumbSize;
              return (
                <Text key={i} color={theme.palette.border}>
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
