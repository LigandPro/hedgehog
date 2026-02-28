import React, { useEffect, useMemo, useState } from 'react';
import { Box, Text, useApp, useInput } from 'ink';
import { useStore } from '../store/index.js';
import { useTerminalSize } from '../hooks/useTerminalSize.js';
import { useTheme } from '../theme/context.js';
import { APP_VERSION } from '../constants/version.js';
import { buildExitSummary } from '../utils/exit-summary.js';

interface SlashCommand {
  id: string;
  command: string;
  title: string;
  subtitle: string;
  keywords: string[];
  action: () => void;
}

const MAX_VISIBLE_COMMANDS = 12;
function truncate(value: string, maxLen: number): string {
  if (maxLen <= 0) return '';
  if (value.length <= maxLen) return value;
  if (maxLen <= 1) return '…';
  return `${value.slice(0, maxLen - 1)}…`;
}

export function CommandMenu(): React.ReactElement | null {
  const { exit } = useApp();
  const { theme } = useTheme();
  const { width, height } = useTerminalSize();

  const showCommandMenu = useStore((state) => state.showCommandMenu);
  const setShowCommandMenu = useStore((state) => state.setShowCommandMenu);
  const commandMenuSeed = useStore((state) => state.commandMenuSeed);
  const setCommandMenuSeed = useStore((state) => state.setCommandMenuSeed);
  const setShowThemeMenu = useStore((state) => state.setShowThemeMenu);
  const setShowHelp = useStore((state) => state.setShowHelp);
  const showHelp = useStore((state) => state.showHelp);
  const setScreen = useStore((state) => state.setScreen);
  const screen = useStore((state) => state.screen);
  const goBack = useStore((state) => state.goBack);
  const isRunning = useStore((state) => state.isRunning);
  const currentJobId = useStore((state) => state.currentJobId);
  const selectedJobId = useStore((state) => state.selectedJobId);
  const jobHistory = useStore((state) => state.jobHistory);
  const setExitSummary = useStore((state) => state.setExitSummary);
  const showToast = useStore((state) => state.showToast);

  const [selectedIndex, setSelectedIndex] = useState(0);
  const [query, setQuery] = useState('/');

  const closeMenu = () => {
    setShowCommandMenu(false);
    setSelectedIndex(0);
    setQuery('/');
    setCommandMenuSeed('/');
  };

  const commands = useMemo<SlashCommand[]>(
    () => [
      {
        id: 'new',
        command: '/new',
        title: 'New session',
        subtitle: 'Start a new pipeline run',
        keywords: ['run', 'start', 'wizard'],
        action: () => {
          setScreen('wizardInputSelection');
          closeMenu();
        },
      },
      {
        id: 'history',
        command: '/history',
        title: 'History',
        subtitle: 'Open recent jobs',
        keywords: ['jobs', 'runs', 'previous'],
        action: () => {
          setScreen('history');
          closeMenu();
        },
      },
      {
        id: 'config',
        command: '/config',
        title: 'Configuration',
        subtitle: 'Open main configuration',
        keywords: ['settings', 'options'],
        action: () => {
          setScreen('configMain');
          closeMenu();
        },
      },
      {
        id: 'theme',
        command: '/theme',
        title: 'Switch theme',
        subtitle: 'Open theme selector',
        keywords: ['theme', 'colors', 'style'],
        action: () => {
          closeMenu();
          setShowThemeMenu(true);
        },
      },
      {
        id: 'help',
        command: '/help',
        title: showHelp ? 'Hide help' : 'Show help',
        subtitle: 'Toggle shortcuts overlay',
        keywords: ['shortcuts', 'keys'],
        action: () => {
          setShowHelp(!showHelp);
          closeMenu();
        },
      },
      {
        id: 'status',
        command: '/status',
        title: 'Status',
        subtitle: 'Open running pipeline status',
        keywords: ['pipeline', 'progress'],
        action: () => {
          if (isRunning) {
            setScreen('pipelineRunner');
          } else {
            showToast('info', 'No active pipeline run');
          }
          closeMenu();
        },
      },
      {
        id: 'back',
        command: '/back',
        title: 'Back',
        subtitle: 'Return to previous screen',
        keywords: ['prev', 'navigation'],
        action: () => {
          goBack();
          closeMenu();
        },
      },
      {
        id: 'exit',
        command: '/exit',
        title: 'Exit app',
        subtitle: 'Return to shell',
        keywords: ['quit', 'close'],
        action: () => {
          closeMenu();
          setExitSummary(buildExitSummary({
            currentJobId,
            selectedJobId,
            jobHistory,
            screen,
          }));
          exit();
        },
      },
    ],
    [
      currentJobId,
      exit,
      goBack,
      isRunning,
      jobHistory,
      selectedJobId,
      screen,
      setExitSummary,
      setScreen,
      setShowHelp,
      setShowThemeMenu,
      showHelp,
      showToast,
    ],
  );

  const filteredCommands = useMemo(() => {
    const normalized = (query.startsWith('/') ? query.slice(1) : query).trim().toLowerCase();
    if (!normalized) return commands;

    return commands.filter((item) => {
      const commandName = item.command.slice(1).toLowerCase();
      if (commandName.includes(normalized)) return true;
      if (item.title.toLowerCase().includes(normalized)) return true;
      if (item.subtitle.toLowerCase().includes(normalized)) return true;
      return item.keywords.some((word) => word.toLowerCase().includes(normalized));
    });
  }, [commands, query]);

  useEffect(() => {
    if (!showCommandMenu) return;
    setSelectedIndex(0);
    setQuery(commandMenuSeed.startsWith('/') ? commandMenuSeed : `/${commandMenuSeed}`);
    setCommandMenuSeed('/');
  }, [commandMenuSeed, setCommandMenuSeed, showCommandMenu]);

  useEffect(() => {
    if (filteredCommands.length === 0) {
      setSelectedIndex(0);
      return;
    }
    setSelectedIndex((prev) => {
      if (prev < filteredCommands.length) return prev;
      return filteredCommands.length - 1;
    });
  }, [filteredCommands.length]);

  useInput(
    (input, key) => {
      if (!showCommandMenu) return;

      if (key.escape) {
        closeMenu();
        return;
      }

      if (key.upArrow || input === 'k') {
        if (filteredCommands.length === 0) return;
        setSelectedIndex((prev) => (prev <= 0 ? filteredCommands.length - 1 : prev - 1));
        return;
      }

      if (key.downArrow || input === 'j') {
        if (filteredCommands.length === 0) return;
        setSelectedIndex((prev) => (prev >= filteredCommands.length - 1 ? 0 : prev + 1));
        return;
      }

      if (key.tab) {
        const selected = filteredCommands[selectedIndex];
        if (selected) {
          setQuery(selected.command);
        }
        return;
      }

      if (key.return || key.rightArrow) {
        const selected = filteredCommands[selectedIndex];
        if (selected) {
          selected.action();
        }
        return;
      }

      if (key.backspace || key.delete) {
        setQuery((prev) => {
          if (prev.length <= 1) return '/';
          return prev.slice(0, -1);
        });
        return;
      }

      if (key.ctrl || key.meta) {
        return;
      }

      if (!input || input.length === 0) {
        return;
      }

      if (!/^[a-zA-Z0-9:_/-]+$/.test(input)) {
        return;
      }

      setQuery((prev) => {
        const base = prev.startsWith('/') ? prev : `/${prev}`;
        let next = base;
        for (const ch of input.toLowerCase()) {
          if (ch === '/') {
            next = '/';
            continue;
          }
          if (/^[a-z0-9:_-]$/.test(ch)) {
            next += ch;
          }
        }
        return next;
      });
      setSelectedIndex(0);
    },
    { isActive: showCommandMenu },
  );

  if (!showCommandMenu) {
    return null;
  }

  const panelWidth = Math.max(1, Math.min(120, width - 4));
  const listWidth = Math.max(1, panelWidth - 2);
  const rowContentWidth = Math.max(1, listWidth - 2);
  const visibleCommands = filteredCommands.slice(0, MAX_VISIBLE_COMMANDS);
  const commandColumnWidth = Math.max(1, Math.min(14, Math.max(1, rowContentWidth - 1)));
  const descriptionColumnWidth = Math.max(0, rowContentWidth - commandColumnWidth - 1);
  const inputValue = query.startsWith('/') ? query.slice(1) : query;

  return (
    <Box
      position="absolute"
      width={width}
      height={height}
      flexDirection="column"
    >
      <Box flexGrow={1} />

      <Box justifyContent="center" paddingX={2}>
        <Box width={panelWidth} flexDirection="column">
          {visibleCommands.length > 0 ? (
            visibleCommands.map((item, index) => {
              const active = index === selectedIndex;
              const rowBackground = active ? theme.palette.panelStrong : undefined;
              const commandLabel = truncate(item.command, commandColumnWidth);
              const description = truncate(
                `${item.title} - ${item.subtitle}`,
                descriptionColumnWidth,
              );

              return (
                <Box key={item.id} paddingX={1}>
                  <Text
                    color={active ? theme.palette.text : theme.palette.textMuted}
                    backgroundColor={rowBackground}
                  >
                    {truncate(commandLabel.padEnd(commandColumnWidth, ' '), commandColumnWidth)}
                  </Text>
                  <Text
                    color={theme.palette.textMuted}
                    backgroundColor={rowBackground}
                  >
                    {descriptionColumnWidth > 0 ? ` ${description}` : ''}
                  </Text>
                </Box>
              );
            })
          ) : (
            <Box>
              <Text color={theme.palette.textMuted}>No commands match "{query}"</Text>
            </Box>
          )}

          <Box marginTop={1} paddingX={1}>
            <Text color={theme.palette.textMuted}>/</Text>
            <Text color={theme.palette.text}>{inputValue}</Text>
            <Text color={theme.palette.textMuted}>█</Text>
          </Box>
        </Box>
      </Box>

      <Box flexGrow={1} />

      <Box width={width} justifyContent="center">
        <Text color="white">ctrl+p</Text>
        <Text color={theme.palette.textMuted}> commands</Text>
      </Box>
      <Box width={width} justifyContent="flex-end" paddingRight={2}>
        <Text color={theme.palette.textMuted}>{APP_VERSION}</Text>
      </Box>
    </Box>
  );
}

export default CommandMenu;
