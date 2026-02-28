#!/usr/bin/env node
import React from 'react';
import { render } from 'ink';
import { App } from './App.js';
import { logger } from './utils/logger.js';
import { ThemeProvider } from './theme/context.js';
import { useStore } from './store/index.js';

const ENTER_ALT_SCREEN = '\x1B[?1049h';
const EXIT_ALT_SCREEN = '\x1B[?1049l';
const RESET_SGR = '\x1B[0m';
const DEFAULT_BG = '\x1B[48;2;0;0;0m';
const DEFAULT_FG = '\x1B[38;2;230;230;230m';

let restoredTerminal = false;
const restoreTerminal = () => {
  if (restoredTerminal) return;
  restoredTerminal = true;
  if (process.stdout.isTTY) {
    process.stdout.write(`${RESET_SGR}${EXIT_ALT_SCREEN}`);
  }
};

if (process.stdout.isTTY) {
  process.stdout.write(`${ENTER_ALT_SCREEN}${DEFAULT_BG}${DEFAULT_FG}`);
}

process.once('SIGINT', restoreTerminal);
process.once('SIGTERM', restoreTerminal);
process.once('exit', restoreTerminal);

logger.info('Starting Hedgehog TUI');

const { waitUntilExit } = render(
  <ThemeProvider>
    <App />
  </ThemeProvider>,
  { exitOnCtrlC: false },
);

waitUntilExit().then(() => {
  restoreTerminal();
  const state = useStore.getState();
  if (state.exitSummary) {
    process.stdout.write(`${state.exitSummary}\n`);
    state.setExitSummary(null);
  }
  logger.info('Hedgehog TUI exited');
  process.exit(0);
}).catch((error) => {
  restoreTerminal();
  logger.error('TUI error:', error);
  process.exit(1);
});
