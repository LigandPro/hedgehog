#!/usr/bin/env node
import React from 'react';
import { render } from 'ink';
import { App } from './App.js';
import { logger } from './utils/logger.js';

// Clear terminal screen to hide npm startup messages
process.stdout.write('\x1B[2J\x1B[0f');

logger.info('Starting Hedgehog TUI');

const { waitUntilExit } = render(<App />);

waitUntilExit().then(() => {
  logger.info('Hedgehog TUI exited');
  process.exit(0);
}).catch((error) => {
  logger.error('TUI error:', error);
  process.exit(1);
});
