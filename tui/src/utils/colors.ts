import chalk from 'chalk';

export const colors = {
  // Primary: Cyan (Hedgehog brand)
  primary: chalk.cyan,
  primaryBold: chalk.bold.cyan,

  // Secondary: Magenta
  secondary: chalk.magenta,

  // Accent: Yellow
  accent: chalk.yellow,

  // Status
  success: chalk.green,
  error: chalk.red,
  warning: chalk.yellow,
  info: chalk.blue,

  // UI
  border: chalk.gray,
  muted: chalk.dim,
  highlight: chalk.bgCyan.black,

  // Brand
  brand: chalk.bold.cyan,
  version: chalk.dim.cyan,
};

export const icons = {
  success: chalk.green('✓'),
  error: chalk.red('✗'),
  pending: chalk.gray('○'),
  running: chalk.yellow('●'),
  arrow: chalk.cyan('→'),
  bullet: chalk.gray('•'),
  hedgehog: '🦔',
};

export default colors;
