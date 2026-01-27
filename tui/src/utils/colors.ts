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
  version: chalk.dim.gray,

  // Debug mode
  debug: chalk.yellow,

  // Footer shortcuts
  shortcutKey: chalk.bold.cyan,
  shortcutLabel: chalk.white,
  shortcutDisabled: chalk.dim.gray,
};

export const icons = {
  success: chalk.green('✓'),
  error: chalk.red('✗'),
  pending: chalk.gray('○'),
  running: chalk.yellow('●'),
  completed: chalk.green('●'),
  skipped: chalk.gray('○'),
  cancelled: chalk.gray('○'),
  arrow: chalk.cyan('→'),
  bullet: chalk.gray('•'),
  hedgehog: '*',
  check: '✓',
  cross: '✗',
  circle: '○',
  filledCircle: '●',
};

export default colors;
