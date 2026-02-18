# Hedgehog TUI

Terminal User Interface for the Hedgehog molecular design pipeline.

## Features

- **Configuration Editor**: Edit pipeline settings for all stages (descriptors, filters, synthesis, docking)
- **Pipeline Wizard**: Step-by-step guided setup for new runs
- **Pipeline Runner**: Execute pipeline with real-time progress tracking
- **File Browser**: Navigate and select input files
- **History**: View past pipeline runs and results

## Requirements

- Node.js 18+
- Python 3.10+ with Hedgehog installed

## Installation

```bash
cd tui
npm install
```

## Running

### Development mode (with auto-reload)

```bash
npm run dev
```

### Production build

```bash
npm run build
npm start
```

### From project root

```bash
cd tui
npm run tui
```

### Smoke checks

From project root:

```bash
# Build + startup smoke (starts TUI in PTY and exits automatically)
uv run python scripts/check_pipeline.py --mode quick
```

This is the recommended way to verify that the TUI can start and connect to
the Python backend in automation-friendly environments.

## Project Structure

```
tui/
├── src/
│   ├── components/     # Reusable UI components
│   │   ├── Header.tsx
│   │   ├── Footer.tsx
│   │   ├── FileBrowser.tsx
│   │   ├── ProgressBar.tsx
│   │   └── ...
│   ├── screens/        # Main application screens
│   │   ├── Welcome.tsx
│   │   ├── ConfigMain.tsx
│   │   ├── ConfigDescriptors.tsx
│   │   ├── ConfigFilters.tsx
│   │   ├── ConfigSynthesis.tsx
│   │   ├── ConfigDocking.tsx
│   │   ├── PipelineRunner.tsx
│   │   ├── History.tsx
│   │   └── wizard/     # Wizard flow screens
│   ├── hooks/          # Custom React hooks
│   ├── services/       # Backend communication
│   │   ├── python-bridge.ts
│   │   └── rpc-client.ts
│   ├── store/          # Zustand state management
│   ├── types/          # TypeScript types
│   ├── utils/          # Utility functions
│   ├── App.tsx         # Main application component
│   └── index.tsx       # Entry point
├── bin/                # CLI entry point
├── package.json
└── tsconfig.json
```

## Key Bindings

| Key | Action |
|-----|--------|
| `↑/↓` | Navigate lists |
| `Enter` | Select/Confirm |
| `Esc` | Go back |
| `Tab` | Next field |
| `/` | Search (in lists) |
| `q` | Quit |
| `?` | Help |

## Architecture

The TUI communicates with a Python backend (`src/hedgehog/tui_backend/`) via JSON-RPC over stdin/stdout. The backend handles:

- Configuration file read/write
- Pipeline execution
- File system operations
- Run history management

## Technologies

- [Ink](https://github.com/vadimdemedes/ink) - React for CLI
- [Zustand](https://github.com/pmndrs/zustand) - State management
- TypeScript
