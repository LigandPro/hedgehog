# 🦔 HEDGEHOG
**Hierarchical Evaluation of Drug GEnerators tHrOugh riGorous filtration**

[![PyPI version](https://badge.fury.io/py/hedgehog.svg)](https://pypi.org/project/hedgehog/)
[![CI](https://github.com/LigandPro/hedgehog/actions/workflows/ci.yaml/badge.svg)](https://github.com/LigandPro/hedgehog/actions/workflows/ci.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

![HEDGEHOG Pipeline](docs/public/pipeline_structure.png)

HEDGEHOG is a benchmark pipeline for evaluating generated molecules.
It applies a staged filtering workflow (standardization, descriptors, structural filters, synthesis checks, docking, and post-docking validation) and generates an interactive HTML report.

## Quick Start

### Install from PyPI

```bash
pip install hedgehog
hedgehog --help
```

### Install from source

```bash
git clone https://github.com/LigandPro/hedgehog.git
cd hedgehog
uv sync
```

### First run

```bash
uv run hedgehog run

# Optional: auto-install missing optional tools without prompts
uv run hedgehog run --auto-install
```

## Common Commands

```bash
# Full pipeline
uv run hedgehog run

# Run with your own molecules
uv run hedgehog run --mols input/my_molecules.csv

# Run a single stage
uv run hedgehog run --stage descriptors

# Run docking with a live progress bar
uv run hedgehog run --stage docking --progress

# Run docking without progress bar (default)
uv run hedgehog run --stage docking

# Regenerate report for an existing run
uv run hedgehog report results/run_10

# Show stages / version
uv run hedgehog info
uv run hedgehog version

# Launch terminal UI
uv run hedgehog tui
```

Progress bar behavior in CLI runs:
- Enabled: add `--progress`
- Disabled: omit `--progress` (default)

## Documentation

For full details, use the documentation instead of this README:

- [Introduction](docs/content/index.mdx)
- [Getting Started](docs/content/getting-started.mdx)
- [CLI Reference](docs/content/cli.mdx)
- [Pipeline Stages](docs/content/pipeline/index.mdx)
- [Configuration](docs/content/configuration/index.mdx)
- [Reporting](docs/content/reporting/index.mdx)
- [Advanced Topics](docs/content/advanced/architecture.mdx)
- [TUI README](tui/README.md)

To run the docs site locally:

```bash
cd docs
pnpm install
pnpm dev
```

## License

MIT
