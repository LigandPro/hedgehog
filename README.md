<p align="center">
  <img src="docs/public/hedgehog-title.png" alt="HEDGEHOG" width="40%" />
</p>
<p align="center">HEDGEHOG: Hierarchical Evaluation of Drug GEnerators tHrOugh riGorous filtration.</p>
<p align="center">
  <a href="https://pypi.org/project/hedgehog/"><img alt="PyPI" src="https://img.shields.io/pypi/v/hedgehog?style=flat-square" /></a>
  <a href="https://github.com/LigandPro/hedgehog/actions/workflows/ci.yaml"><img alt="CI" src="https://img.shields.io/github/actions/workflow/status/LigandPro/hedgehog/ci.yaml?style=flat-square&branch=main" /></a>
  <a href="https://opensource.org/licenses/MIT"><img alt="License: MIT" src="https://img.shields.io/badge/license-MIT-yellow.svg?style=flat-square" /></a>
  <a href="https://www.python.org/downloads/"><img alt="Python 3.10+" src="https://img.shields.io/badge/python-3.10%2B-blue.svg?style=flat-square" /></a>
</p>

<p align="center">
  <img src="docs/public/hedgehog-tui-home.png" alt="HEDGEHOG Terminal UI" width="90%" />
</p>



## Quick Start

HEDGEHOG is a stage-based molecular design evaluation pipeline for:

- molecule preparation
- descriptor calculation
- structural filtering
- retrosynthesis filtering
- docking
- docking pose filtering
- final reports

The full pipeline can require optional external tools and receptor inputs. Start
with the safe smoke run below to verify the Python environment, bundled example
molecules, descriptor calculation, and structural filters before enabling
retrosynthesis or docking.

### Recommended install: source checkout

```bash
git clone https://github.com/LigandPro/hedgehog.git
cd hedgehog
uv sync
```

This is the recommended way to run HEDGEHOG end to end. The repository checkout contains the editable configs, bundled examples, TUI sources, and the `modules/` workspace used by setup commands.

Requirements:

- Python 3.10+
- `uv`
- optional: Node.js >= 18 and npm for the TUI
- optional: AiZynthFinder for retrosynthesis
- optional: GNINA, SMINA, or Matcha for docking

### PyPI install

```bash
python -m pip install hedgehog
hedgehog --help
```

Use the PyPI package only if you already manage your own config files and input paths. The default quick start, `hedgehog setup ...` workflows, and TUI usage are designed around a source checkout.

### First safe run

```bash
uv run hedgehog --stage descriptors --stage struct_filters --force-new
```

This avoids docking and retrosynthesis. Use it as the first validation that the
local environment, bundled examples, descriptor calculation, and structural
filters are working.

### Full pipeline

```bash
uv run hedgehog setup aizynthfinder
uv run hedgehog --auto-install
```

Full pipeline execution may require AiZynthFinder, GNINA/SMINA/Matcha, valid
receptor structures, reference ligands, and enough CPU/GPU resources.

## Input Format

Recommended molecule input is CSV/TSV with a `smiles` header:

```csv
smiles,model_name
CCO,demo
CCN,demo
c1ccccc1,demo
```

Required:

- `smiles`

Optional:

- `model_name` or `name`
- `mol_idx`

If `mol_idx` is missing, HEDGEHOG assigns a stable ID and uses it to join stage
outputs, docking scores, and report data.

## Common Commands

```bash
# Safe smoke run
uv run hedgehog --stage descriptors --stage struct_filters --force-new

# Full pipeline after optional tools are available
uv run hedgehog --auto-install

# Run with your own molecules
uv run hedgehog --mols input/my_molecules.csv

# Run a single stage
uv run hedgehog --stage descriptors

# Run multiple selected stages
uv run hedgehog --stage descriptors --stage struct_filters

# Run docking with a live progress bar
uv run hedgehog --stage docking --progress

# Run docking without progress bar (default)
uv run hedgehog --stage docking

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

## Results

Results are written under the configured output directory, usually as an
auto-numbered run folder:

```text
results/run_N/
├── stages/
├── output/
└── report.html
```

## Documentation

![HEDGEHOG Pipeline](docs/public/pipeline_structure_v2.png)

For full details, use the documentation instead of this README:

- [Introduction](docs/content/index.mdx)
- [Getting Started](docs/content/getting-started.mdx)
- [CLI Reference](docs/content/cli.mdx)
- [TUI](docs/content/tui.mdx)
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
