# 🦔 HEDGEHOG
**Hierarchical Evaluation of Drug GEnerators tHrOugh riGorous filtration**

[![PyPI version](https://badge.fury.io/py/hedgehog.svg)](https://pypi.org/project/hedgehog/)
[![CI](https://github.com/LigandPro/hedgehog/actions/workflows/ci.yaml/badge.svg)](https://github.com/LigandPro/hedgehog/actions/workflows/ci.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

![HEDGEHOG Pipeline](data/imgs/pipeline_structure.png)


Comprehensive benchmark pipeline for evaluating generative models in molecular design.

### Pipeline Stages:

Each stage takes the output of the previous one, progressively filtering the molecule set:

1) **Mol Prep (Datamol)**: salts/solvents & fragments cleanup, largest-fragment selection, metal disconnection, uncharging, tautomer canonicalization, stereochemistry removal → produces standardized “clean” molecules ([molPrep folder](src/hedgehog/stages/molPrep/))
2) **Molecular Descriptors**: 22 physicochemical descriptors (logP, HBD/HBA, TPSA, QED, etc.) → molecules outside thresholds are removed ([descriptors folder](src/hedgehog/stages/descriptors/))
3) **Structural Filters**: 6 criteria with ~2500 SMARTS patterns (PAINS, Glaxo, NIBR, Bredt, etc.) → flagged molecules are removed ([structural filters folder](src/hedgehog/stages/structFilters/))
4) **Synthesis Evaluation**: SA score, SYBA score, AiZynthFinder retrosynthesis → unsynthesizable molecules are removed ([synthesis folder](src/hedgehog/stages/synthesis/))
5) **Molecular Docking**: SMINA and/or GNINA → binding affinity scoring ([docking folder](src/hedgehog/stages/docking/))
6) **Docking Filters**: post-docking pose quality filtering → poor binders are removed
7) **Final Descriptors**: recalculation on the filtered set

Post-pipeline analysis: MolEval generative metrics

## Setup & Run

### Install from PyPI

```bash
python -m pip install hedgehog
hedgehog --help
```

Base install is intentionally lightweight and works on modern Python versions
(including Python 3.13) without optional heavy docking extras.

Optional extras:

```bash
# Legacy PoseCheck backend for docking filters
python -m pip install 'hedgehog[docking-legacy]'

# Shepherd-Score Python dependency only (may be unavailable on some ABIs, e.g. cp313)
python -m pip install 'hedgehog[shepherd]'
```

Recommended Shepherd setup is an isolated worker environment:

```bash
uv run hedgehog setup shepherd-worker --yes
```

### Install from source (recommended for development)

```bash
# Clone repository
git clone https://github.com/LigandPro/hedgehog.git
cd hedgehog


# Install AiZynthFinder (for synthesis stage) - recommended CLI flow
uv run hedgehog setup aizynthfinder

# Legacy helper script (alternative)
./modules/install_aizynthfinder.sh

# Install package with uv
uv sync
```

You are ready to use **🦔 HEDGEHOG** for your purpose!

**Usage**

```bash
# Run full pipeline on a proposed small test data from `data/test/`
uv run hedgehog run

# Alternatively, using the short-name alias:
uv run hedge run

# Run full pipeline on your own molecule file
uv run hedge run --mols data/my_molecules.csv

# Run full pipeline on your own molecule files via glob
uv run hedge run --mols "data/generated/*.csv"

# Run specific stage
uv run hedge run --stage descriptors

# Run a specific stage on your own molecule file
uv run hedge run --stage descriptors --mols data/my_molecules.csv

# Auto-install missing optional external tools during a run
uv run hedge run --auto-install

# Reuse the existing results folder
uv run hedge run --reuse

# Force a fresh results folder for stage reruns
uv run hedge run --stage docking --force-new

# Enable live progress bar in CLI
uv run hedge run --progress

# Regenerate HTML report from an existing run
uv run hedge report results/run_10

# Show pipeline stages and current version
uv run hedge info
uv run hedge version

# Get help
uv run hedge --help
```

**GNINA (CPU/GPU) notes**

GNINA is auto-resolved during the docking stage:
- If `gnina` is already on `PATH`, HEDGEHOG uses it.
- Otherwise it auto-downloads a compatible Linux GNINA binary to `~/.hedgehog/bin/gnina`.

By default, auto-install uses `HEDGEHOG_GNINA_VARIANT=auto` behavior
(prefer CUDA build when NVIDIA GPU is detected, otherwise CPU fallback).
If needed, you can override explicitly:

```bash
export HEDGEHOG_GNINA_VARIANT=auto   # or: gpu
uv run hedge run --stage docking --auto-install
```

GNINA runtime now auto-discovers CUDA/PyTorch libraries from common locations
(including `site-packages/nvidia/*/lib`, active conda env, and `~/miniforge/lib`)
so manual `gnina_ld_library_path` is usually not required.

**Terminal UI (TUI)**

For interactive configuration and pipeline management, use the TUI:
```bash
uv run hedgehog tui
```

If the TUI has not been built yet, the CLI will install/build it automatically on first launch.
You can also launch it directly from the TUI package:

```bash
cd tui
npm run tui
```

See [tui/README.md](tui/README.md) for details and developer workflow.

**Unified verification pipeline**

Use one command entry point for local/CI checks:

```bash
# Quick local smoke (CLI + TUI build + TUI startup/quit in PTY)
uv run python scripts/check_pipeline.py --mode quick

# CI smoke profile (same checks, no full production run)
uv run python scripts/check_pipeline.py --mode ci

# Full local verification (quick checks + full production pipeline run)
uv run python scripts/check_pipeline.py --mode full
```

`--mode full` runs `uv run hedgehog run` with the default production config,
so it can be long-running and requires external stage dependencies (for example
docking/synthesis tooling) to be installed in your local environment.

**Git hooks with Lefthook (recommended)**

Use Lefthook to block commits/pushes that would fail CI:

```bash
# Install Lefthook (macOS)
brew install lefthook

# Register git hooks from lefthook.yml
lefthook install

# Optional: run hooks manually
lefthook run pre-commit
lefthook run pre-push
```

Current local gates:

- `pre-commit`: staged Python formatting/lint (`ruff`) and whitespace checks.
- `pre-push`: repository-wide `ruff` checks, `pytest`, pipeline smoke (`scripts/check_pipeline.py --mode ci`), and docs build (`docs && pnpm build`).

If you need to skip a specific hook command once (not recommended), use `SKIP`:

```bash
SKIP=docs-build git push
```

**Documentation Site**

```bash
cd docs && pnpm install && pnpm dev
```

The docs site is built with [Nextra](https://nextra.site) and available at `http://localhost:3000`.

**HTML Reports**

After each pipeline run, an interactive HTML report is automatically generated as `report.html` in the results folder. The report includes:
- Pipeline summary and molecule retention funnel
- Per-stage statistics and visualizations
- Descriptor distributions
- Filter pass/fail breakdowns
- Synthesis scores and docking results

**Configure your run**
Edit config for each stage in [configs folder](src/hedgehog/configs/) based on metrics you want to calculate.


<!-- ## REINVENT4 fine-tune
To fine-tune REINVENT4 follow these steps:
1. Clone REINVENT4 repository and setup environment:
    ```bash 
    git clone https://github.com/MolecularAI/REINVENT4.git
    ```
2. Configure transfer learning (aka fine-tuning)
   1. Adjust `configs/toml/transfer_learning.toml` following provided `configs/toml/README.md` instructions, 
   2. Set input model file for Mol2Mol generator as provided by authors `priors/reinvent.prior`.
   3. Set the following parameters:
     ```ini
     num_epochs = 1000
     save_every_n_epochs = 10
     batch_size = [adjust appropriately to reduce training time]
     ```
3) Train the model:
     Run the training using the modified configuration file. It takes approximetely 72 hours to train a model on ~750 samples with that setup. 
     
     Once trained, fine-tuned model can be used for downstream evaluation and benchmarking tasks.
--- -->
