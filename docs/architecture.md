# Architecture

## Project Overview

`hedgehog` is a Python-first molecular evaluation pipeline with a Typer CLI and an Ink-based TUI.
It processes generated molecules through staged filtering, docking, and reporting workflows.
The repository also contains helper modules and a docs website.

## Domain Map

| Domain | Description | Key Paths |
| --- | --- | --- |
| CLI and orchestration | CLI commands, argument handling, run lifecycle, setup commands | `src/hedgehog/main.py`, `src/hedgehog/pipeline.py` |
| Pipeline stages | Stage implementations for preparation, descriptors, filtering, synthesis, docking | `src/hedgehog/stages/**` |
| Reporting and metrics | HTML report generation and MolEval-based metrics | `src/hedgehog/reporting/**`, `src/hedgehog/vendor/moleval/**` |
| TUI + backend RPC | Node/Ink TUI and Python JSON-RPC backend for interactive runs | `tui/**`, `src/hedgehog/tui_backend/**` |
| Setup and external tools | Installation helpers for optional external dependencies | `src/hedgehog/setup/**`, `modules/**` |

## Layers

## Presentation Layer

- CLI entry points via Typer (`hedgehog`, `hedge`).
- TUI entry point and screen flow in Ink/TypeScript.

## Service Layer

- Pipeline orchestration (`calculate_metrics`) and stage sequencing.
- TUI backend JSON-RPC handlers for config, files, validation, history, and pipeline control.

## Persistence and Filesystem Layer

- YAML configuration files under `src/hedgehog/configs/`.
- Pipeline run artifacts under results folders (`input/`, `stages/`, `output/`, `configs/`).

## Integrations Layer

- RDKit/datamol/medchem for molecular processing.
- Optional docking and synthesis tooling (GNINA/SMINA/AiZynthFinder/worker setup).
- Vendored `moleval` metrics implementation.

## Key Interfaces

- `hedgehog.main:app`: CLI contract for run/report/info/version/setup/tui commands.
- `hedgehog.pipeline.calculate_metrics(...)`: pipeline orchestration contract used by CLI and backend.
- `hedgehog.tui_backend.server.JsonRpcServer`: JSON-RPC transport contract for the TUI.

## External Dependencies

| Dependency | Purpose | Surface |
| --- | --- | --- |
| RDKit, datamol, medchem | Molecule parsing and transformations | pipeline stages and input preparation |
| Typer, Rich | CLI UX and command parsing | `src/hedgehog/main.py` |
| Node + Ink | Interactive terminal UI | `tui/**` |
| Optional docking/synthesis tools | Extended scoring and retrosynthesis | setup commands and stage configs |

## Known Constraints

- Some workflows depend on heavy optional tools and system-specific binaries.
- TUI smoke execution requires a TTY; non-TTY environments should at least validate build success.
- Vendored third-party code exists under `src/hedgehog/vendor/` and should be changed only with strong justification.

Last reviewed: 2026-02-24.
