@AGENTS.md

# HEDGEHOG

Python CLI pipeline for evaluating generative models in molecular design (drug discovery).
7 stages: molPrep → descriptors → structFilters → synthesis → docking → dockingFilters → finalDescriptors.

## Context Graph

Load these as needed — do NOT read all at once:

- `agent-docs/architecture.md` — Pipeline architecture, stage contracts, data flow, key classes
- `agent-docs/conventions.md` — Code style (Ruff), naming, imports, testing patterns
- `agent-docs/quality.md` — CI pipeline, local checks, test architecture
- `agent-docs/tools.md` — CLI commands, dev tools, SLURM scripts, config templates
- `agent-docs/domain-glossary.md` — Cheminformatics terminology, tools, internal terms

## Critical Paths

| What | Path |
|------|------|
| CLI entry | `src/hedgehog/main.py` |
| Pipeline orchestration | `src/hedgehog/pipeline.py` |
| Stage modules | `src/hedgehog/stages/{stageName}/main.py` |
| Master config | `src/hedgehog/configs/config.yml` + per-stage YAMLs |
| TUI backend | `src/hedgehog/tui_backend/server.py` |
| User docs | `docs/content/` (Nextra MDX) |

## Common Commands

```bash
uv run hedgehog run                  # run full pipeline
uv run hedgehog run --stage <name>   # run single stage
uv run ruff check .                  # lint
uv run ruff format .                 # format
uv run pytest                        # tests
cd docs && pnpm dev                  # docs dev server
cd tui && npm run build              # build TUI
```
