# Agent Dispatch

This file is a dispatcher for agent work. Detailed guidance lives under `docs/`.

## Working Agreements

- Write all file content and terminal output in English.
- Do not speculate about code you did not open.
- If a file is referenced, read it before making claims.
- Keep changes focused; do not add unrelated refactors.
- If behavior changes, update docs in the same PR.

## Context Graph

### Core Nodes

- [`docs/architecture.md`](docs/architecture.md): domain map, layers, boundaries.
- [`docs/conventions.md`](docs/conventions.md): naming, error/logging, tests, review.
- [`docs/quality.md`](docs/quality.md): A-F quality baseline and improvement plan.
- [`docs/tools/index.md`](docs/tools/index.md): helper scripts and when to use them.
- [`docs/tools/agent-mode.md`](docs/tools/agent-mode.md): agent runtime workflow.

### Runbooks

- [`docs/runbooks/ci-verification.md`](docs/runbooks/ci-verification.md): CI expectations.
- [`docs/runbooks/server-final-verification.md`](docs/runbooks/server-final-verification.md): end-to-end validation flow.

## Task Routing

- Architecture or cross-cutting change: read `docs/architecture.md` first.
- Coding style, test, or review questions: read `docs/conventions.md`.
- Risky or fragile area changes: read `docs/quality.md` and apply mitigation steps.
- New/changed scripts: update `docs/tools/index.md` and add tool docs.

## Required Smoke Checks

### CLI (when touching `src/hedgehog/**`, `pyproject.toml`, or CLI behavior)

```bash
uv run hedgehog --help
uv run hedge --help
uv run hedgehog run --help
uv run hedgehog version
```

### TUI (when touching `tui/**`, `src/hedgehog/tui_backend/**`, or `hedgehog tui`)

```bash
cd tui && npm ci
cd tui && npm run build
timeout 3 node tui/dist/index.js
```

In non-TTY environments, treat successful build as the minimum signal.

## Maintenance Rules

- Keep this file under 100 lines.
- Keep details in `docs/*`; keep this file pointer-only.
- Validate links and remove stale pointers after major refactors.
