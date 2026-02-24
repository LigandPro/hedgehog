# Agent Mode

This repository does not currently expose a single dedicated `--agent-mode` flag.
Use the workflow below to achieve low-noise, fail-fast agent runs.

## Runtime Flags and Defaults

- `HEDGEHOG_PLAIN_OUTPUT=1`: disables rich formatting for machine-friendly output.
- `uv run hedgehog run --stage <stage>`: narrows execution scope for targeted changes.
- `uv run hedgehog run --auto-install`: allows optional dependency bootstrap during runs.
- `uv run hedgehog run --progress`: enables live progress events for long stage runs.

For TUI smoke/diagnostics (`scripts/check_pipeline.py`):

- `HEDGEHOG_TUI_LOG=<path>`: captures TUI lifecycle logs.
- `HEDGEHOG_TUI_DEBUG=1`: enables extra TUI debug traces.

## Security Guardrails

- Keep any relaxed/diagnostic runtime behavior in local development only.
- Do not commit credentials, tokens, or host-specific secrets in scripts or docs.
- Preserve normal validation boundaries in production-like runs.

## Verification Commands

```bash
uv run hedgehog --help
uv run hedge --help
uv run hedgehog run --help
uv run hedgehog version
cd tui && npm ci
cd tui && npm run build
timeout 3 node tui/dist/index.js
```

## CI Integration Notes

- CI should execute the same smoke contract used locally.
- Prefer a single orchestrator command where possible:
  - `uv run python scripts/check_pipeline.py --mode ci`
- If CI fails, treat it as active work, not a post-merge cleanup item.
