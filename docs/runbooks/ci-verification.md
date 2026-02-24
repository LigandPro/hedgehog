# CI Verification Runbook

## Policy

- Always check CI status before considering work complete.
- For PR branches, inspect the latest remote checks for the current branch/PR.
- Treat failing checks as active issues to fix or clearly report with blockers.

## Local Pre-CI Baseline

```bash
uv run python scripts/check_pipeline.py --mode ci
```

This runs CLI smoke, TUI build, and TUI smoke checks in one flow.

## Behavior Changes and Docs Sync

- CLI behavior source of truth is `uv run hedgehog ... --help`.
- When CLI flags/defaults/help output changes, update:
  - `docs/content/cli.mdx`
  - `README.md` examples if affected
