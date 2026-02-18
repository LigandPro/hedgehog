## Agent Rules (Repository-Wide)

These instructions apply to **all** changes in this repository.

### Language

- **All content written to files must be in English.**
- Keep console output/logging in English.

### No Guessing

- Do not speculate about code you have not opened.
- If a user references a specific file, read it before answering.

### Required Checks When Adding/Changing Code

When you add or change code, you must run the relevant smoke checks and keep documentation consistent.

#### CLI (Python / Typer)

Run these commands whenever you touch:
- `src/hedgehog/**`
- `pyproject.toml`
- `docs/content/cli.mdx`
- anything affecting CLI behavior

Commands:

```bash
uv run hedgehog --help
uv run hedge --help
uv run hedgehog run --help
uv run hedgehog version
```

#### TUI (Node / Ink) + Python Backend

Run these commands whenever you touch:
- `tui/**`
- `src/hedgehog/tui_backend/**`
- the `hedgehog tui` command implementation

Commands:

```bash
cd tui && npm ci
cd tui && npm run build
timeout 3 node tui/dist/index.js
```

Notes:
- `node tui/dist/index.js` requires a TTY (Ink). Run it from an interactive terminal.
- In non-TTY environments, treat build success as the minimum signal.

### Documentation Sync (Source of Truth)

- The source of truth for CLI behavior is `uv run hedgehog ... --help`.
- If you add/change CLI commands, options, defaults, or error behavior, update `docs/content/cli.mdx` accordingly.
- If README examples become inaccurate, update `README.md` as well.

### CI Verification

- Always check CI status before considering a task complete.
- For PR work, verify the latest remote CI checks for the current branch/PR.
- If CI is failing, treat it as an active issue and fix it or clearly report the blocker.
