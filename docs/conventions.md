# Conventions

## Naming

- Python modules, files, and functions: `snake_case`.
- Python classes: `PascalCase`.
- Constants: `UPPER_SNAKE_CASE`.
- TypeScript files in `tui/src`: `PascalCase.tsx` for components, `camelCase.ts` for hooks/utils.
- CLI options and flags: kebab-case, implemented through Typer options.

## Error Handling

- Raise explicit exceptions instead of silent fallbacks in internal logic.
- At CLI and RPC boundaries, convert exceptions into actionable user-facing messages.
- Do not swallow exceptions without logging context.
- Prefer fail-fast behavior for invalid configuration and missing required files.

## Logging

- Keep logs concise and in English.
- Use structured, contextual messages (`what failed`, `where`, `next action`).
- Respect plain-output mode (`HEDGEHOG_PLAIN_OUTPUT=1`) for machine-friendly output.
- Avoid noisy debug logs in default paths unless diagnosing failures.

## Testing

- Test framework: `pytest` with tests under `tests/test_*.py`.
- Keep test names behavior-oriented (`test_<expected_behavior>`).
- For CLI-affecting changes, run required smoke checks:
  - `uv run hedgehog --help`
  - `uv run hedge --help`
  - `uv run hedgehog run --help`
  - `uv run hedgehog version`
- For TUI/backend changes, run:
  - `cd tui && npm ci`
  - `cd tui && npm run build`
  - `timeout 3 node tui/dist/index.js` (TTY preferred)

## Code Review Expectations

- Include tests or smoke evidence for behavior changes.
- Update docs when command behavior, defaults, or operational flow changes.
- Keep diffs focused; avoid opportunistic refactors in unrelated domains.
- For risky integration changes, include rollback or mitigation notes in PR description.
