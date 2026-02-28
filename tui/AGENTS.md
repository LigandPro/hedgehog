# AGENTS (TUI)

This file applies to `tui/**` and augments the repository-level [AGENTS.md](../AGENTS.md).
If rules conflict, follow repository-level rules first.

- ALWAYS USE PARALLEL TOOLS WHEN APPLICABLE.
- Prefer automation: execute requested actions directly unless blocked by missing info or safety/irreversibility.
- The default branch in this repo is `main`.
- Keep changes focused. Do not perform unrelated refactors.

## Style Guide

### General Principles

- Keep things in one function unless composable or reusable.
- Avoid `try`/`catch` where possible.
- Avoid using the `any` type.
- Prefer single word variable names where possible.
- Use Node APIs when possible (`node:fs`, `node:path`, `node:os`, etc.).
- Rely on type inference when possible; avoid explicit annotations unless necessary for exports/contracts or clarity.
- Prefer functional array methods (`flatMap`, `filter`, `map`) over manual loops where readability stays good.

### Naming

Prefer short names for variables and functions. Use multiple words only when necessary for clarity.

```ts
// Good
const job = state.currentJob
function render(list: Row[]) {}

// Bad
const currentJobFromState = state.currentJob
function renderTableRowsForFooter(list: Row[]) {}
```

Reduce total variable count by inlining values used once.

```ts
// Good
const file = await readFile(path.join(dir, "state.json"), "utf8")

// Bad
const statePath = path.join(dir, "state.json")
const file = await readFile(statePath, "utf8")
```

### Destructuring

Avoid unnecessary destructuring. Use dot notation to preserve context.

```ts
// Good
state.screen
state.searchActive

// Bad
const { screen, searchActive } = state
```

### Variables

Prefer `const` over `let`. Use ternaries or early returns instead of reassignment.

```ts
// Good
const mode = running ? "running" : "idle"

// Bad
let mode
if (running) mode = "running"
else mode = "idle"
```

### Control Flow

Avoid `else` statements where early returns are clearer.

```ts
// Good
function title(screen: string) {
  if (screen === "welcome") return "Home"
  return "Other"
}

// Bad
function title(screen: string) {
  if (screen === "welcome") return "Home"
  else return "Other"
}
```

## Testing

- Avoid mocks as much as possible.
- Test actual behavior; do not duplicate implementation logic into tests.
- For `tui/**` changes, run checks from the `tui` directory context.

## Required Verification (for `tui/**` changes)

```bash
cd tui && npm ci
cd tui && npm run build
timeout 3 node tui/dist/index.js
```

Notes:
- Ink runtime requires a TTY for full interaction.
- In non-TTY environments, treat successful build as minimum required signal.
- On macOS, if `timeout` is unavailable, use an equivalent timed runner.

## TUI UX Rules

- Shortcut hints must be truthful: only show actions that exist in current context.
- Keep shortcut style consistent across screens (`key` segment + `label` segment).
- Validate narrow terminal behavior:
  - no broken ASCII/logo rendering
  - graceful truncation of long hint lines
  - scrolling for long selectable lists

## Agent Roles (for multi-step TUI work)

- `impl-ui`: components/screens/rendering changes in `tui/src/components/**`, `tui/src/screens/**`
- `impl-state`: store/hooks/keybindings logic in `tui/src/store/**`, `tui/src/hooks/**`, app wiring
- `rv-build`: run checks and report build/runtime failures with command output

Rules:
- Lead coordinates and merges results.
- No phase is done until verification gates pass.
