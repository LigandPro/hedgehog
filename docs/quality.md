# Quality

## Grading Rubric

- `A`: explicit contracts, strong tests, low-risk modifications.
- `B`: mostly clear behavior with minor ambiguity or integration risk.
- `C`: workable but with notable fragility, sparse docs, or heavy implicit assumptions.
- `D`: high ambiguity and frequent breakage risk.
- `F`: unsafe to change without deep tribal knowledge.

## Domain Grades

| Domain | Grade | Evidence | Owner | Target |
| --- | --- | --- | --- | --- |
| CLI and pipeline orchestration | B | Strong command surface and broad pytest coverage for core flows | Maintainers | A- |
| Pipeline stage implementations | B- | Many stage-specific tests, but external tool coupling increases variability | Maintainers | B+ |
| TUI + Python backend | B- | Dedicated tests and build/smoke script exist; runtime depends on TTY and Node/Python bridge | Maintainers | B+ |
| Setup and external integrations | C | Environment-dependent installs and optional binary dependencies increase failure modes | Maintainers | B |
| Vendored metrics code | C | Useful but partly third-party and less ergonomic to modify safely | Maintainers | B- |

## Current Risks

- Optional external dependencies can fail due to environment drift.
- TUI runtime diagnostics rely on PTY behavior and cross-process communication.
- Vendored code updates can introduce regressions if changed without targeted validation.

## Improvement Plan

1. Add explicit integration checklists per setup path and keep them versioned in docs.
2. Strengthen failure-mode tests around setup/install and missing toolchain scenarios.
3. Expand TUI backend contract tests for JSON-RPC error shapes and lifecycle events.
4. Track vendor update policy and add focused regression smoke commands.

Last reviewed: 2026-02-24.
