# AGENTS Dispatch Migration

Status: Implemented
Date: 2026-02-24

## Objective

Migrate root `AGENTS.md` from instruction monolith to thin dispatcher format and move details into dedicated docs nodes.

## Delivered

- Added:
  - `docs/architecture.md`
  - `docs/conventions.md`
  - `docs/quality.md`
  - `docs/tools/index.md`
  - `docs/tools/agent-mode.md`
  - `docs/runbooks/ci-verification.md`
  - `docs/runbooks/server-final-verification.md`
- Rewrote root `AGENTS.md` as pointer-first table of contents.

## Follow-ups

- Keep `docs/quality.md` grades updated as tooling/tests evolve.
- Add per-tool docs under `docs/tools/` when new helper scripts are introduced.
