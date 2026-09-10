# Hedgehog TUI

Terminal UI for configuring and running the Hedgehog pipeline without hand-editing
every YAML file.

If the CLI feels opaque, start here — then open the stage configs only when you
need a precise knob.

---

## What you can do

| Screen / flow | Purpose |
| --- | --- |
| Config editors | Edit mol prep, descriptors, structural filters, synthesis, docking, … |
| Wizard | Guided setup for a new run |
| Pipeline runner | Launch a run and watch progress |
| File browser | Pick molecule / receptor inputs |
| History | Reopen past runs |

Backend lives in `src/hedgehog/tui_backend/` (JSON-RPC over stdin/stdout).

---

## Requirements

- Node.js **18+** and npm  
- Python **3.10+** with a Hedgehog checkout (`uv sync` from repo root)

---

## Install & run

From the repo:

```bash
cd tui
npm install
```

**Dev (auto-reload)**

```bash
npm run dev
```

**Production build**

```bash
npm run build
npm start
```

**From project root**

```bash
uv run hedgehog tui
# or
cd tui && npm run tui
```

**Smoke check** (build + can the TUI talk to Python?)

```bash
uv run python scripts/check_pipeline.py --mode quick
```

---

## Everyday recipes

```bash
# 1) Open TUI, pick moses_1000 / your CSV, enable stages you care about
uv run hedgehog tui

# 2) Prefer a safe first run: mol prep + descriptors + structural filters
#    (leave synthesis/docking off until optional tools are set up)

# 3) If a stage looks “dead”, check its `run: true/false` in the matching YAML
#    under src/hedgehog/configs/
```

Configs the TUI edits are the same files the CLI uses
(`config_mol_prep.yml`, `config_structFilters.yml`, …).

For Common Alerts SMARTS / rulesets, see  
`src/hedgehog/struct_filters/data/README.md`.

---

## Key bindings

| Key | Action |
| --- | --- |
| `↑` / `↓` | Move in lists |
| `Enter` | Select / confirm |
| `→` / `e` | Edit path in file browser |
| `Space` | Toggle option / quick search in file browser |
| `Esc` / `←` | Back (most screens) |
| `Esc` or `q` | Quit on Welcome |
| `Ctrl+C` | Quit anywhere |
| `Ctrl+F` | Search / filter |
| `/` | Command palette |
| `?` | Help overlay |

---

## Project map (short)

```text
tui/
├── src/
│   ├── screens/      # Welcome, configs, runner, history, wizard
│   ├── components/   # Shared UI bits
│   ├── services/     # python-bridge / RPC
│   ├── store/        # Zustand state
│   └── types/
├── bin/              # CLI entry
└── package.json
```

---

## Troubleshooting

| Problem | Fix |
| --- | --- |
| TUI starts but configs don’t save | Run from a full repo checkout; configs are relative paths |
| Backend / RPC errors | `uv sync` at repo root, then `npm run tui` again |
| Node too old | Upgrade to Node 18+ |
| “It works in CLI, not in TUI” | Confirm both use the same `src/hedgehog/configs/*.yml` |

More detail: [docs TUI page](../docs/content/tui.mdx) and [public docs](https://hedgehog.ligandpro.ru).
