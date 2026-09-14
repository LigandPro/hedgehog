# Synthesis data

Small helper files for the synthesis / SYNC scoring path.

You normally **don’t point YAML at these by hand** — the code knows the default
locations. This README is so you know what they are if you open the folder.

---

## What’s here

| File | Used by pipeline? | Role |
| --- | --- | --- |
| `sync_smiles_vocab_with_aro.txt` | **Yes** (SYNC) | Atom / token vocabulary for SYNC featurization |
| `common_reagents.smi` | Reference list | Conservative terminal reagents (notes for humans / experiments) |

---

## 1. `sync_smiles_vocab_with_aro.txt`

Plain text, one token per line (`H`, `C`, `Cpi`, `N`, …).

SYNC loads it from this folder when you enable the `sync` scorer:

```yaml
# in config_synthesis.yml
enabled_scores:
  - sa
  - syba
  - rascore
  - sync          # <-- needs the vocab + model checkpoint

sync_auto_install: true
sync_device: cpu    # or cuda:0
```

**If SYNC fails with a vocab / tokenization error**, check that this file is
present and hasn’t been truncated. Don’t reorder tokens casually — the model
expects this vocabulary.

Model checkpoint itself lives under `modules/sync/` (see `modules/README.md`).

---

## 2. `common_reagents.smi`

Short SMILES list with comments, e.g.:

```text
# Conservative terminal reagents commonly emitted as reactants by USPTO templates.
O
N
F
```

Handy as a **reference** for common building-block-like terminals. It is **not**
a drop-in replacement for a real vendor stock and is not the main AiZynthFinder
config.

---

## Related knobs (not in this folder)

| Goal | Where |
| --- | --- |
| Turn synthesis stage on/off | `config_synthesis.yml` → `run` |
| Scores vs retrosynthesis | `enabled_scores`, `run_retrosynthesis` |
| Optional scorers auto-setup | `uv run hedgehog --auto-install` / `setup sync|fsscore|gasa` |
| AiZynthFinder assets | `modules/aizynthfinder/` |

Full synthesis docs: [https://hedgehog.ligandpro.ru](https://hedgehog.ligandpro.ru)  
Modules workspace: `modules/README.md`.
