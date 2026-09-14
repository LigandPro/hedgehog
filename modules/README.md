# Modules workspace

This folder is where Hedgehog keeps **optional external tools and assets**
(AiZynthFinder data, Lilly binaries, SYBA checkout, Sync checkpoint, FSScore/GASA
checkouts, docking helpers, …).

You usually don’t edit code here by hand. Prefer:

```bash
uv run hedgehog setup <thing>
# or a full pipeline with auto-bootstrap
uv run hedgehog --auto-install
```

---

## What’s typically here

| Path | What it is | How it usually arrives |
| --- | --- | --- |
| `aizynthfinder/` | Retrosynthesis engine + public data | `uv run hedgehog setup aizynthfinder` |
| `lilly_medchem_rules/bin/` | Lilly CLI binaries for demerit scoring | Vendored with the repo |
| `syba/` | SYBA synthesizability scorer sources | Vendored / env install |
| `sync/` | SYNC model checkpoint | `setup sync` or `sync_auto_install` |
| `fsscore/` / `gasa/` | Optional score checkouts | `--auto-install` when enabled in synthesis config |
| `mce18.py` | MCE-18 descriptor helper | Vendored from [Tong-Du/MCE-18](https://github.com/Tong-Du/MCE-18) |
| `docking/` / `dockingTools/` / `matcha_remote/` / `sigmatcha/` | Docking-related assets | Stage/setup dependent |

Exact contents can grow as you enable optional scorers. That’s normal.

---

## AiZynthFinder (most common setup)

**Recommended**

```bash
uv run hedgehog setup aizynthfinder
```

That installs the optional retrosynthesis extra, downloads public data, and
wires logging.

**Legacy script** (same idea, older entrypoint):

```bash
cd modules
./install_aizynthfinder.sh
```

**Manual sketch** (only if you need full control):

```bash
uv sync --extra retrosynthesis
mkdir -p modules/aizynthfinder/public modules/aizynthfinder/aizynthfinder/data
uv run python -m aizynthfinder.tools.download_public_data modules/aizynthfinder/public
cp src/hedgehog/synthesis/logging.yml modules/aizynthfinder/aizynthfinder/data/logging.yml
```

Then in `src/hedgehog/configs/config_synthesis.yml`:

```yaml
run: true
run_retrosynthesis: true
n_jobs: 32
# aizynthfinder_* knobs live in the same file
```

---

## Other optional pieces (short recipes)

```bash
uv run hedgehog setup sync
uv run hedgehog setup fsscore
uv run hedgehog setup gasa
```

For synthesis optional scorers (`sync` / `fsscore` / `gasa` / `nonpher`), put them
in `enabled_scores` and prefer `--auto-install` so paths don’t have to be typed
by hand. Override paths only when you bring your own env.

---

## Troubleshooting

| Problem | Likely fix |
| --- | --- |
| Synthesis stage skips retrosynthesis | `run_retrosynthesis: true` + AiZynthFinder setup done |
| FSScore/GASA are all `NaN` | Enable scorer in `enabled_scores` and run with `--auto-install` |
| Lilly filter fails to start | Check `modules/lilly_medchem_rules/bin` is on PATH (Hedgehog usually sets this) |
| “Where do I put custom models?” | Prefer config/ENV overrides; otherwise `$HEDGEHOG_OPTIONAL_ENV_ROOT/...` |

Main project setup: see the root [README.md](../README.md).
Deep stage docs: [https://hedgehog.ligandpro.ru](https://hedgehog.ligandpro.ru).
