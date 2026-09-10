# Structural filters data

This folder is the **cheat sheet + source of truth** for Common Alerts and protecting-group SMARTS.

If something looks wrong in Stage 3, start here: check which file you need, copy the exact SMARTS/ruleset name, then paste it into `config_structFilters.yml` (or a strict/exploration preset).

---

## What’s in this folder


| File                           | Loaded by pipeline? | Use it for                             |
| ------------------------------ | ------------------- | -------------------------------------- |
| `common_alerts_collection.csv` | **Yes**             | Running Common Alerts                  |
| `common_alerts_catalog.json`   | No                  | Browsing rules / copying SMARTS safely |
| `protecting_groups.csv`        | **Yes**             | Protecting-groups filter               |


Configs that point here:

- `src/hedgehog/configs/config_structFilters.yml`
- `…_strict.yml` / `…_exploration.yml`

---



## 1. `common_alerts_collection.csv` (the real catalog)

~2 458 rules across **23** rulesets (PAINS, Glaxo, Dundee, Toxicophore, …).

Hedgehog loads this when:

```yaml
calculate_common_alerts: true
alerts_data_path: src/hedgehog/struct_filters/data/common_alerts_collection.csv
```



### Columns (with real examples)


| Column                | Meaning                               | Example                    |
| --------------------- | ------------------------------------- | -------------------------- |
| `rule_id`             | Unique ID                             | `340`                      |
| `rule_set`            | Numeric family id                     | `15`                       |
| `rule_set_name`       | Family name (use this in YAML)        | `PAINS`                    |
| `description`         | Alert name (can be empty or numeric!) | `ene_six_het_A(483)`       |
| `smarts`              | Pattern that actually matches         | `[#6]-1(-[#6](~[!#6&!#1]…` |
| `priority`            | Rule priority                         | `5`                        |
| `mincount`            | Min hits required                     | `1`                        |
| `source`              | Origin                                | `ChEMBL` or `Litterature`  |
| `catalog_description` | Short family blurb                    | (text)                     |


**Tricky rows you’ll meet**


| Situation                        | Example                                           | What to do                                              |
| -------------------------------- | ------------------------------------------------- | ------------------------------------------------------- |
| Normal named alert               | PAINS / `ene_six_het_A(483)`                      | Fine to talk about by description                       |
| Empty description                | DNABinder / `description=""` / SMARTS `[NH][NH2]` | Exclude via **SMARTS**, not name                        |
| Numeric “name”                   | LINT / `description=26` / SMARTS `N-C(=S)-N`      | Prefer **SMARTS**; the number is not a chemical name    |
| Same description, several SMARTS | happens inside some families                      | Don’t key only by description — use SMARTS or `rule_id` |




### Mini example rows

```text
rule_id,rule_set_name,description,smarts,source
340,PAINS,ene_six_het_A(483),[#6]-1(-[#6](~[!#6&!#1]~[#6]-[!#6&!#1]-[#6]-1=[!#6&!#1])~[!#6&!#1])=[#6;!R],ChEMBL
1388,DNABinder,,[NH][NH2],Litterature
1219,LINT,26,N-C(=S)-N,Litterature
```

---



## 2. How to configure Common Alerts without shooting yourself in the foot

Think in **two layers**:

1. **Calculate** — which alerts are scored / written to diagnostics
2. **Filter** — which of those become a hard reject (`filter_common_alerts`)



### Layer A — what gets calculated

```yaml
calculate_common_alerts: true
alerts_data_path: src/hedgehog/struct_filters/data/common_alerts_collection.csv

# Which families to load from the CSV:
include_rulesets: all          # or null/[] for none, or an explicit list:
# include_rulesets:
#   - PAINS
#   - Glaxo
#   - Dundee

# Drop exact SMARTS strings (copy-paste from CSV/catalog, character-perfect):
exclude_smarts:
  - "[NH][NH2]"
  - "N-C(=S)-N"
```



### Layer B — what becomes a hard gate

```yaml
filter_common_alerts: true

# When filtering is on, optionally narrow further:
common_alerts_filter_include_rulesets:   # empty = every calculated ruleset
  - PAINS
common_alerts_filter_exclude_rulesets: []  # remove whole families from the gate
```

**Mental model**

- `include_rulesets` / `exclude_smarts` → shape the **calculation**  
- `filter_common_alerts` + `common_alerts_filter_*` → shape the **reject gate**  
- Diagnostics (e.g. `hits_long.csv` under the common_alerts outputs) show *what hit* — they don’t configure exclusions by themselves



### Typical recipes

**“Run everything for diagnostics, gate only on PAINS”**

```yaml
calculate_common_alerts: true
include_rulesets: all
exclude_smarts: []
filter_common_alerts: true
common_alerts_filter_include_rulesets:
  - PAINS
```

**“Keep PAINS+Glaxo, but silence two noisy SMARTS”**

```yaml
calculate_common_alerts: true
include_rulesets:
  - PAINS
  - Glaxo
exclude_smarts:
  - "cN=[N+]=[N-]"
  - "[$(N#[N+]-[N-]),$([N-]=[N+]=N)]"
filter_common_alerts: true
common_alerts_filter_include_rulesets: []   # gate on whatever you calculated
```

**“Calculate alerts but never reject on them”**

```yaml
calculate_common_alerts: true
include_rulesets: all
filter_common_alerts: false
```

---



## 3. `common_alerts_catalog.json` (for humans)

Same rules as the CSV, nested by ruleset:

```json
{
  "_meta": {
    "purpose": "...",
    "pitfalls": ["description is NOT unique", "..."],
    "how_to_use": ["Prefer rule_id or smarts for lookups", "..."]
  },
  "rulesets": {
    "PAINS": {
      "n_rules": 481,
      "rules": [
        {
          "rule_id": "340",
          "description": "ene_six_het_A(483)",
          "smarts": "[#6]-1(-[#6](~[!#6&!#1]~[#6]-[!#6&!#1]-[#6]-1=[!#6&!#1])~[!#6&!#1])=[#6;!R]"
        }
      ]
    }
  }
}
```

**Workflow**

1. Open the catalog (or CSV).
2. Find the ruleset → rule.
3. Copy the `smarts` field into `exclude_smarts`.
4. Do **not** build a dict keyed only by `description` — collisions and empty names will bite you.

---



## 4. `protecting_groups.csv`

Used when `calculate_protecting_groups: true`. Loaded from this folder automatically.


| Column      | Example                            |
| ----------- | ---------------------------------- |
| `name`      | `fmoc`                             |
| `smiles`    | `O=COCC1C2=C(C3=C1C=CC=C3)C=CC=C2` |
| `smarts`    | `[#8]=[#6]-[#8]-[#6]-[#6]1-…`      |
| `group`     | `protecting_groups`                |
| `hierarchy` | `hedgehog.protecting_groups`       |


Hard reject for this filter is controlled separately with `filter_protecting_groups`.

---



## 5. Quick troubleshooting


| Problem                                     | Likely cause                             | Fix                                                    |
| ------------------------------------------- | ---------------------------------------- | ------------------------------------------------------ |
| Alert family never appears                  | Not in `include_rulesets`                | Add the family or use `all`                            |
| Molecule still hits an alert you “excluded” | Excluded by description / typo in SMARTS | Paste **exact** `smarts` into `exclude_smarts`         |
| Gate feels softer/harder than plots         | `calculate_*` vs `filter_*` mismatch     | Plots/diagnostics ≠ hard reject                        |
| Empty or weird alert names                  | DNABinder / LINT-style rows              | Use SMARTS / `rule_id` from catalog                    |
| Wrong file edited                           | Editing the JSON catalog                 | Pipeline reads the **CSV** (and protecting_groups CSV) |


---



## 6. Related knobs outside this folder

These aren’t in `data/`, but people confuse them with alerts:

- `filter_NIBR` / `nibr_max_severity`
- `filter_lilly` / `lilly_demerit_cutoff` (soft demerit style)
- `filter_molgraph_stats` / `molgraph_max_severity`
- `filter_undefined_stereo_center` / `stereo_max_undefined`
- output toggles: `write_per_filter_outputs`, `write_structural_liability_profile`, `generate_plots`, `generate_failure_analysis`

Full stage policy lives in the YAML configs under `src/hedgehog/configs/`.