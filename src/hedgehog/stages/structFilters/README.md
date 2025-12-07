## Structural filters module

Run medicinal chemistry and structural alert filters on generated molecules, before and/or after descriptors. 

The stage reads [config_structFilters](/src/hedgehog/configs/config_structFilters.yml) file.

Structural Filters:
  - Common structural alerts using `data/common_alerts_collection.csv` with selectable `include_rulesets` and `exclude_descriptions` (e.g., Dundee, BMS, Glaxo, PAINS, etc.).
  - Molecular Graph severity filter (severities 1–11).
  - Molecular complexity filters set.
  - NIBR structural filters, [NIBR](https://iwatobipen.github.io/is_life_worth_living-/jupyter/2021/11/19/NIBR_SF.html).
  - Bredt bridgehead filter, [Bredt rule](https://en.wikipedia.org/wiki/Bredt%27s_rule).
  - Lilly demerits, [Eli Lilly Medchem Rules](https://pubs.acs.org/doi/10.1021/jm301008n).


**Run structural filters stage only:**
```bash
uv run hedgehog run --stage struct_filters
```
Alternatively, using the short-name alias:
```bash
uv run hedge run --stage struct_filters
```
