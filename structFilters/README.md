## Structural filters module

Run medicinal chemistry and structural alert filters on generated molecules, before and/or after descriptors. This stage aggregates several rule sets and algorithmic filters and outputs per‑filter summaries, extended details, and pass lists for downstream evaluation.

The stage reads `configs/config_structFilters.yml` 

### What it does
- Filters:
  - Common structural alerts using `data/common_alerts_collection.csv` with selectable `include_rulesets` and `exclude_descriptions` (e.g., Dundee, BMS, Glaxo, PAINS, etc.).
  - Molecular Graph severity filter (`molgraph_stats`, severities 1–11).
  - Molecular complexity filters set (`molcomplexity`).
  - NIBR structural filters (`NIBR`).
  - Bredt bridgehead filter (`bredt`).
  - Lilly demerits (`lilly`).
- Runs either on original sampled molecules (`beforeDescriptors`) or on `Descriptors/passDescriptorsSMILES.csv` (`Descriptors`).
- Produces metrics, extended per‑molecule details, pass and fail files, and plots.