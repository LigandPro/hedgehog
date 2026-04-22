# Data Directory

This directory contains structural alert definitions used by the **🦔 HEDGEHOG** structural filters stage.


**`common_alerts_collection.csv`** (442KB, 2,459 rows)
**Purpose**: Collection of structural alerts for drug discovery

**Usage**: Used by the structural alerts filter (`calculate_common_alerts: True`) in `config_structFilters.yml`

**Columns**:
- `rule_id`: Unique identifier for each rule
- `rule_set`: Numeric identifier for the rule set
- `description`: Human-readable description of the structural alert
- `smarts`: SMARTS pattern defining the problematic substructure
- `rule_set_name`: Name of the rule set (e.g., "Glaxo", "Dundee", "BMS", "PAINS")
- `priority`: Priority level of the rule
- `mincount`: Minimum count threshold
- `source`: Data source (typically "ChEMBL")
- `catalog_description`: Description of the catalog/source

**Rule Sets Included**:
- **Glaxo**: Glaxo Wellcome Hard filters for reactive/unstable compounds
- **Dundee**: University of Dundee NTD Screening Library filters
- **BMS**: Bristol-Myers Squibb HTS Deck filters  
- **PAINS**: Pan-Assay Interference Compounds
- **Inpharmatica**: Commercial filter collection
- **And others**: MLSMR, LINT, GST-Hitters, etc.

**Configuration**: Specific rule sets and descriptions can be included/excluded via the `include_rulesets` and `exclude_descriptions` sections in `config_structFilters.yml`. Additional ready-to-use profiles are available in `config_structFilters_strict.yml`, `config_structFilters_balanced.yml`, and `config_structFilters_exploration.yml`.

```yaml
# Enable/disable structural alerts filtering
calculate_common_alerts: True

# Path to alerts collection
alerts_data_path: "src/hedgehog/struct_filters/data/common_alerts_collection.csv"

# Select which rule sets to include
include_rulesets:
  - "Dundee"
  - "BMS" 
  - "PAINS"
  # ... etc

# Exclude specific descriptions within rule sets
exclude_descriptions:
  "Dundee":
    - "Aliphatic long chain"
    - "isolated alkene"
    - "triple bond"
  "Inpharmatica":
    - "Filter82_pyridinium"
  "LD50-Oral":
    - "phenylpiperazine"
  "SureChEMBL":
    - "aminothiazole"
  "LINT":
    - "high halogen content (>3)"
  "HIS-Hitters":
    - "Picolylamines_A"
```
