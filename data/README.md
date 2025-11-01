# Data Directory Documentation

This directory contains reference datasets and filter definitions used by the HEDGE framework for molecular generation evaluation and structural filtering.

## Files Overview

### 1. `common_alerts_collection.csv` (442KB, 2,459 rows)
**Purpose**: Comprehensive collection of structural alerts for drug discovery

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

**Configuration**: Specific rule sets and descriptions can be included/excluded via the `include_rulesets` and `exclude_descriptions` sections in `config_structFilters.yml`

---

### 2. `KRAS_target_mols.csv` (54KB, 742 rows)
**Purpose**: Collection of molecules targeting the KRAS protein

**Usage**: Reference dataset for benchmarking molecular generation models on KRAS-targeting compounds

**Columns**:
- `SMILES`: Canonical SMILES representations of KRAS-targeting molecules

**Notes**: 
- Contains 741 unique KRAS inhibitors and modulators
- Can be used as a target dataset for conditional generation tasks
- Useful for evaluating model performance on specific protein targets

---

### 3. `mcf.csv` (739B, 24 rows) 
**Purpose**: Medicinal Chemistry Filters (MCF) definitions

**Usage**: Defines SMARTS patterns for common medicinal chemistry filters

**Columns**:
- `names`: Filter identifier (MCF1-MCF22)
- `smarts`: SMARTS pattern defining the molecular filter

**Filter Types Include**:
- MCF1: α,β-unsaturated nitriles
- MCF2: α,β-unsaturated sulfones  
- MCF3: α,β-unsaturated amides
- MCF4: Alkyl halides
- MCF5: Epoxides
- MCF6: Isocyanates
- MCF7: Aldehydes
- MCF8: Imines
- MCF9: Aziridines
- MCF10-22: Various other problematic substructures

---

### 4. `wehi_pains.csv` (203KB, 480 rows)
**Purpose**: WEHI (Walter and Eliza Hall Institute) PAINS filter collection

**Usage**: Specialized PAINS (Pan-Assay Interference Compounds) filters for identifying problematic compounds in screening

**Format**: Single column containing SMARTS patterns with embedded identifiers

**Content**: 
- Each row contains a SMARTS pattern for identifying PAINS
- Patterns are annotated with rule identifiers (e.g., `<regId=anil_di_alk_F(14)>`)
- Covers various classes of interfering compounds like:
  - Aniline derivatives
  - Hydrazones  
  - Quinones
  - Thio compounds
  - Reactive heterocycles
  - Michael acceptors

---

## Configuration
`KRAS_target_mols.csv` is main orient as target molecules. 

`wehi_pains.csv` and `mcf.csv` are used by default to evaluate models with `MolScore` module. 

The usage of `common_alerts_collection.csv` file is controlled through `configs/config_structFilters.yml`:

```yaml
# Enable/disable structural alerts filtering
calculate_common_alerts: True

# Path to alerts collection
alerts_data_path: "data/common_alerts_collection.csv"

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
    # ... etc
```

> **!CAUTION** \
> Please do **NOT** completely remove `wehi_pains.csv`, `mcf.csv` and `common_alerts_collection.csv` files, in sake of code runnig. 
