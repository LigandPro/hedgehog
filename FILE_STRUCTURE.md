# HEDGE File Structure Documentation

This document describes the improved file organization structure for HEDGE pipeline outputs.

## Table of Contents

- [Overview](#overview)
- [New Structure](#new-structure)
- [Migration from Legacy Structure](#migration-from-legacy-structure)
- [Detailed Breakdown](#detailed-breakdown)
- [File Naming Conventions](#file-naming-conventions)

## Overview

The new file structure provides:
- **Hierarchical organization**: Clear separation of stages, outputs, configs, and inputs
- **Numbered stages**: Easy-to-understand execution order (01, 02, 03, ...)
- **Consistent naming**: snake_case throughout, standardized file names
- **Logical grouping**: Related files organized in subdirectories
- **Backward compatibility**: Legacy structure still supported

## New Structure

```
results/
├── input/
│   └── sampled_molecules.csv              # Input molecules (sampled from config)
│
├── stages/                                # All pipeline stages
│   ├── 01_descriptors_initial/
│   │   ├── metrics/
│   │   │   ├── descriptors_all.csv        # All computed descriptors
│   │   │   └── skipped_molecules.csv      # Failed to parse
│   │   ├── filtered/
│   │   │   ├── filtered_molecules.csv     # Molecules that passed filters
│   │   │   ├── failed_molecules.csv       # Molecules that failed filters
│   │   │   ├── descriptors_passed.csv     # Detailed metrics for passed
│   │   │   ├── descriptors_failed.csv     # Detailed metrics for failed
│   │   │   └── pass_flags.csv             # Pass/fail flags per descriptor
│   │   └── plots/
│   │       └── descriptors_distribution.png  # Visualization
│   │
│   ├── 02_structural_filters_pre/         # Pre-descriptors filters (optional)
│   │   ├── common_alerts/
│   │   │   ├── metrics.csv
│   │   │   ├── extended.csv
│   │   │   └── filtered_molecules.csv
│   │   ├── bredt/
│   │   ├── NIBR/
│   │   ├── filtered_molecules.csv         # Combined passed molecules
│   │   ├── failed_molecules.csv           # Combined failed molecules
│   │   ├── molecule_counts_comparison.png
│   │   └── restriction_ratios_comparison.png
│   │
│   ├── 03_structural_filters_post/        # Post-descriptors filters
│   │   ├── common_alerts/
│   │   │   ├── metrics.csv
│   │   │   ├── extended.csv
│   │   │   ├── filtered_molecules.csv
│   │   │   └── CommonAlertsBreakdown/    # Detailed analysis
│   │   │       ├── filter_failures_plot.png
│   │   │       ├── all_filters_reasons_plot.png
│   │   │       ├── complete_reasons_breakdown.csv
│   │   │       └── ...
│   │   ├── bredt/
│   │   ├── NIBR/
│   │   ├── lilly/
│   │   ├── molgraph_stats/
│   │   ├── molcomplexity/
│   │   ├── filtered_molecules.csv         # Combined result
│   │   ├── failed_molecules.csv           # Failed molecules
│   │   ├── molecule_counts_comparison.png
│   │   └── restriction_ratios_comparison.png
│   │
│   ├── 04_synthesis/
│   │   ├── synthesis_scores.csv           # RAScore, SAScore, etc.
│   │   ├── synthesis_extended.csv         # With retrosynthesis results
│   │   ├── filtered_molecules.csv         # Synthesizable molecules
│   │   ├── input_smiles.smi               # Input for AiZynthFinder
│   │   └── retrosynthesis_results.json    # AiZynthFinder output
│   │
│   ├── 05_docking/
│   │   ├── ligands.csv                    # Prepared ligands
│   │   ├── smina/
│   │   │   ├── poses/                     # Individual docking poses
│   │   │   │   ├── mol_0_pose.pdbqt
│   │   │   │   └── ...
│   │   │   └── scores.csv                 # SMINA scores
│   │   └── gnina/
│   │       ├── output.sdf                 # GNINA output
│   │       └── scores.csv                 # GNINA scores
│   │
│   └── 06_descriptors_final/              # Re-calculate descriptors on final set
│       ├── metrics/
│       │   └── descriptors_all.csv
│       ├── filtered/
│       │   └── filtered_molecules.csv
│       └── plots/
│           └── descriptors_distribution.png
│
├── output/                                # Final pipeline results
│   ├── final_molecules.csv                # Final filtered molecules
│   └── pipeline_summary.json              # Summary statistics
│
├── configs/                               # Configuration snapshot
│   ├── master_config_resolved.yml         # Resolved configuration
│   ├── config_descriptors.yml
│   ├── config_structFilters.yml
│   ├── config_synthesis.yml
│   └── config_docking.yml
│
└── logs/                                  # Pipeline logs
    └── pipeline_YYYYMMDD_HHMMSS.log
```

## Migration from Legacy Structure

### Legacy Structure (Old)

```
results/test/
├── beforeDescriptors_StructFilters/       # Confusing name
│   ├── Pains_metrics.csv                  # CamelCase
│   └── Pains_filteredMols.csv             # Mixed naming
├── Descriptors/                           # CamelCase
│   ├── perMoleculeDescriptors.csv         # camelCase
│   ├── passDescriptorsSMILES.csv          # Unclear
│   └── filteredMetricsDistribution.png
├── StructFilters/
│   ├── Brenk_metrics.csv
│   └── passStructFiltersSMILES.csv        # Repetitive pattern
├── Synthesis/
│   └── passSynthesisSMILES.csv
├── Docking/
│   ├── smina_results/                     # Mixed convention
│   └── gnina_results/
├── finalDescriptors/                      # camelCase (inconsistent!)
├── run_configs/
├── sampledMols.csv                        # camelCase
└── finalMolecules.csv
```

### Mapping: Old → New

| Legacy Path | New Path | Notes |
|-------------|----------|-------|
| `Descriptors/passDescriptorsSMILES.csv` | `stages/01_descriptors_initial/filtered/filtered_molecules.csv` | Standardized name |
| `Descriptors/perMoleculeDescriptors.csv` | `stages/01_descriptors_initial/metrics/descriptors_all.csv` | Clearer organization |
| `beforeDescriptors_StructFilters/` | `stages/02_structural_filters_pre/` | Numbered, clear name |
| `StructFilters/` | `stages/03_structural_filters_post/` | Numbered, explicit timing |
| `StructFilters/Pains_metrics.csv` | `stages/03_structural_filters_post/common_alerts/metrics.csv` | Grouped by filter |
| `Synthesis/passSynthesisSMILES.csv` | `stages/04_synthesis/filtered_molecules.csv` | Standardized |
| `Docking/smina_results/` | `stages/05_docking/smina/` | Cleaner hierarchy |
| `finalDescriptors/` | `stages/06_descriptors_final/` | Consistent naming |
| `finalMolecules.csv` | `output/final_molecules.csv` | Explicit output folder |
| `sampledMols.csv` | `input/sampled_molecules.csv` | Explicit input folder |
| `run_configs/` | `configs/` | Shorter, clearer |

## Detailed Breakdown

### Stage 01: Descriptors Initial

**Purpose**: Calculate 22 physicochemical descriptors and filter molecules

**Outputs**:
- `metrics/descriptors_all.csv`: All computed descriptors for all molecules
- `filtered/filtered_molecules.csv`: Molecules passing descriptor thresholds
- `filtered/failed_molecules.csv`: Molecules failing descriptor thresholds
- `filtered/descriptors_passed.csv`: Full metrics for passed molecules
- `filtered/descriptors_failed.csv`: Full metrics for failed molecules
- `filtered/pass_flags.csv`: Boolean flags showing which descriptors passed/failed
- `plots/descriptors_distribution.png`: Distribution visualization

### Stage 02: Structural Filters (Pre-descriptors)

**Purpose**: Apply structural filters before descriptor calculation (optional)

**Outputs**:
- `{filter_name}/metrics.csv`: Summary statistics
- `{filter_name}/extended.csv`: Detailed results with all columns
- `{filter_name}/filtered_molecules.csv`: Molecules passing this specific filter
- `filtered_molecules.csv`: Combined result from all filters
- `failed_molecules.csv`: Molecules failing any filter
- `molecule_counts_comparison.png`: Count comparison plot
- `restriction_ratios_comparison.png`: Heatmap of restriction ratios

**Available filters**: common_alerts, bredt, NIBR, lilly, molgraph_stats, molcomplexity

### Stage 03: Structural Filters (Post-descriptors)

**Purpose**: Apply comprehensive structural filters after descriptors

**Outputs**: Same as Stage 02, plus:
- `common_alerts/CommonAlertsBreakdown/`: Detailed failure analysis
  - `filter_failures_plot.png`: Failures per filter
  - `all_filters_reasons_plot.png`: Multi-panel reason breakdown
  - `complete_reasons_breakdown.csv`: Full reason catalog
  - `comprehensive_reasons_overview.png`: Top failure reasons
  - `filter_summary_table.csv`: Summary statistics

### Stage 04: Synthesis

**Purpose**: Evaluate synthetic accessibility

**Outputs**:
- `synthesis_scores.csv`: RAScore, SAScore, SCScore, SYBA
- `synthesis_extended.csv`: Includes retrosynthesis results
- `filtered_molecules.csv`: Synthesizable molecules
- `input_smiles.smi`: Input for AiZynthFinder
- `retrosynthesis_results.json`: AiZynthFinder raw output

### Stage 05: Docking

**Purpose**: Molecular docking with SMINA and/or GNINA

**Outputs**:
- `ligands.csv`: Prepared ligands metadata
- `smina/poses/`: Individual docking poses (.pdbqt files)
- `smina/scores.csv`: SMINA docking scores
- `gnina/output.sdf`: GNINA output (SDF format)
- `gnina/scores.csv`: GNINA CNN scores

### Stage 06: Descriptors Final

**Purpose**: Recalculate descriptors on final filtered set

**Outputs**: Same structure as Stage 01

## File Naming Conventions

### Principles

1. **snake_case everywhere**: `filtered_molecules.csv`, not `filteredMolecules.csv`
2. **Descriptive names**: `descriptors_all.csv`, not `perMoleculeDescriptors.csv`
3. **Standardized outputs**: All stages use `filtered_molecules.csv` for passed molecules
4. **Consistent prefixes**:
   - `descriptors_*` for descriptor files
   - `synthesis_*` for synthesis files
   - No mixed prefixes like `pass*SMILES`

### Standard File Names

| Purpose | File Name | Old Name |
|---------|-----------|----------|
| Molecules that passed | `filtered_molecules.csv` | `pass{Stage}SMILES.csv` |
| Molecules that failed | `failed_molecules.csv` | `fail{Stage}SMILES.csv` |
| All computed metrics | `descriptors_all.csv` | `perMoleculeDescriptors.csv` |
| Metrics for passed | `descriptors_passed.csv` | `passDescriptorsMetrics.csv` |
| Pass/fail flags | `pass_flags.csv` | `descriptorsPassFlags.csv` |
| Distribution plot | `descriptors_distribution.png` | `filteredMetricsDistribution.png` |
| Input molecules | `sampled_molecules.csv` | `sampledMols.csv` |
| Final output | `final_molecules.csv` | `finalMolecules.csv` |

## Implementation Notes

### Backward Compatibility

The code maintains full backward compatibility:
- Legacy paths are checked as fallbacks
- Old file names are recognized
- Gradual migration supported

### Configuration

No configuration file changes needed! The new structure is applied automatically based on constants in `pipeline.py`.

### Key Constants

```python
# New structure directories
DIR_INPUT = 'input'
DIR_STAGES = 'stages'
DIR_OUTPUT = 'output'
DIR_CONFIGS = 'configs'
DIR_LOGS = 'logs'

# Stage subdirectories
DIR_DESCRIPTORS_INITIAL = 'stages/01_descriptors_initial'
DIR_STRUCT_FILTERS_PRE = 'stages/02_structural_filters_pre'
DIR_STRUCT_FILTERS_POST = 'stages/03_structural_filters_post'
DIR_SYNTHESIS = 'stages/04_synthesis'
DIR_DOCKING = 'stages/05_docking'
DIR_DESCRIPTORS_FINAL = 'stages/06_descriptors_final'

# Standard file names
FILE_SAMPLED_MOLECULES = 'sampled_molecules.csv'
FILE_FINAL_MOLECULES = 'final_molecules.csv'
FILE_FILTERED_MOLECULES = 'filtered_molecules.csv'
```

## Benefits

1. **Clarity**: Anyone can understand the pipeline flow at a glance
2. **Organization**: Related files grouped logically
3. **Consistency**: Uniform naming conventions reduce confusion
4. **Maintainability**: Easier to find and debug issues
5. **Scalability**: Easy to add new stages or outputs
6. **Professional**: Follows industry best practices

## Related Files

- `src/hedge/pipeline.py`: Main constants and path handling
- `src/hedge/stages/*/main.py`: Stage entry points
- `src/hedge/stages/*/utils.py`: Stage utilities

---

**Version**: 1.0
**Last Updated**: 2025-11-05
**Status**: Implemented with backward compatibility
