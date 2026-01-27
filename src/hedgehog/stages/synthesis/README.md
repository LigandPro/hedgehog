## Synthesis Module

Evaluates synthetic accessibility of molecules using multiple metrics:
- **SA Score**: Synthetic Accessibility score (1-10, lower is better). Uses RDKit's SA Score calculator.
- **RA Score**: Retrosynthetic Accessibility score (0-1 probability, higher is better). Uses MolScore with RAScore.json config.
- **SYBA**: Synthetic Bayesian Accessibility (higher is better). Uses SYBA package.
- **Retrosynthesis**: Uses AiZynthFinder to find retrosynthetic paths and verify synthetic feasibility.

All scores are calculated for each molecule, and filters can be applied based on configurable thresholds.

### Configuration

Configure the synthesis stage in `src/hedgehog/configs/config_synthesis.yml`:

- `run`: Set to `true` to enable the synthesis stage
- `run_retrosynthesis`: Set to `true` to run AiZynthFinder retrosynthesis analysis. If `false`, only SA/SYBA/RA scores will be calculated (faster, no route search).
- `filter_solved_only`: If `true`, only keep molecules that were solved (retrosynthesis path found). If `false`, keep all molecules. Only applies when `run_retrosynthesis` is `true`.

**Score thresholds:**
- `sa_score_min`/`sa_score_max`: SA Score range (default: 0-4.5)
- `syba_score_min`/`syba_score_max`: SYBA Score range (default: 0-inf)
- `ra_score_min`/`ra_score_max`: RA Score range (default: 0.5-1)

### Usage

**Run synthesis stage within entire pipeline:**

Set `run: true` in `src/hedgehog/configs/config_synthesis.yml` and run:
```bash
uv run hedgehog run
```

**Run synthesis stage only:**
```bash
uv run hedgehog run --stage synthesis

# Alternatively, using the short-name alias:
uv run hedge run --stage synthesis
```

**Run without retrosynthesis (scores only):**

Set `run_retrosynthesis: false` in config to skip AiZynthFinder analysis and only calculate SA/SYBA/RA scores.

### Output Files

- `synthesis_scores.csv` - SA Score, RA Score, SYBA scores for all molecules
- `synthesis_extended.csv` - Scores with retrosynthesis results (when retrosynthesis enabled)
- `filtered_molecules.csv` - Molecules passing all filters
- `input_smiles.smi` - Input file for AiZynthFinder
- `retrosynthesis_results.json` - Raw AiZynthFinder output
