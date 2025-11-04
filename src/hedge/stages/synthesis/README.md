## Synthesis Module

Evaluates synthetic accessibility of molecules using multiple metrics:
- **SA Score**: Synthetic Accessibility score (1-10, lower is better). Uses RDKit's SA Score calculator.
- **RA Score**: Retrosynthetic Accessibility score (0-1 probability, higher is better). Uses MolScore with RAScore.json config.
- **SYBA**: Synthetic Bayesian Accessibility (higher is better). Uses SYBA package.
- **Retrosynthesis**: Uses AiZynthFinder to find retrosynthetic paths and verify synthetic feasibility.

All scores are calculated for each molecule, and filters can be applied based on configurable thresholds.

### Configuration

Configure the synthesis stage in `src/hedge/configs/config_synthesis.yml`:

- `run`: Set to `True` to enable the synthesis stage
- `filter_solved_only`: If `True`, only keep molecules that were solved (retrosynthesis path found). If `False`, keep all molecules.

### Usage

**Run synthesis stage within entire pipeline:**
Set `run: True` in `src/hedge/configs/config_synthesis.yml` and run:
```bash
uv run hedge run
```

**Run synthesis stage only:**
```bash
uv run hedge run --stage synthesis
```

