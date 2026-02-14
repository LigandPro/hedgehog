from __future__ import annotations

from pathlib import Path

from hedgehog.configs.logger import load_config, logger
from hedgehog.stages.molPrep.utils import run_mol_prep
from hedgehog.stages.structFilters.utils import process_path


def main(data, config: dict, subfolder: str | None = None):
    """Run Datamol-based molecule preparation (MolPrep).

    This stage standardizes molecules before descriptor calculation:
    - salts/solvents removal + largest fragment selection
    - metal disconnection
    - uncharging + normalization/reionization
    - tautomer canonicalization (standardize_smiles)
    - stereochemistry removal
    - strict filtering (allowed atoms, radicals, isotopes, single-fragment, neutrality)

    Output:
      - filtered_molecules.csv
      - failed_molecules.csv
      - metrics.csv
      - duplicates_removed.csv (optional)
    """
    if data is None or len(data) == 0:
        logger.warning("No molecules provided for MolPrep. Skipping.")
        return None

    folder_to_save = Path(process_path(config["folder_to_save"]))
    subfolder = subfolder or str(Path("stages") / "00_mol_prep")
    out_dir = folder_to_save / subfolder
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg = load_config(config["config_mol_prep"])
    return run_mol_prep(data, cfg, out_dir)

