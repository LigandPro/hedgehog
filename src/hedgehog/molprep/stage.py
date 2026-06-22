from __future__ import annotations

from pathlib import Path

from hedgehog._constants import KEY_FOLDER_TO_SAVE
from hedgehog.configs.logger import load_config, logger
from hedgehog.large_dataset import is_large_dataset_mode
from hedgehog.molprep.large import run_large
from hedgehog.molprep.orchestrator import run_mol_prep
from hedgehog.utils.paths import process_path


def run(data, config: dict, subfolder: str | None = None, reporter=None):
    """Run Datamol-based molecule preparation (MolPrep).

    This stage standardizes molecules before descriptor calculation:
    - salts/solvents removal + largest fragment selection
    - metal disconnection
    - uncharging + normalization/reionization
    - tautomer canonicalization (standardize_smiles)
    - stereochemistry removal
    - strict filtering (allowed atoms, radicals, isotopes, single-fragment)

    Output:
      - filtered_molecules.csv
      - failed_molecules.csv
      - metrics.csv
      - molprep_detail.csv (per-molecule step/filter pass flags)
      - duplicates_removed.csv (optional)
    """
    folder_to_save = Path(process_path(config[KEY_FOLDER_TO_SAVE]))
    subfolder = subfolder or str(Path("stages") / "01_mol_prep")
    out_dir = folder_to_save / subfolder
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg = load_config(config["config_mol_prep"])
    if is_large_dataset_mode(config):
        return run_large(config, cfg, out_dir, reporter=reporter)

    if data is None or len(data) == 0:
        logger.warning("No molecules provided for MolPrep. Skipping.")
        return None
    return run_mol_prep(data, cfg, out_dir, reporter=reporter)
