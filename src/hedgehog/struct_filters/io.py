"""File I/O utilities for structural filters stage."""

from pathlib import Path

import pandas as pd

from hedgehog.utils.paths import process_path as _shared_process_path


def process_path(folder_to_save, key_word=None):
    """Backward-compatible wrapper around shared process_path helper."""
    return _shared_process_path(folder_to_save, key_word)


def inject_identity_columns_to_all_csvs(config, stage_dir):
    """Ensure identity columns are ordered consistently in all CSVs."""
    from hedgehog._constants import KEY_FOLDER_TO_SAVE

    base_folder = Path(process_path(config[KEY_FOLDER_TO_SAVE]))
    target_folder = base_folder / stage_dir

    csv_paths = [str(p) for p in target_folder.glob("*.csv")]
    for path in csv_paths:
        try:
            df = pd.read_csv(path)
            if "smiles" not in df.columns:
                continue

            identity_order = ["smiles", "model_name", "mol_idx"]
            ordered = [c for c in identity_order if c in df.columns] + [
                c for c in df.columns if c not in identity_order
            ]
            df = df[ordered]
            df.to_csv(path, index=False)
        except Exception:
            continue
