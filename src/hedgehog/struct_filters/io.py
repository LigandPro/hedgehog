"""File I/O utilities for structural filters stage."""

import os
from pathlib import Path

import pandas as pd

from hedgehog.configs.logger import logger


def process_path(folder_to_save, key_word=None):
    """Ensure path ends with '/' and create directory if needed."""
    # Accept both str and Path-like inputs.
    folder_to_save = str(folder_to_save)
    if not folder_to_save.endswith("/"):
        folder_to_save += "/"

    if key_word:
        folder_to_save += f"{key_word}/"

    Path(folder_to_save).mkdir(parents=True, exist_ok=True)
    return folder_to_save


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
