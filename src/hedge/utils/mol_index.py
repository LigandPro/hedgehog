import json
from pathlib import Path

import pandas as pd

# Constants
MODEL_NAME_COLUMN = "model_name"
MOL_IDX_COLUMN = "mol_idx"
DEFAULT_MODEL_NAME = "single"

RUN_CONFIGS_DIR = "run_configs"
MODEL_INDEX_MAP_FILE = "model_index_map.json"

MOL_IDX_FORMAT = "LP-{:04d}-{:05d}"
MODEL_NUMBER_WIDTH = 4
COUNTER_WIDTH = 5

# Internal temporary columns
TEMP_MODEL_NUMBER = "_model_number"
TEMP_COUNTER = "_per_model_counter"


def assign_mol_idx(df: pd.DataFrame, run_base: Path, logger: object | None = None) -> pd.DataFrame:
    """
    Assign a stable mol_idx for each row in the dataframe.
    
    Format: LP-<model_number>-<per_model_counter>
    - model_number: 4-digit zero-padded integer starting from 1 for each unique model
    - per_model_counter: 5-digit zero-padded integer starting from 1 for each molecule within a model
    
    Args:
        df: Input dataframe with molecules
        run_base: Base directory for the run (used to store model index mapping)
        logger: Optional logger instance
    
    Returns
    -------
        Dataframe with mol_idx column added
    """
    df = df.copy()

    if MODEL_NAME_COLUMN not in df.columns:
        df[MODEL_NAME_COLUMN] = DEFAULT_MODEL_NAME

    present_models = [str(m) for m in pd.unique(df[MODEL_NAME_COLUMN].dropna().astype(str))]
    model_map = {model: idx + 1 for idx, model in enumerate(present_models)}
    _save_model_index_map(run_base, model_map)

    df[TEMP_MODEL_NUMBER] = df[MODEL_NAME_COLUMN].map(model_map)
    df[TEMP_COUNTER] = df.groupby(MODEL_NAME_COLUMN).cumcount() + 1

    df[MOL_IDX_COLUMN] = df.apply(lambda row: MOL_IDX_FORMAT.format(int(row[TEMP_MODEL_NUMBER]), int(row[TEMP_COUNTER])), axis=1)
    df = df.drop(columns=[TEMP_MODEL_NUMBER, TEMP_COUNTER])

    return df


def _save_model_index_map(run_base: Path, model_map: dict) -> None:
    """Save model index mapping to JSON file for this run."""
    try:
        dest_dir = run_base / RUN_CONFIGS_DIR
        dest_dir.mkdir(parents=True, exist_ok=True)

        map_file = dest_dir / MODEL_INDEX_MAP_FILE
        with open(map_file, "w") as f:
            json.dump(model_map, f, indent=2, sort_keys=True)
    except Exception:
        pass


