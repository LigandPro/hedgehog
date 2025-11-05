import json
from pathlib import Path
from typing import Optional

import polars as pl

# Constants
MODEL_NAME_COLUMN = 'model_name'
MOL_IDX_COLUMN = 'mol_idx'
DEFAULT_MODEL_NAME = 'single'

RUN_CONFIGS_DIR = 'run_configs'
MODEL_INDEX_MAP_FILE = 'model_index_map.json'

MOL_IDX_FORMAT = 'LP-{:04d}-{:05d}'
MODEL_NUMBER_WIDTH = 4
COUNTER_WIDTH = 5

# Internal temporary columns
TEMP_MODEL_NUMBER = '_model_number'
TEMP_COUNTER = '_per_model_counter'


def assign_mol_idx(df: pl.DataFrame, run_base: Path, logger: Optional[object] = None) -> pl.DataFrame:
    """
    Assign a stable mol_idx for each row in the dataframe.

    Format: LP-<model_number>-<per_model_counter>
    - model_number: 4-digit zero-padded integer starting from 1 for each unique model
    - per_model_counter: 5-digit zero-padded integer starting from 1 for each molecule within a model

    Args:
        df: Input dataframe with molecules
        run_base: Base directory for the run (used to store model index mapping)
        logger: Optional logger instance

    Returns:
        Dataframe with mol_idx column added
    """
    df = df.clone()

    if MODEL_NAME_COLUMN not in df.columns:
        df = df.with_columns(pl.lit(DEFAULT_MODEL_NAME).alias(MODEL_NAME_COLUMN))

    present_models = df[MODEL_NAME_COLUMN].drop_nulls().cast(pl.Utf8).unique().to_list()
    model_map = {model: idx + 1 for idx, model in enumerate(present_models)}
    _save_model_index_map(run_base, model_map)

    df = df.with_columns(
        pl.col(MODEL_NAME_COLUMN)
        .replace(model_map)
        .cast(pl.Int64)
        .alias(TEMP_MODEL_NUMBER),
        (pl.int_range(pl.len()).over(MODEL_NAME_COLUMN) + 1).alias(TEMP_COUNTER),
    )

    df = df.with_columns(
        pl.format(
            "LP-{}-{}",
            pl.col(TEMP_MODEL_NUMBER).cast(pl.Utf8).str.zfill(MODEL_NUMBER_WIDTH),
            pl.col(TEMP_COUNTER).cast(pl.Utf8).str.zfill(COUNTER_WIDTH),
        ).alias(MOL_IDX_COLUMN)
    )
    df = df.drop([TEMP_MODEL_NUMBER, TEMP_COUNTER])

    return df


def _save_model_index_map(run_base: Path, model_map: dict) -> None:
    """Save model index mapping to JSON file for this run."""
    try:
        dest_dir = run_base / RUN_CONFIGS_DIR
        dest_dir.mkdir(parents=True, exist_ok=True)
        
        map_file = dest_dir / MODEL_INDEX_MAP_FILE
        with open(map_file, 'w') as f:
            json.dump(model_map, f, indent=2, sort_keys=True)
    except Exception:
        pass

