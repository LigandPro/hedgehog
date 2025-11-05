import glob
import os
from pathlib import Path
from typing import List, Tuple, Optional

import polars as pl

from hedge.utils.mol_index import assign_mol_idx


# Constants
SMILES_COLUMN = 'smiles'
MODEL_NAME_COLUMN = 'model_name'
NAME_COLUMN = 'name'
MOL_IDX_COLUMN = 'mol_idx'
DEFAULT_MODEL_NAME = 'single'

SUPPORTED_EXTENSIONS = ['csv', 'tsv', 'txt']
MODE_SINGLE = 'single_comparison'
MODE_MULTI = 'multi_comparison'

RUN_CONFIGS_DIR = 'run_configs'
RUN_MODELS_MAPPING_FILE = 'run_models_mapping.csv'
MODEL_INDEX_MAP_FILE = 'model_index_map.json'


def _find_smiles_column(df: pl.DataFrame) -> str:
    """Find the SMILES column in a dataframe (case-insensitive)."""
    return next((col for col in df.columns if col.lower() == SMILES_COLUMN), None)


def _normalize_smiles_column(df: pl.DataFrame, path: str) -> pl.DataFrame:
    """Normalize SMILES column in dataframe."""
    smiles_col = _find_smiles_column(df)
    if smiles_col:
        return df.rename({smiles_col: SMILES_COLUMN})
    return df


def _normalize_model_name_column(df: pl.DataFrame, path: str) -> pl.DataFrame:
    """Normalize model_name column: check for model_name/name, extract from path if missing."""
    lower_cols = {c.lower(): c for c in df.columns}
    model_col = lower_cols.get(MODEL_NAME_COLUMN) or lower_cols.get(NAME_COLUMN)

    if model_col:
        if model_col.lower() != MODEL_NAME_COLUMN:
            df = df.rename({model_col: MODEL_NAME_COLUMN})
    else:
        df = df.with_columns(pl.lit(_extract_model_name_from_path(path)).alias(MODEL_NAME_COLUMN))

    return df


def _apply_sampling(df: pl.DataFrame, sample_size: Optional[int], logger, model_name: Optional[str] = None) -> pl.DataFrame:
    """Apply sampling to dataframe if sample_size is specified."""
    if sample_size is None:
        return df

    if len(df) < sample_size:
        model_info = f" for {model_name} model" if model_name else ""
        logger.warning(f'Sample size {sample_size} exceeds data size {len(df)}{model_info}')
        return df

    return df.sample(n=sample_size, seed=42)


def _detect_mode_and_paths(generated_mols_path: str) -> Tuple[str, List[str]]:
    """Detect if input is single or multi-model comparison based on file pattern."""
    matched = glob.glob(generated_mols_path)
    if not matched:
        if os.path.exists(generated_mols_path) and os.path.isfile(generated_mols_path):
            matched = [generated_mols_path]
        else:
            raise FileNotFoundError(f"No files matched pattern: {generated_mols_path}")

    if len(matched) > 1:
        return MODE_MULTI, matched

    single_path = matched[0]
    ext = Path(single_path).suffix.lower().lstrip('.')

    if ext in SUPPORTED_EXTENSIONS:
        try:
            df = pl.read_csv(single_path)
            lower_cols = {c.lower(): c for c in df.columns}

            candidate_col = lower_cols.get(MODEL_NAME_COLUMN) or lower_cols.get(NAME_COLUMN)

            if candidate_col:
                n_distinct = df[candidate_col].n_unique()
                if n_distinct and n_distinct > 1:
                    return MODE_MULTI, [single_path]
        except Exception:
            pass

    return MODE_SINGLE, [single_path]


def _extract_model_name_from_path(path: str) -> str:
    """Extract model name from file path using exact format: path.split('/')[-1].split('.')[0]"""
    return path.split('/')[-1].split('.')[0]


def _load_multi_comparison_data(paths: List[str], sample_size: int, logger) -> pl.DataFrame:
    """Load and merge data from multiple model files."""
    dataframes = []

    for path in paths:
        try:
            df = pl.read_csv(path)
        except Exception:
            df = pl.read_csv(path, has_header=False)
            new_cols = [SMILES_COLUMN] + [f'col_{i}' for i in range(1, len(df.columns))]
            df = df.rename(dict(zip(df.columns, new_cols)))

        df = _normalize_smiles_column(df, path)
        df = _normalize_model_name_column(df, path)
        model_name = df[MODEL_NAME_COLUMN][0] if MODEL_NAME_COLUMN in df.columns else "unknown"
        df = _apply_sampling(df, sample_size, logger, model_name)
        dataframes.append(df)

    return pl.concat(dataframes, how="vertical")


def _load_single_comparison_data(single_path: str, sample_size: int, logger) -> pl.DataFrame:
    """Load data from a single file (single model comparison)."""
    try:
        data = pl.read_csv(single_path)
    except Exception:
        data = pl.read_csv(single_path, has_header=False)
        new_cols = [SMILES_COLUMN] + [f'col_{i}' for i in range(1, len(data.columns))]
        data = data.rename(dict(zip(data.columns, new_cols)))

    data = _normalize_smiles_column(data, single_path)
    data = _normalize_model_name_column(data, single_path)

    if MODEL_NAME_COLUMN in data.columns:
        initial_count = len(data)
        data = data.unique(subset=[SMILES_COLUMN, MODEL_NAME_COLUMN])
        duplicates_removed = initial_count - len(data)
        if duplicates_removed > 0:
            logger.info(f"Removed {duplicates_removed} duplicate molecules within models")
    else:
        initial_count = len(data)
        data = data.unique(subset=[SMILES_COLUMN])
        duplicates_removed = initial_count - len(data)
        if duplicates_removed > 0:
            logger.info(f"Removed {duplicates_removed} duplicate molecules")

    data = _apply_sampling(data, sample_size, logger)
    return data


def _load_multi_file_with_model_column(single_path: str, sample_size: int, logger) -> pl.DataFrame:
    """Load a single file containing multiple models (with model_name column)."""
    df = pl.read_csv(single_path)
    df = _normalize_smiles_column(df, single_path)

    lower_cols = {c.lower(): c for c in df.columns}
    if MODEL_NAME_COLUMN.lower() not in lower_cols and NAME_COLUMN.lower() not in lower_cols:
        raise ValueError(f"Expected a '{MODEL_NAME_COLUMN}' or '{NAME_COLUMN}' column for multi-comparison detection.")

    df = _normalize_model_name_column(df, single_path)
    cols_to_keep = [SMILES_COLUMN, MODEL_NAME_COLUMN]
    if MOL_IDX_COLUMN in df.columns:
        cols_to_keep.append(MOL_IDX_COLUMN)

    df = df.select(cols_to_keep).drop_nulls(subset=[SMILES_COLUMN])
    initial_count = len(df)
    df = df.unique(subset=[SMILES_COLUMN, MODEL_NAME_COLUMN])
    duplicates_removed = initial_count - len(df)
    if duplicates_removed > 0:
        logger.info(f"Removed {duplicates_removed} duplicate molecules within models")

    if sample_size is not None:
        sampled = []
        for model, grp in df.group_by(MODEL_NAME_COLUMN):
            sampled.append(_apply_sampling(grp, sample_size, logger, model))
        df = pl.concat(sampled, how="vertical")

    return df


def _save_run_model_mapping(data: pl.DataFrame, folder_to_save: Path):
    """Save per-run model mapping for provenance."""
    try:
        run_configs_dir = folder_to_save / RUN_CONFIGS_DIR
        run_configs_dir.mkdir(parents=True, exist_ok=True)

        if MOL_IDX_COLUMN in data.columns and MODEL_NAME_COLUMN in data.columns:
            tmp = data.select([MODEL_NAME_COLUMN, MOL_IDX_COLUMN]).drop_nulls()
            tmp = tmp.with_columns(
                pl.col(MOL_IDX_COLUMN).cast(pl.Utf8).str.split('-').list.get(1).cast(pl.Int64).alias('model_index')
            )
            run_map = tmp.group_by(MODEL_NAME_COLUMN).agg(pl.col('model_index').first())
            run_map = run_map.select(['model_index', MODEL_NAME_COLUMN]).sort('model_index')
            run_map.write_csv(run_configs_dir / RUN_MODELS_MAPPING_FILE)
    except Exception:
        pass


def prepare_input_data(config: dict, logger) -> pl.DataFrame:
    """
    Prepare input molecular data from config-specified sources.

    Detects single vs multi-model comparison modes, loads and normalizes data,
    applies sampling if configured, and assigns molecular indices.
    """
    generated_mols_path = config['generated_mols_path']
    folder_to_save = Path(config['folder_to_save'])
    save_sampled_mols = config.get('save_sampled_mols', False)
    sample_size = config.get('sample_size') if save_sampled_mols else None

    detected_mode, matched_paths = _detect_mode_and_paths(generated_mols_path)

    logger.info(f'Loading generated molecules from {generated_mols_path}...')

    if detected_mode == MODE_SINGLE:
        data = _load_single_comparison_data(matched_paths[0], sample_size, logger)
    elif detected_mode == MODE_MULTI:
        if len(matched_paths) > 1:
            logger.info(f'Loading multi-model comparison from {len(matched_paths)} file(s)...')
            data = _load_multi_comparison_data(matched_paths, sample_size, logger)
        else:
            data = _load_multi_file_with_model_column(matched_paths[0], sample_size, logger)
    else:
        raise ValueError(f'Invalid detected mode: {detected_mode}')

    if MODEL_NAME_COLUMN not in data.columns:
        data = _normalize_model_name_column(data, matched_paths[0])
        logger.warning(f"model_name was missing and extracted from path: {matched_paths[0]}")

    id_cols = [SMILES_COLUMN, MODEL_NAME_COLUMN]
    if MOL_IDX_COLUMN in data.columns:
        id_cols.append(MOL_IDX_COLUMN)
    ordered_cols = id_cols + [c for c in data.columns if c not in id_cols]
    data = data.select(ordered_cols)

    _save_run_model_mapping(data, folder_to_save)

    return data


