import glob
import os
from pathlib import Path
from typing import List, Tuple

import pandas as pd

from utils.mol_index import assign_mol_idx


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


def _find_smiles_column(df: pd.DataFrame) -> str:
    """Find the SMILES column in a dataframe (case-insensitive)."""
    smiles_col = next((col for col in df.columns if col.lower() == SMILES_COLUMN), None)
    return smiles_col


def _normalize_dataframe_columns(df: pd.DataFrame, has_header: bool = True) -> pd.DataFrame:
    """Normalize dataframe to have a 'smiles' column."""
    if has_header:
        smiles_col = _find_smiles_column(df)
        if smiles_col:
            return df[[smiles_col]].copy().rename(columns={smiles_col: SMILES_COLUMN})
    
    # No header or no SMILES column found
    df = df.copy()
    df.columns = [SMILES_COLUMN] + [f'col_{i}' for i in range(1, len(df.columns))]
    return df[[SMILES_COLUMN]]


def _detect_mode_and_paths(generated_mols_path: str) -> Tuple[str, List[str]]:
    """Detect if input is single or multi-model comparison based on file pattern."""
    matched = glob.glob(generated_mols_path)
    if not matched:
        # Try treating as direct file path if glob pattern didn't match
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
            df = pd.read_csv(single_path)
            lower_cols = {c.lower(): c for c in df.columns}
            
            # Check for model_name or name column
            candidate_col = lower_cols.get(MODEL_NAME_COLUMN) or lower_cols.get(NAME_COLUMN)
            
            if candidate_col:
                n_distinct = df[candidate_col].nunique(dropna=True)
                if n_distinct and n_distinct > 1:
                    return MODE_MULTI, [single_path]
        except Exception:
            pass

    return MODE_SINGLE, [single_path]


def _extract_model_name_from_path(path: str) -> str:
    """Extract model name from file path using exact format: path.split('/')[-1].split('.')[0]"""
    return path.split('/')[-1].split('.')[0]


def _load_multi_comparison_data(paths: List[str], sample_size: int, logger) -> pd.DataFrame:
    """Load and merge data from multiple model files."""
    dataframes = []
    
    for path in paths:
        try:
            df = pd.read_csv(path)
            # Preserve mol_idx and model_name if present - don't normalize away
            smiles_col = _find_smiles_column(df)
            if smiles_col:
                df = df.rename(columns={smiles_col: SMILES_COLUMN})
            else:
                df = _normalize_dataframe_columns(df, has_header=True)
        except Exception:
            df = pd.read_csv(path, header=None)
            df = _normalize_dataframe_columns(df, has_header=False)

        # Check for model_name or name column first, then extract from path if missing
        lower_cols = {c.lower(): c for c in df.columns}
        model_col = lower_cols.get(MODEL_NAME_COLUMN) or lower_cols.get(NAME_COLUMN)
        
        if model_col:
            # Rename to MODEL_NAME_COLUMN if it was NAME_COLUMN
            if model_col.lower() != MODEL_NAME_COLUMN:
                df = df.rename(columns={model_col: MODEL_NAME_COLUMN})
        else:
            # Extract model_name from path if not present in file
            model_name = _extract_model_name_from_path(path)
            df[MODEL_NAME_COLUMN] = model_name
        
        # Apply sampling if needed
        if sample_size is not None:
            if len(df) < sample_size:
                logger.warning(
                    f'Sample size {sample_size} exceeds data size {len(df)} for {df[MODEL_NAME_COLUMN].iloc[0] if MODEL_NAME_COLUMN in df.columns else "unknown"} model'
                )
            else:
                df = df.sample(sample_size, random_state=42)

        dataframes.append(df)

    return pd.concat(dataframes, axis=0, ignore_index=True)


def _load_single_comparison_data(single_path: str, sample_size: int, logger) -> pd.DataFrame:
    """Load data from a single file (single model comparison)."""
    try:
        data = pd.read_csv(single_path)
        # Preserve mol_idx and model_name if present - don't normalize away
        smiles_col = _find_smiles_column(data)
        if smiles_col:
            data = data.rename(columns={smiles_col: SMILES_COLUMN})
        else:
            data = _normalize_dataframe_columns(data, has_header=True)
    except Exception:
        data = pd.read_csv(single_path, header=None)
        data = _normalize_dataframe_columns(data, has_header=False)

    # Check for model_name or name column first, then extract from path if missing
    lower_cols = {c.lower(): c for c in data.columns}
    model_col = lower_cols.get(MODEL_NAME_COLUMN) or lower_cols.get(NAME_COLUMN)
    
    if model_col:
        # Rename to MODEL_NAME_COLUMN if it was NAME_COLUMN
        if model_col.lower() != MODEL_NAME_COLUMN:
            data = data.rename(columns={model_col: MODEL_NAME_COLUMN})
    else:
        # Extract model_name from path if not present in file
        model_name = _extract_model_name_from_path(single_path)
        data[MODEL_NAME_COLUMN] = model_name

    data = data.drop_duplicates(subset=SMILES_COLUMN).reset_index(drop=True)
    
    if sample_size is not None:
        if len(data) < sample_size:
            logger.warning(f'Sample size {sample_size} exceeds data size {len(data)}')
        else:
            data = data.sample(sample_size, random_state=42)
    
    return data


def _load_multi_file_with_model_column(single_path: str, sample_size: int, logger) -> pd.DataFrame:
    """Load a single file containing multiple models (with model_name column)."""
    df = pd.read_csv(single_path)
    
    # Normalize SMILES column
    smiles_col = _find_smiles_column(df)
    if smiles_col:
        other_cols = [c for c in df.columns if c != smiles_col]
        df = df[[smiles_col] + other_cols].copy()
        df = df.rename(columns={smiles_col: SMILES_COLUMN})
    else:
        df = pd.read_csv(single_path, header=None)
        df.columns = [SMILES_COLUMN] + [f'col_{i}' for i in range(1, len(df.columns))]

    # Find model column
    lower_cols = {c.lower(): c for c in df.columns}
    model_col = lower_cols.get(MODEL_NAME_COLUMN) or lower_cols.get(NAME_COLUMN)
    
    if not model_col:
        raise ValueError(
            f"Expected a '{MODEL_NAME_COLUMN}' or '{NAME_COLUMN}' column for multi-comparison detection."
        )
    
    df = df.rename(columns={model_col: MODEL_NAME_COLUMN})
    
    # Preserve mol_idx if present, keep all identity columns
    cols_to_keep = [SMILES_COLUMN, MODEL_NAME_COLUMN]
    if MOL_IDX_COLUMN in df.columns:
        cols_to_keep.append(MOL_IDX_COLUMN)
    
    df = df[cols_to_keep].dropna(subset=[SMILES_COLUMN])
    df = df.drop_duplicates(subset=[SMILES_COLUMN, MODEL_NAME_COLUMN]).reset_index(drop=True)

    # Sample per model if needed
    if sample_size is not None:
        sampled = []
        for model, grp in df.groupby(MODEL_NAME_COLUMN):
            if len(grp) < sample_size:
                logger.warning(
                    f'Sample size {sample_size} exceeds data size {len(grp)} for {model} model'
                )
                sampled.append(grp)
            else:
                sampled.append(grp.sample(sample_size, random_state=42))
        df = pd.concat(sampled, axis=0, ignore_index=True)
    
    return df


def _save_run_model_mapping(data: pd.DataFrame, folder_to_save: Path):
    """Save per-run model mapping for provenance."""
    try:
        run_configs_dir = folder_to_save / RUN_CONFIGS_DIR
        run_configs_dir.mkdir(parents=True, exist_ok=True)
        
        if MOL_IDX_COLUMN in data.columns and MODEL_NAME_COLUMN in data.columns:
            tmp = data[[MODEL_NAME_COLUMN, MOL_IDX_COLUMN]].dropna().copy()
            tmp['model_index'] = tmp[MOL_IDX_COLUMN].astype(str).str.split('-').str[1].astype(int)
            run_map = tmp.groupby(MODEL_NAME_COLUMN, as_index=False)['model_index'].first()
            run_map = run_map[['model_index', MODEL_NAME_COLUMN]].sort_values('model_index')
            run_map.to_csv(run_configs_dir / RUN_MODELS_MAPPING_FILE, index=False)
    except Exception:
        pass


def prepare_input_data(config: dict, logger) -> pd.DataFrame:
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

    logger.info(f"Detected mode: {detected_mode}")
    logger.info(f'Loading generated molecules from {generated_mols_path}...')

    # Load data based on detected mode
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

    # Extract model_name from path if not present in input file
    # This should already be handled in the load functions, but double-check as safety
    if MODEL_NAME_COLUMN not in data.columns:
        # Check for name column as fallback
        lower_cols = {c.lower(): c for c in data.columns}
        if NAME_COLUMN.lower() in lower_cols:
            data = data.rename(columns={lower_cols[NAME_COLUMN.lower()]: MODEL_NAME_COLUMN})
        else:
            # Extract from path (should only happen if somehow missed in load functions)
            model_name = _extract_model_name_from_path(matched_paths[0])
            data[MODEL_NAME_COLUMN] = model_name
            logger.warning(f"model_name extracted from path in prepare_input_data: {model_name}")

    # Preserve mol_idx from input file - only assign if missing
    if MOL_IDX_COLUMN not in data.columns or data[MOL_IDX_COLUMN].isna().all():
        try:
            data = assign_mol_idx(data, run_base=folder_to_save, logger=logger)
        except Exception as e:
            logger.warning(f'Could not assign mol_idx: {e}')
    else:
        logger.debug('mol_idx already present in input data, preserving original values')

    # Ensure consistent column order
    id_cols = [SMILES_COLUMN, MODEL_NAME_COLUMN, MOL_IDX_COLUMN]
    ordered_cols = id_cols + [c for c in data.columns if c not in id_cols]
    data = data[ordered_cols]

    # Save model mapping for this run
    _save_run_model_mapping(data, folder_to_save)

    return data


