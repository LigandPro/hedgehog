import logging
from pathlib import Path
from typing import cast

import pandas as pd

# Constants
SMILES_COLUMN = "smiles"
MODEL_NAME_COLUMN = "model_name"
NAME_COLUMN = "name"
MOL_IDX_COLUMN = "mol_idx"
DEFAULT_MODEL_NAME = "single"

SUPPORTED_EXTENSIONS = ["csv", "tsv", "txt"]
MODE_SINGLE = "single_comparison"
MODE_MULTI = "multi_comparison"

RUN_CONFIGS_DIR = "run_configs"
RUN_MODELS_MAPPING_FILE = "run_models_mapping.csv"
MODEL_INDEX_MAP_FILE = "model_index_map.json"


def _find_column_case_insensitive(
    df: pd.DataFrame, column_name: str
) -> str | None:
    """Find a column by name (case-insensitive)."""
    lower_cols = {c.lower(): c for c in df.columns}
    return lower_cols.get(column_name.lower())


def _normalize_smiles_column(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize SMILES column in dataframe."""
    smiles_col = _find_column_case_insensitive(df, SMILES_COLUMN)
    if smiles_col and smiles_col != SMILES_COLUMN:
        return df.rename(columns={smiles_col: SMILES_COLUMN})
    return df


def _normalize_model_name_column(
    df: pd.DataFrame, path: str
) -> pd.DataFrame:
    """Normalize model_name column. Extract from path if missing."""
    model_col = _find_column_case_insensitive(
        df, MODEL_NAME_COLUMN
    ) or _find_column_case_insensitive(df, NAME_COLUMN)

    if model_col:
        if model_col != MODEL_NAME_COLUMN:
            df = df.rename(columns={model_col: MODEL_NAME_COLUMN})
    else:
        df[MODEL_NAME_COLUMN] = _extract_model_name_from_path(path)

    return df


def _normalize_columns(df: pd.DataFrame, path: str) -> pd.DataFrame:
    """Apply all column normalizations."""
    df = _normalize_smiles_column(df)
    df = _normalize_model_name_column(df, path)
    return df


def _apply_sampling(
    df: pd.DataFrame,
    sample_size: int | None,
    model_name: str | None = None,
) -> tuple[pd.DataFrame, dict | None]:
    """Apply sampling to dataframe if sample_size is specified."""
    if sample_size is None:
        return df, None

    if len(df) < sample_size:
        return df, {
            "model_name": model_name or "unknown",
            "requested": sample_size,
            "available": len(df),
        }

    return df.sample(sample_size, random_state=42), None


def _log_sampling_warnings(
    warnings: list[dict], logger: logging.Logger
) -> None:
    """Log sampling warnings in a formatted manner."""
    if not warnings:
        return

    if len(warnings) == 1:
        warn = warnings[0]
        logger.warning(
            "Sample size %d exceeds data size %d",
            warn["requested"],
            warn["available"],
        )
        return

    logger.warning("")
    logger.warning(
        "[yellow]Sample size exceeded for %d model(s):[/yellow]",
        len(warnings),
    )
    for warn in warnings:
        logger.warning(
            "  [dim]•[/dim] [bold]%s[/bold]: %d requested, %d available",
            warn["model_name"],
            warn["requested"],
            warn["available"],
        )


def _remove_duplicates(
    df: pd.DataFrame, logger: logging.Logger
) -> pd.DataFrame:
    """Remove duplicate molecules, logging if any were removed."""
    initial_count = len(df)

    if MODEL_NAME_COLUMN in df.columns:
        df = df.drop_duplicates(
            subset=[SMILES_COLUMN, MODEL_NAME_COLUMN]
        ).reset_index(drop=True)
        msg = "Removed %s duplicate molecules within models"
    else:
        df = df.drop_duplicates(subset=SMILES_COLUMN).reset_index(drop=True)
        msg = "Removed %s duplicate molecules"

    duplicates_removed = initial_count - len(df)
    if duplicates_removed > 0:
        logger.info(msg, duplicates_removed)

    return df


def _read_csv_with_fallback(path: str) -> pd.DataFrame:
    """Read CSV, falling back to headerless format if parsing fails."""
    try:
        return pd.read_csv(path)
    except (pd.errors.ParserError, ValueError):
        df = pd.read_csv(path, header=None)
        df.columns = [SMILES_COLUMN] + [
            f"col_{i}" for i in range(1, len(df.columns))
        ]
        return df


def _detect_mode_and_paths(
    generated_mols_path: str,
) -> tuple[str, list[str]]:
    """Detect if input is single or multi-model comparison based on file pattern."""
    path_obj = Path(generated_mols_path)

    matched = (
        [str(p) for p in path_obj.parent.glob(path_obj.name)]
        if "*" in generated_mols_path
        else []
    )

    if not matched:
        if path_obj.exists() and path_obj.is_file():
            matched = [generated_mols_path]
        else:
            msg = f"No files matched pattern: {generated_mols_path}"
            raise FileNotFoundError(msg)

    if len(matched) > 1:
        return MODE_MULTI, matched

    single_path = matched[0]
    ext = Path(single_path).suffix.lower().lstrip(".")

    if ext not in SUPPORTED_EXTENSIONS:
        return MODE_SINGLE, [single_path]

    if _file_has_multiple_models(single_path):
        return MODE_MULTI, [single_path]

    return MODE_SINGLE, [single_path]


def _file_has_multiple_models(path: str) -> bool:
    """Check if a file contains multiple distinct models."""
    try:
        df = pd.read_csv(path)
        candidate_col = _find_column_case_insensitive(
            df, MODEL_NAME_COLUMN
        ) or _find_column_case_insensitive(df, NAME_COLUMN)

        if candidate_col:
            n_distinct = df[candidate_col].nunique(dropna=True)
            return n_distinct > 1
    except (pd.errors.ParserError, ValueError, KeyError):
        pass
    return False


def _extract_model_name_from_path(path: str) -> str:
    """Extract model name from file path."""
    return path.split("/")[-1].split(".")[0]


def _load_multi_comparison_data(
    paths: list[str], sample_size: int | None, logger: logging.Logger
) -> pd.DataFrame:
    """Load and merge data from multiple model files."""
    dataframes = []
    sampling_warnings = []

    for path in paths:
        df = _read_csv_with_fallback(path)
        df = _normalize_columns(df, path)

        model_name = (
            df[MODEL_NAME_COLUMN].iloc[0]
            if MODEL_NAME_COLUMN in df.columns
            else "unknown"
        )
        df, warning_info = _apply_sampling(df, sample_size, model_name)
        if warning_info:
            sampling_warnings.append(warning_info)
        dataframes.append(df)

    _log_sampling_warnings(sampling_warnings, logger)
    return pd.concat(dataframes, axis=0, ignore_index=True)


def _load_single_comparison_data(
    single_path: str, sample_size: int | None, logger: logging.Logger
) -> pd.DataFrame:
    """Load data from a single file (single model comparison)."""
    data = _read_csv_with_fallback(single_path)
    data = _normalize_columns(data, single_path)
    data = _remove_duplicates(data, logger)

    data, warning_info = _apply_sampling(data, sample_size)
    if warning_info:
        _log_sampling_warnings([warning_info], logger)

    return data


def _load_multi_file_with_model_column(
    single_path: str, sample_size: int | None, logger: logging.Logger
) -> pd.DataFrame:
    """Load a single file containing multiple models (must have model_name column)."""
    df = pd.read_csv(single_path)
    df = _normalize_smiles_column(df)

    if not (
        _find_column_case_insensitive(df, MODEL_NAME_COLUMN)
        or _find_column_case_insensitive(df, NAME_COLUMN)
    ):
        msg = (
            f"Expected a '{MODEL_NAME_COLUMN}' or '{NAME_COLUMN}' "
            "column for multi-comparison detection."
        )
        raise ValueError(msg)

    df = _normalize_model_name_column(df, single_path)

    cols_to_keep = [SMILES_COLUMN, MODEL_NAME_COLUMN]
    if MOL_IDX_COLUMN in df.columns:
        cols_to_keep.append(MOL_IDX_COLUMN)

    df = df[cols_to_keep].dropna(subset=[SMILES_COLUMN])
    df = _remove_duplicates(df, logger)

    if sample_size is None:
        return df

    sampled = []
    sampling_warnings = []
    for model, grp in df.groupby(MODEL_NAME_COLUMN):
        sampled_grp, warning_info = _apply_sampling(grp, sample_size, model)
        sampled.append(sampled_grp)
        if warning_info:
            sampling_warnings.append(warning_info)

    _log_sampling_warnings(sampling_warnings, logger)
    return pd.concat(sampled, axis=0, ignore_index=True)


def _save_run_model_mapping(
    data: pd.DataFrame, folder_to_save: Path
) -> None:
    """Save per-run model mapping for provenance."""
    try:
        run_configs_dir = folder_to_save / RUN_CONFIGS_DIR
        run_configs_dir.mkdir(parents=True, exist_ok=True)

        if (
            MOL_IDX_COLUMN in data.columns
            and MODEL_NAME_COLUMN in data.columns
        ):
            tmp = data[[MODEL_NAME_COLUMN, MOL_IDX_COLUMN]].dropna().copy()
            tmp["model_index"] = (
                tmp[MOL_IDX_COLUMN]
                .astype(str)
                .str.split("-")
                .str[1]
                .astype(int)
            )
            run_map = tmp.groupby(MODEL_NAME_COLUMN, as_index=False)[
                "model_index"
            ].first()
            run_map = run_map[
                ["model_index", MODEL_NAME_COLUMN]
            ].sort_values("model_index")
            run_map.to_csv(
                run_configs_dir / RUN_MODELS_MAPPING_FILE, index=False
            )
    except (ValueError, KeyError, IndexError) as e:
        # Silently fail if mapping cannot be saved
        # This is not critical for the main workflow
        del e  # Avoid unused variable warning


def prepare_input_data(
    config: dict, logger: logging.Logger
) -> pd.DataFrame:
    """
    Prepare input molecular data from config-specified sources.

    Detects single vs multi-model comparison modes, loads and normalizes data,
    applies sampling if configured, and assigns molecular indices.
    """
    generated_mols_path = config["generated_mols_path"]
    folder_to_save = Path(config["folder_to_save"])
    save_sampled_mols = config.get("save_sampled_mols", False)
    sample_size: int | None = (
        cast("int | None", config.get("sample_size"))
        if save_sampled_mols
        else None
    )

    detected_mode, matched_paths = _detect_mode_and_paths(
        generated_mols_path
    )

    logger.info(
        "Loading generated molecules from %s...", generated_mols_path
    )

    if detected_mode == MODE_SINGLE:
        data = _load_single_comparison_data(
            matched_paths[0], sample_size, logger
        )
    elif detected_mode == MODE_MULTI:
        if len(matched_paths) > 1:
            logger.info(
                "Loading multi-model comparison from %s file(s)...",
                len(matched_paths),
            )
            data = _load_multi_comparison_data(
                matched_paths, sample_size, logger
            )
        else:
            data = _load_multi_file_with_model_column(
                matched_paths[0], sample_size, logger
            )
    else:
        msg = f"Invalid detected mode: {detected_mode}"
        raise ValueError(msg)

    if MODEL_NAME_COLUMN not in data.columns:
        data = _normalize_model_name_column(data, matched_paths[0])
        logger.warning(
            "model_name was missing and extracted from path: %s",
            matched_paths[0],
        )

    id_cols = [SMILES_COLUMN, MODEL_NAME_COLUMN]
    if MOL_IDX_COLUMN in data.columns:
        id_cols.append(MOL_IDX_COLUMN)
    ordered_cols = id_cols + [c for c in data.columns if c not in id_cols]
    data = data[ordered_cols]

    _save_run_model_mapping(data, folder_to_save)

    return data
