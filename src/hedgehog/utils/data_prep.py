import logging
from collections import Counter
from pathlib import Path
from typing import cast

import pandas as pd

from hedgehog._constants import KEY_FOLDER_TO_SAVE

# Constants
SMILES_COLUMN = "smiles"
MODEL_NAME_COLUMN = "model_name"
NAME_COLUMN = "name"
MOL_IDX_COLUMN = "mol_idx"
DEFAULT_MODEL_NAME = "single"

SUPPORTED_EXTENSIONS = ["csv", "tsv", "txt", "sdf", "smi", "smiles"]
MODE_SINGLE = "single_comparison"
MODE_MULTI = "multi_comparison"

RUN_CONFIGS_DIR = "configs"
RUN_MODELS_MAPPING_FILE = "run_models_mapping.csv"
MODEL_INDEX_MAP_FILE = "model_index_map.json"


def _find_column_case_insensitive(df: pd.DataFrame, column_name: str) -> str | None:
    """Find a column by name (case-insensitive)."""
    lower_cols = {c.lower(): c for c in df.columns}
    return lower_cols.get(column_name.lower())


def _normalize_smiles_column(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize SMILES column in dataframe."""
    smiles_col = _find_column_case_insensitive(df, SMILES_COLUMN)
    if smiles_col and smiles_col != SMILES_COLUMN:
        return df.rename(columns={smiles_col: SMILES_COLUMN})
    return df


def _normalize_model_name_column(df: pd.DataFrame, path: str) -> pd.DataFrame:
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


def _log_sampling_warnings(warnings: list[dict], logger: logging.Logger) -> None:
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


def _remove_duplicates(df: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
    """Remove duplicate molecules within each model, logging if any were removed.

    Important: Deduplication is done per-model (smiles + model_name) to preserve
    model-specific statistics. If model_name is missing, it's added as a fallback.
    """
    initial_count = len(df)

    # Ensure model_name column exists for proper per-model deduplication
    if MODEL_NAME_COLUMN not in df.columns:
        logger.warning(
            "model_name column missing during deduplication. Adding default model_name='single'. "
            "This should not happen in normal pipeline execution."
        )
        df[MODEL_NAME_COLUMN] = DEFAULT_MODEL_NAME

    # Always deduplicate within models (smiles + model_name)
    df = df.drop_duplicates(subset=[SMILES_COLUMN, MODEL_NAME_COLUMN]).reset_index(
        drop=True
    )
    msg = "Removed %s duplicate molecules within models"

    duplicates_removed = initial_count - len(df)
    if duplicates_removed > 0:
        logger.info(msg, duplicates_removed)

    return df


def _finalize_identity_columns(df: pd.DataFrame, path: str) -> pd.DataFrame:
    """Validate and trim input data to the identity columns used by the pipeline."""
    if SMILES_COLUMN not in df.columns:
        msg = f"Input file {path} is missing required '{SMILES_COLUMN}' column."
        raise ValueError(msg)

    normalized = df.copy()
    for col in (SMILES_COLUMN, MODEL_NAME_COLUMN, MOL_IDX_COLUMN):
        if col in normalized.columns:
            series = normalized[col].astype("string").str.strip()
            normalized[col] = series.mask(series.eq(""), pd.NA)

    missing_smiles = normalized[SMILES_COLUMN].isna()
    if missing_smiles.any():
        if MODEL_NAME_COLUMN in normalized.columns:
            model_counts = Counter(
                normalized.loc[missing_smiles, MODEL_NAME_COLUMN]
                .fillna("unknown")
                .astype(str)
            )
            preview = ", ".join(
                f"{model}={count}" for model, count in model_counts.most_common(5)
            )
        else:
            preview = "model_name unavailable"

        msg = (
            f"Input file {path} contains {int(missing_smiles.sum())} row(s) with "
            f"empty '{SMILES_COLUMN}' values ({preview})."
        )
        raise ValueError(msg)

    cols_to_keep = [SMILES_COLUMN, MODEL_NAME_COLUMN]
    if MOL_IDX_COLUMN in normalized.columns:
        cols_to_keep.append(MOL_IDX_COLUMN)
    return normalized[cols_to_keep].copy()


def _read_csv_with_fallback(path: str) -> pd.DataFrame:
    """Read CSV, recognizing a one-column headerless SMILES file."""
    try:
        df = pd.read_csv(path)
    except (pd.errors.ParserError, ValueError):
        df = pd.read_csv(path, header=None)
        df.columns = [SMILES_COLUMN] + [f"col_{i}" for i in range(1, len(df.columns))]
        return df

    if (
        _find_column_case_insensitive(df, SMILES_COLUMN) is None
        and len(df.columns) == 1
    ):
        headerless = pd.read_csv(path, header=None)
        headerless.columns = [SMILES_COLUMN]
        return headerless
    return df


def _read_sdf(path: str) -> pd.DataFrame:
    """Read SDF into a dataframe with identity columns and SDF properties."""
    try:
        from rdkit import Chem
    except ImportError as err:
        raise RuntimeError("RDKit is required to read SDF inputs") from err

    model_name_default = _extract_model_name_from_path(path)
    rows: list[dict] = []
    supplier = Chem.SDMolSupplier(path, removeHs=False)
    for mol in supplier:
        if mol is None:
            continue

        try:
            bare = Chem.RemoveHs(mol)
            smiles = Chem.MolToSmiles(bare)
        except Exception:
            continue
        if not smiles:
            continue

        row: dict[str, object] = {SMILES_COLUMN: smiles}

        props = mol.GetPropsAsDict(includePrivate=False, includeComputed=False)
        for key, value in props.items():
            row[str(key)] = value

        model_name_val = row.get(MODEL_NAME_COLUMN)
        if model_name_val is None or str(model_name_val).strip() == "":
            row[MODEL_NAME_COLUMN] = model_name_default

        mol_idx_val = row.get(MOL_IDX_COLUMN)
        if mol_idx_val is None or str(mol_idx_val).strip() == "":
            title = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
            if title:
                row[MOL_IDX_COLUMN] = title

        rows.append(row)

    if not rows:
        raise ValueError(f"No readable molecules found in SDF: {path}")
    return pd.DataFrame(rows)


def _read_smi(path: str) -> pd.DataFrame:
    """Read SMI/SMILES text file into a dataframe with identity columns."""
    model_name_default = _extract_model_name_from_path(path)
    rows: list[dict] = []
    with open(path, encoding="utf-8", errors="ignore") as handle:
        for idx, line in enumerate(handle):
            text = line.strip()
            if not text:
                continue
            parts = text.split()
            smiles = parts[0].strip()
            if not smiles:
                continue
            mol_idx = parts[1].strip() if len(parts) > 1 else f"{idx + 1}"
            rows.append(
                {
                    SMILES_COLUMN: smiles,
                    MODEL_NAME_COLUMN: model_name_default,
                    MOL_IDX_COLUMN: mol_idx,
                }
            )
    if not rows:
        raise ValueError(f"No readable SMILES found in file: {path}")
    return pd.DataFrame(rows)


def _read_input_file(path: str) -> pd.DataFrame:
    """Read supported input file types into a dataframe."""
    ext = Path(path).suffix.lower().lstrip(".")
    if ext == "sdf":
        return _read_sdf(path)
    if ext in {"smi", "smiles"}:
        return _read_smi(path)
    return _read_csv_with_fallback(path)


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
        elif path_obj.exists() and path_obj.is_dir():
            # Handle directory: find all supported files in the directory
            all_extensions = SUPPORTED_EXTENSIONS

            matched = [
                str(p)
                for p in path_obj.iterdir()
                if p.is_file() and p.suffix.lower().lstrip(".") in all_extensions
            ]
            if not matched:
                msg = f"No supported files found in directory: {generated_mols_path}"
                raise FileNotFoundError(msg)
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


def resolve_molecule_input_paths(
    generated_mols_path: str | None,
    generated_mols_paths: list[str] | None = None,
) -> tuple[str, list[str]]:
    """Resolve molecule inputs from an explicit path list or a path/glob/dir."""
    explicit_paths = [
        str(Path(path).expanduser())
        for path in (generated_mols_paths or [])
        if str(path).strip()
    ]
    if explicit_paths:
        missing = [path for path in explicit_paths if not Path(path).exists()]
        if missing:
            raise FileNotFoundError(
                "No such molecule input file(s): " + ", ".join(missing)
            )
        if len(explicit_paths) == 1:
            return _detect_mode_and_paths(explicit_paths[0])
        return MODE_MULTI, explicit_paths

    if not generated_mols_path:
        raise FileNotFoundError("No molecule input path provided")
    return _detect_mode_and_paths(generated_mols_path)


def _file_has_multiple_models(path: str) -> bool:
    """Check if a file contains multiple distinct models."""
    try:
        df = _read_input_file(path)
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
    return Path(path).stem


def _load_multi_comparison_data(
    paths: list[str], sample_size: int | None, logger: logging.Logger
) -> pd.DataFrame:
    """Load and merge data from multiple model files."""
    dataframes = []
    sampling_warnings = []

    for path in paths:
        df = _read_input_file(path)
        df = _normalize_columns(df, path)
        df = _finalize_identity_columns(df, path)

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
    data = _read_input_file(single_path)
    data = _normalize_columns(data, single_path)
    data = _finalize_identity_columns(data, single_path)
    data = _remove_duplicates(data, logger)

    data, warning_info = _apply_sampling(data, sample_size)
    if warning_info:
        _log_sampling_warnings([warning_info], logger)

    return data


def _load_multi_file_with_model_column(
    single_path: str, sample_size: int | None, logger: logging.Logger
) -> pd.DataFrame:
    """Load a single file containing multiple models (must have model_name column)."""
    df = _read_input_file(single_path)
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
    df = _finalize_identity_columns(df, single_path)
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


def _save_run_model_mapping(data: pd.DataFrame, folder_to_save: Path) -> None:
    """Save per-run model mapping for provenance."""
    try:
        run_configs_dir = folder_to_save / RUN_CONFIGS_DIR
        run_configs_dir.mkdir(parents=True, exist_ok=True)

        if MOL_IDX_COLUMN in data.columns and MODEL_NAME_COLUMN in data.columns:
            tmp = data[[MODEL_NAME_COLUMN, MOL_IDX_COLUMN]].dropna().copy()
            tmp["model_index"] = (
                tmp[MOL_IDX_COLUMN].astype(str).str.split("-").str[1].astype(int)
            )
            run_map = tmp.groupby(MODEL_NAME_COLUMN, as_index=False)[
                "model_index"
            ].first()
            run_map = run_map[["model_index", MODEL_NAME_COLUMN]].sort_values(
                "model_index"
            )
            run_map.to_csv(run_configs_dir / RUN_MODELS_MAPPING_FILE, index=False)
    except (ValueError, KeyError, IndexError) as e:
        # Silently fail if mapping cannot be saved
        # This is not critical for the main workflow
        del e  # Avoid unused variable warning


def prepare_input_data(config: dict, logger: logging.Logger) -> pd.DataFrame:
    """
    Prepare input molecular data from config-specified sources.

    Detects single vs multi-model comparison modes, loads and normalizes data,
    applies sampling if configured, and assigns molecular indices.
    """
    generated_mols_path = config.get("generated_mols_path")
    folder_to_save = Path(config[KEY_FOLDER_TO_SAVE])
    sample_size = cast(int | None, config.get("sample_size"))

    detected_mode, matched_paths = resolve_molecule_input_paths(
        generated_mols_path,
        config.get("generated_mols_paths"),
    )
    config["generated_mols_paths"] = matched_paths
    if not generated_mols_path:
        config["generated_mols_path"] = (
            matched_paths[0] if len(matched_paths) == 1 else str(Path(matched_paths[0]).parent)
        )
        generated_mols_path = config["generated_mols_path"]

    logger.info("Loading generated molecules from %s...", generated_mols_path)

    if detected_mode == MODE_SINGLE:
        data = _load_single_comparison_data(matched_paths[0], sample_size, logger)
    elif detected_mode == MODE_MULTI:
        if len(matched_paths) > 1:
            logger.info(
                "Loading multi-model comparison from %s file(s)...",
                len(matched_paths),
            )
            data = _load_multi_comparison_data(matched_paths, sample_size, logger)
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


def materialize_named_sdf_from_sources(
    source_paths: list[str],
    identity_df: pd.DataFrame,
    output_sdf: Path,
) -> Path:
    """Write a combined SDF preserving source coordinates and LP molecule names."""
    try:
        from rdkit import Chem
    except ImportError as err:
        raise RuntimeError(
            "RDKit is required to materialize docking SDF inputs"
        ) from err

    if identity_df.empty:
        raise ValueError("Cannot materialize SDF from an empty identity table")

    required = {SMILES_COLUMN, MODEL_NAME_COLUMN, MOL_IDX_COLUMN}
    missing = required.difference(identity_df.columns)
    if missing:
        raise ValueError(
            "Identity table is missing required columns: "
            + ", ".join(sorted(missing))
        )

    from collections import defaultdict, deque

    pools: dict[tuple[str, str], deque] = defaultdict(deque)
    for path in source_paths:
        if Path(path).suffix.lower() != ".sdf":
            raise ValueError(f"Expected SDF source path, got: {path}")
        model_name = _extract_model_name_from_path(path)
        for mol in Chem.SDMolSupplier(path, removeHs=False):
            if mol is None:
                continue
            try:
                bare = Chem.RemoveHs(mol)
            except Exception:
                bare = mol
            try:
                smiles = Chem.MolToSmiles(bare)
            except Exception:
                continue
            if not smiles:
                continue
            pools[(model_name, smiles)].append(mol)

    output_sdf.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(output_sdf))
    written = 0
    try:
        for _, row in identity_df.iterrows():
            smiles = str(row[SMILES_COLUMN]).strip()
            model_name = str(row[MODEL_NAME_COLUMN]).strip()
            mol_idx = str(row[MOL_IDX_COLUMN]).strip()
            pool = pools.get((model_name, smiles))
            if not pool:
                raise ValueError(
                    "Could not rematerialize SDF molecule "
                    f"model_name={model_name!r} smiles={smiles!r}"
                )
            mol = pool.popleft()
            mol.SetProp("_Name", mol_idx)
            mol.SetProp("mol_idx", mol_idx)
            mol.SetProp("model_name", model_name)
            mol.SetProp("input_smiles", smiles)
            mol.SetProp("smiles", smiles)
            writer.write(mol)
            written += 1
    finally:
        writer.close()

    if written == 0:
        raise ValueError(f"No molecules written to docking SDF: {output_sdf}")
    return output_sdf
