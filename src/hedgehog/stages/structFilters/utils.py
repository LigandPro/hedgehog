import contextlib
import importlib
import io
import logging
import os
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

warnings.filterwarnings("ignore", category=FutureWarning)

from rdkit import Chem

# Add local Lilly binaries to PATH before importing LillyDemeritsFilters
_LILLY_BIN_PATH = (
    Path(__file__).parent
    / ".."
    / ".."
    / ".."
    / ".."
    / "modules"
    / "lilly_medchem_rules"
    / "bin"
).resolve()
if _LILLY_BIN_PATH.exists():
    current_path = os.environ.get("PATH", "")
    lilly_str = str(_LILLY_BIN_PATH)
    if lilly_str not in current_path:
        os.environ["PATH"] = f"{lilly_str}:{current_path}"

_VALID_SCHEDULERS = {"threads", "processes"}
_COMPLEXITY_FILTERS_CACHE = {}
_MOLCOMPLEXITY_ALERT_NAMES: list[str] | None = None
_ALERT_COMPILED_SMARTS: list[tuple[str, object, str]] | None = None
_ALERT_RULESET_NAMES: list[str] | None = None
_SUPPRESS_PANDASTOOLS_ENV = "HEDGEHOG_SHOW_PANDASTOOLS_WARNINGS"


def _suppress_pandastools_warning() -> None:
    if os.environ.get(_SUPPRESS_PANDASTOOLS_ENV, "").strip() == "1":
        return
    for logger_name in ("rdkit.Chem.PandasTools", "rdkit.Chem.PandasPatcher"):
        rdkit_logger = logging.getLogger(logger_name)
        rdkit_logger.setLevel(logging.ERROR)
        rdkit_logger.propagate = False
        if not any(
            isinstance(handler, logging.NullHandler)
            for handler in rdkit_logger.handlers
        ):
            rdkit_logger.addHandler(logging.NullHandler())


def _import_medchem_quietly():
    """Import medchem while muting native stdout/stderr noise."""
    stdout_fd = os.dup(1)
    stderr_fd = os.dup(2)
    devnull_fd = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull_fd, 1)
        os.dup2(devnull_fd, 2)
        return importlib.import_module("medchem")
    finally:
        os.dup2(stdout_fd, 1)
        os.dup2(stderr_fd, 2)
        os.close(devnull_fd)
        os.close(stdout_fd)
        os.close(stderr_fd)


def _resolve_scheduler(config_structFilters, scheduler_key, default="processes"):
    """Resolve scheduler name for medchem parallel APIs."""
    scheduler_raw = config_structFilters.get(scheduler_key)
    if scheduler_raw is None:
        scheduler_raw = config_structFilters.get("parallel_scheduler", default)
    scheduler = str(scheduler_raw).strip().lower()
    if scheduler not in _VALID_SCHEDULERS:
        logger.warning(
            "Invalid scheduler '%s' for '%s'. Using '%s'.",
            scheduler_raw,
            scheduler_key,
            default,
        )
        return default
    return scheduler


def _get_complexity_filter(metric_name):
    """Get cached complexity filter instance for worker process."""
    complexity_filter = _COMPLEXITY_FILTERS_CACHE.get(metric_name)
    if complexity_filter is None:
        complexity_filter = mc.complexity.ComplexityFilter(complexity_metric=metric_name)
        _COMPLEXITY_FILTERS_CACHE[metric_name] = complexity_filter
    return complexity_filter


def _compute_molcomplexity_one(args):
    """Compute all complexity metrics for one molecule."""
    mol_idx, mol = args
    alert_names = _MOLCOMPLEXITY_ALERT_NAMES or []
    row = {"_mol_idx": mol_idx}

    if mol is None:
        row.update({f"pass_{name}": False for name in alert_names})
        row["pass"] = False
        row["pass_any"] = False
        return row

    try:
        prepared_mol = Chem.RemoveHs(mol)
    except Exception:
        prepared_mol = mol

    passed_all = True
    passed_any = False

    for name in alert_names:
        complexity_filter = _get_complexity_filter(name)
        if name == "smcm":
            with (
                contextlib.redirect_stdout(io.StringIO()),
                contextlib.redirect_stderr(io.StringIO()),
            ):
                try:
                    passed = bool(complexity_filter(prepared_mol))
                except Exception:
                    passed = False
        else:
            try:
                passed = bool(complexity_filter(prepared_mol))
            except Exception:
                passed = False

        row[f"pass_{name}"] = passed
        passed_all = passed_all and passed
        passed_any = passed_any or passed

    row["pass"] = passed_all
    row["pass_any"] = passed_any
    return row


def _init_molcomplexity_worker(alert_names: list[str]) -> None:
    global _MOLCOMPLEXITY_ALERT_NAMES
    _MOLCOMPLEXITY_ALERT_NAMES = alert_names


def _split_indexed_mols(indexed_mols, n_chunks):
    """Split indexed molecules into near-equal chunks."""
    if not indexed_mols:
        return []
    n_chunks = max(1, min(n_chunks, len(indexed_mols)))
    chunk_size = (len(indexed_mols) + n_chunks - 1) // n_chunks
    return [
        indexed_mols[i : i + chunk_size]
        for i in range(0, len(indexed_mols), chunk_size)
    ]


def _process_nibr_chunk(args):
    """Run NIBR filter for one chunk with single-worker medchem call."""
    indexed_chunk, keep_details = args
    mol_indices = [item[0] for item in indexed_chunk]
    mol_chunk = [item[1] for item in indexed_chunk]

    nibr_filters = mc.structural.NIBRFilters()
    chunk_df = nibr_filters(
        mols=mol_chunk,
        n_jobs=1,
        scheduler="threads",
        keep_details=keep_details,
    )
    if len(chunk_df) != len(indexed_chunk):
        logger.warning(
            "NIBR chunk returned %d rows for %d molecules.",
            len(chunk_df),
            len(indexed_chunk),
        )
        template = chunk_df.iloc[-1].to_dict() if len(chunk_df) > 0 else None
        chunk_df = _ensure_dataframe_length(chunk_df, len(indexed_chunk), template)

    chunk_df = chunk_df.reset_index(drop=True)
    if "mol" in chunk_df.columns:
        chunk_df = chunk_df.drop(columns=["mol"])
    chunk_df["_mol_idx"] = mol_indices
    return chunk_df.to_dict("records")


_suppress_pandastools_warning()
mc = _import_medchem_quietly()

try:
    from medchem.structural.lilly_demerits import LillyDemeritsFilters

    LILLY_AVAILABLE = True
except ImportError:
    LILLY_AVAILABLE = False
    LillyDemeritsFilters = None

from hedgehog.configs.logger import load_config, logger
from hedgehog.utils.datamol_import import import_datamol_quietly
from hedgehog.utils.parallel import parallel_map, resolve_n_jobs

dm = import_datamol_quietly()

# Default columns for Lilly filter results
_LILLY_DEFAULT_COLUMNS = ["smiles", "status", "pass_filter", "demerit_score", "reasons"]


def _create_failed_row(smiles=None, reason="unsupported_or_missing"):
    """Create a failed molecule row for Lilly filter results."""
    return {
        "smiles": smiles,
        "status": "exclude",
        "pass_filter": False,
        "demerit_score": None,
        "reasons": reason,
    }


def _create_failed_dataframe(count, smiles_list=None, reason="unsupported_or_missing"):
    """Create a DataFrame of failed molecule rows."""
    if smiles_list is None:
        smiles_list = [None] * count
    return pd.DataFrame([_create_failed_row(smi, reason) for smi in smiles_list])


def _pad_dataframe_to_length(df, target_length, template_row=None):
    """Pad a DataFrame to reach target length using failed rows."""
    current_length = len(df)
    if current_length >= target_length:
        return df

    missing_count = target_length - current_length
    if template_row is not None:
        missing_rows = []
        for _ in range(missing_count):
            row = template_row.copy()
            row["smiles"] = None
            row["status"] = "exclude"
            row["pass_filter"] = False
            row["demerit_score"] = None
            if "reasons" in row:
                row["reasons"] = "unsupported_or_missing"
            missing_rows.append(row)
        missing_df = pd.DataFrame(missing_rows)
    else:
        missing_df = _create_failed_dataframe(missing_count)

    return pd.concat([df, missing_df], ignore_index=True)


def _ensure_dataframe_length(df, expected_length, template_row=None):
    """Ensure DataFrame has exactly the expected length by padding or trimming."""
    if len(df) < expected_length:
        return _pad_dataframe_to_length(df, expected_length, template_row)
    if len(df) > expected_length:
        return df.iloc[:expected_length].reset_index(drop=True)
    return df


def camelcase(any_str):
    """Convert underscore-separated string to CamelCase."""
    return "".join(word.capitalize() for word in any_str.split("_"))


def build_identity_map_from_descriptors(config):
    """Build a map of (smiles, model_name) -> mol_idx from descriptors output."""
    base_folder = Path(process_path(config["folder_to_save"]))
    id_path = base_folder / "Descriptors" / "passDescriptorsSMILES.csv"

    try:
        if id_path.exists():
            id_df = pd.read_csv(id_path)
            identity_map = {
                (row["smiles"], row["model_name"]): row["mol_idx"]
                for _, row in id_df.iterrows()
            }
            return identity_map, id_df
    except Exception:
        pass

    return {}, None


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


def sdf_to_mols(sdf_file, subsample):
    """Read molecules from SDF file with subsampling."""
    molecules = dm.read_sdf(sdf_file)
    mols_list = []
    smiles_list = []

    for i, mol in enumerate(molecules):
        if i >= subsample:
            break
        if mol is not None:
            mols_list.append(mol)
            smiles_list.append(dm.to_smiles(mol))

    return mols_list, smiles_list


def dropna(mols, smiles):
    """Remove None values from molecule and SMILES lists."""
    df = pd.DataFrame({"mols": mols, "smiles": smiles}).dropna()
    return df["mols"].tolist(), df["smiles"].tolist()


def format_number(x, pos=None):
    """Format number for display. pos parameter is for matplotlib FuncFormatter compatibility."""
    if x >= 1e6:
        return f"{x / 1e6:.1f}M"
    elif x >= 1e3:
        return f"{x / 1e3:.1f}K"
    return f"{x:.0f}"


def get_model_colors(model_names, cmap=None):
    """Generate color map for models."""
    if cmap is None:
        colors = plt.cm.YlOrRd(np.linspace(1, 0, len(model_names) + 1))
    else:
        colors = plt.colormaps.get_cmap(cmap)(np.linspace(1, 0, len(model_names) + 1))

    return dict(zip(model_names, colors, strict=False))


def clean_name(name):
    """Clean metric names for display."""
    for pattern in ["metrics", "_", ".csv"]:
        name = name.replace(pattern, "")
    return name.strip()


def filter_alerts(config):
    """Filter structural alerts based on configuration."""
    df = pd.read_csv(config["alerts_data_path"])
    mask = df["rule_set_name"].isin(config["include_rulesets"])

    for ruleset in config["exclude_descriptions"]:
        is_other_ruleset = df["rule_set_name"] != ruleset
        is_not_excluded = ~df["description"].isin(
            config["exclude_descriptions"][ruleset]
        )
        mask &= is_other_ruleset | is_not_excluded

    return df[mask]


def common_postprocessing_statistics(filter_results, res_df, stat, extend):
    """Combine filter results with statistics."""
    if stat is not None:
        res_df = pd.concat([stat, res_df])

    filter_results = filter_results.drop(columns="mol")
    if extend is not None:
        filter_extended = pd.concat([extend, filter_results], ignore_index=True)
    else:
        filter_extended = filter_results.copy()

    return res_df, filter_extended


def _apply_subsample_to_dataframe(df, subsample):
    """Apply structFilters subsampling policy to input DataFrame."""
    if subsample is None or subsample <= 0:
        return df

    if "model_name" in df.columns and df["model_name"].nunique(dropna=True) > 1:
        return df.groupby("model_name", group_keys=False).head(subsample)

    return df.head(subsample) if len(df) > subsample else df


def _parse_smiles_item(args):
    """Parse one SMILES string into RDKit Mol."""
    row_idx, smiles_raw, model_name, mol_idx = args
    try:
        mol = dm.to_mol(smiles_raw, sanitize=True)
    except Exception:
        mol = None
    return row_idx, smiles_raw, model_name, mol_idx, mol


def prepare_structfilters_input(df, subsample, parse_n_jobs, progress_cb=None):
    """Prepare structFilters payload once: subsample + SMILES parsing.

    Returns:
        dict with keys:
          - ``mols``: list of valid RDKit mols
          - ``smiles_model_mols``: tuples (smiles, model_name_or_none, mol, mol_idx)
          - ``base_df``: DataFrame with only valid rows
    """
    if df is None or len(df) == 0:
        return {
            "mols": [],
            "smiles_model_mols": [],
            "base_df": pd.DataFrame(columns=["smiles", "model_name", "mol_idx"]),
        }

    if "smiles" not in df.columns:
        raise KeyError("Input DataFrame must contain a 'smiles' column")

    data = _apply_subsample_to_dataframe(df, subsample).copy()
    if "model_name" not in data.columns:
        data["model_name"] = "single"
    if "mol_idx" not in data.columns:
        data["mol_idx"] = range(len(data))

    if len(data) == 0:
        return {
            "mols": [],
            "smiles_model_mols": [],
            "base_df": data.iloc[0:0].copy(),
        }

    is_multi = data["model_name"].nunique(dropna=True) > 1
    model_vals = data["model_name"].tolist() if is_multi else [None] * len(data)
    items = list(
        enumerate(
            zip(
                data["smiles"].tolist(),
                model_vals,
                data["mol_idx"].tolist(),
                strict=False,
            )
        )
    )
    parse_payload = [
        (row_idx, item[0], item[1], item[2]) for row_idx, item in items
    ]

    n_jobs = int(parse_n_jobs)
    if n_jobs <= 0:
        n_jobs = os.cpu_count() or 1
    parsed = parallel_map(
        _parse_smiles_item,
        parse_payload,
        n_jobs,
        progress=progress_cb,
    )

    valid_row_indices = []
    smiles_model_mols = []
    mols = []
    for row_idx, smiles_raw, model_name, mol_idx, mol in parsed:
        if mol is None:
            continue
        valid_row_indices.append(row_idx)
        smiles_model_mols.append((smiles_raw, model_name, mol, mol_idx))
        mols.append(mol)

    if valid_row_indices:
        base_df = data.iloc[valid_row_indices].copy().reset_index(drop=True)
    else:
        base_df = data.iloc[0:0].copy()

    return {
        "mols": mols,
        "smiles_model_mols": smiles_model_mols,
        "base_df": base_df,
    }


def process_prepared_payload(config, prepared_payload, apply_filter, progress_cb=None):
    """Apply one structural filter on pre-parsed payload."""
    mols = prepared_payload["mols"]
    smiles = prepared_payload["smiles_model_mols"]

    if len(mols) == 0:
        return None

    if progress_cb is not None:
        try:
            return apply_filter(config, mols, smiles, progress_cb=progress_cb)
        except TypeError:
            pass
    return apply_filter(config, mols, smiles)


def process_one_dataframe(config, df, apply_filter, subsample, progress_cb=None):
    """Process a DataFrame directly without re-reading from file.

    Args:
        config: Configuration dictionary
        df: DataFrame with 'smiles', 'model_name', and 'mol_idx' columns
        apply_filter: Filter function to apply
        subsample: Maximum number of molecules to process

    Returns:
        DataFrame with filter results or None if no molecules to process
    """
    if df is None or len(df) == 0:
        return None

    if subsample is None or subsample <= 0:
        data = df
    elif "model_name" in df.columns and df["model_name"].nunique(dropna=True) > 1:
        # In multi-model mode, apply the subsample limit per model, not globally.
        # The global sampling is handled earlier by prepare_input_data().
        data = df.groupby("model_name", group_keys=False).head(subsample)
    else:
        data = df.head(subsample) if len(df) > subsample else df
    smiles_col = "smiles"

    is_multi = data["model_name"].nunique(dropna=True) > 1
    if is_multi:
        smiles_str = data[smiles_col].tolist()
        model_names = data["model_name"].tolist()
        mols = [dm.to_mol(x, sanitize=True) for x in smiles_str]
        mol_indices = data["mol_idx"].tolist()
        smiles = list(zip(smiles_str, model_names, mols, mol_indices, strict=False))
    else:
        smiles_list = data[smiles_col].tolist()
        mols = [dm.to_mol(x, sanitize=True) for x in smiles_list]
        mol_indices = data["mol_idx"].tolist()
        smiles = list(
            zip(smiles_list, [None] * len(smiles_list), mols, mol_indices, strict=False)
        )

    # Filter out invalid molecules
    cleaned = []
    for item in smiles:
        if len(item) >= 3:
            smi, model, mol = item[0], item[1], item[2]
            if mol is not None:
                if len(item) >= 4:
                    mol_idx = item[3]
                    cleaned.append((smi, model, mol, mol_idx))
                else:
                    cleaned.append((smi, model, mol))
    smiles = cleaned
    mols = [it[2] for it in smiles]

    if len(mols) == 0:
        return None

    if progress_cb is not None:
        try:
            return apply_filter(config, mols, smiles, progress_cb=progress_cb)
        except TypeError:
            pass
    return apply_filter(config, mols, smiles)


def process_one_file(config, input_path, apply_filter, subsample, progress_cb=None):
    """Process molecules from file through a filter function.

    For CSV files, consider using process_one_dataframe() directly if you
    already have the DataFrame loaded to avoid reading the file twice.

    Args:
        config: Configuration dictionary
        input_path: Path to input file (csv, smi, sdf, or txt)
        apply_filter: Filter function to apply
        subsample: Maximum number of molecules to process

    Returns:
        DataFrame with filter results or None if no molecules to process
    """
    input_type = input_path[input_path.rfind(".") + 1 :]
    assert input_type in {"csv", "smi", "sdf", "txt"}

    if input_type == "csv":
        data = pd.read_csv(input_path)
        return process_one_dataframe(config, data, apply_filter, subsample)

    elif input_type == "smi" or input_type == "txt":
        with open(input_path) as file:
            lines = [line.rstrip("\n") for line in file]

        if subsample <= len(lines):
            lines = np.random.permutation(lines)[:subsample].tolist()

        smiles = []
        model_names = []
        for line in lines:
            parts = line.split(",")
            if len(parts) == 2:
                smi, model = parts
                smiles.append(smi)
                model_names.append(model)
            else:
                smiles.append(parts[0])
        mols = [dm.to_mol(x) for x in smiles]
        if len(model_names) == len(smiles):
            smiles = list(zip(smiles, model_names, mols, strict=False))

    elif input_type == "sdf":
        mols, smiles = sdf_to_mols(input_path, subsample)

    if isinstance(smiles[0], tuple):
        cleaned = []
        for item in smiles:
            if len(item) >= 3:
                smi, model, mol = item[0], item[1], item[2]
                if mol is not None:
                    if len(item) >= 4:
                        mol_idx = item[3]
                        cleaned.append((smi, model, mol, mol_idx))
                    else:
                        cleaned.append((smi, model, mol))
        smiles = cleaned
        mols = [it[2] for it in smiles]
    else:
        mols, smiles = dropna(mols, smiles)

    assert len(mols) == len(smiles), f"{len(mols)}, {len(smiles)}"
    if not isinstance(smiles[0], tuple):
        assert len(mols) <= subsample

    for mol, smi in zip(mols, smiles, strict=False):
        smi_val = smi[0] if isinstance(smi, tuple) else smi
        assert mol is not None, f"{smi_val}"

    final_result = None
    if len(mols) > 0:
        if isinstance(smiles[0], tuple):
            if progress_cb is not None:
                try:
                    final_result = apply_filter(
                        config, mols, smiles, progress_cb=progress_cb
                    )
                except TypeError:
                    final_result = apply_filter(config, mols, smiles)
            else:
                final_result = apply_filter(config, mols, smiles)
        else:
            if progress_cb is not None:
                try:
                    final_result = apply_filter(config, mols, progress_cb=progress_cb)
                except TypeError:
                    final_result = apply_filter(config, mols)
            else:
                final_result = apply_filter(config, mols)

    return final_result


def add_model_name_col(final_result, smiles_with_model):
    """Add smiles, model_name, and mol_idx columns from input data."""
    expected_len = len(smiles_with_model)
    actual_len = len(final_result)

    if actual_len != expected_len:
        logger.error(
            "Length mismatch in add_model_name_col: final_result has %d rows, "
            "but smiles_with_model has %d items.",
            actual_len,
            expected_len,
        )

        if actual_len < expected_len:
            logger.warning(
                "Padding final_result from %d to %d rows.", actual_len, expected_len
            )
            if len(final_result) > 0:
                template_row = final_result.iloc[-1].to_dict()
                final_result = _pad_dataframe_to_length(
                    final_result, expected_len, template_row
                )
                logger.info("Padded to %d rows.", len(final_result))
            else:
                logger.error("final_result is empty. Creating all rows from scratch.")
                smiles_list = [
                    item[0] if isinstance(item, tuple) and len(item) > 0 else None
                    for item in smiles_with_model
                ]
                final_result = _create_failed_dataframe(expected_len, smiles_list)
        else:
            logger.warning(
                "Trimming final_result from %d to %d rows.", actual_len, expected_len
            )
            final_result = final_result.iloc[:expected_len].reset_index(drop=True)

    if len(final_result) != expected_len:
        msg = (
            f"CRITICAL: After padding/trimming, final_result length ({len(final_result)}) "
            f"still doesn't match expected ({expected_len})."
        )
        raise ValueError(msg)

    smiles_vals = [item[0] for item in smiles_with_model]
    model_vals = [
        item[1] if item[1] is not None else "single" for item in smiles_with_model
    ]
    mol_idx_vals = [item[3] if len(item) >= 4 else None for item in smiles_with_model]

    final_result["smiles"] = smiles_vals
    final_result["model_name"] = model_vals
    final_result["mol_idx"] = mol_idx_vals

    return final_result


def filter_function_applier(filter_name):
    filters = {
        "common_alerts": apply_structural_alerts,
        "molgraph_stats": apply_molgraph_stats,
        "molcomplexity": apply_molcomplexity_filters,
        "NIBR": apply_nibr_filter,
        "bredt": apply_bredt_filter,
        "lilly": apply_lilly_filter,
        "protecting_groups": apply_protecting_groups,
        "ring_infraction": apply_ring_infraction,
        "stereo_center": apply_stereo_center,
        "halogenicity": apply_halogenicity,
        "symmetry": apply_symmetry,
    }
    if filter_name not in filters:
        raise ValueError(f"Filter {filter_name} not found")
    return filters[filter_name]


def _check_alerts_single_mol(args):
    """Check all alert SMARTS patterns against a single molecule.

    Args:
        args: Tuple of (mol_idx, mol).

    Returns:
        dict with per-rule-set pass/fail flags and failure reasons.
    """
    mol_idx, mol = args
    compiled_smarts = _ALERT_COMPILED_SMARTS or []
    rule_set_names = _ALERT_RULESET_NAMES or []

    row_data = {"_mol_idx": mol_idx}
    row_data["smiles"] = dm.to_smiles(mol) if mol is not None else None

    if mol is not None:
        reasons_map: dict[str, list[str]] = {}
        for name in rule_set_names:
            row_data[f"pass_{name}"] = True
            reasons_map[name] = []

        for ruleset, patt, desc in compiled_smarts:
            if mol.HasSubstructMatch(patt):
                row_data[f"pass_{ruleset}"] = False
                if desc and desc not in reasons_map[ruleset]:
                    reasons_map[ruleset].append(desc)

        for name in rule_set_names:
            row_data[f"reasons_{name}"] = ";".join(reasons_map[name])
    else:
        for name in rule_set_names:
            row_data[f"pass_{name}"] = False
            row_data[f"reasons_{name}"] = "invalid_molecule"

    return row_data


def _init_alert_worker(
    compiled_smarts: list[tuple[str, object, str]],
    rule_set_names: list[str],
) -> None:
    global _ALERT_COMPILED_SMARTS, _ALERT_RULESET_NAMES
    _ALERT_COMPILED_SMARTS = compiled_smarts
    _ALERT_RULESET_NAMES = rule_set_names


def apply_structural_alerts(config, mols, smiles_modelName_mols=None, progress_cb=None):
    logger.info("Calculating Common Alerts...")

    # Load config and filter alerts ONCE outside
    config_structFilters = load_config(config["config_structFilters"])
    alert_data = filter_alerts(config_structFilters)

    # Get rule set names for column creation
    rule_set_names = alert_data["rule_set_name"].unique().tolist()

    logger.info(
        "Processing %d filtered molecules with %d alert rule sets",
        len(mols),
        len(rule_set_names),
    )

    # Pre-compile SMARTS patterns once to avoid recompiling per molecule.
    compiled_smarts: list[tuple[str, Chem.Mol, str]] = []
    for _, row in alert_data.iterrows():
        smarts = row.get("smarts")
        ruleset = row.get("rule_set_name")
        desc = row.get("description") or ""
        if not isinstance(smarts, str) or not isinstance(ruleset, str):
            continue
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue
        compiled_smarts.append((ruleset, patt, str(desc)))

    # Parallel per-molecule alert checking.
    n_jobs = resolve_n_jobs(config_structFilters, config)
    logger.info("Common Alerts workers: %d (scheduler=processes)", n_jobs)
    items = [(i, mol) for i, mol in enumerate(mols)]
    results_list = parallel_map(
        _check_alerts_single_mol,
        items,
        n_jobs,
        progress=progress_cb,
        initializer=_init_alert_worker,
        initargs=(compiled_smarts, rule_set_names),
    )

    # Restore mol column (RDKit Mol objects are not picklable across processes).
    for row_data in results_list:
        idx = row_data.pop("_mol_idx")
        row_data["mol"] = mols[idx]

    results = pd.DataFrame(results_list)

    # Calculate pass and pass_any columns.
    pass_cols = [
        c
        for c in results.columns
        if c.startswith("pass_") and c not in {"pass", "pass_any"}
    ]
    results["pass"] = results[pass_cols].all(axis=1)
    results["pass_any"] = results[pass_cols].any(axis=1)

    if smiles_modelName_mols is not None:
        results = add_model_name_col(results, smiles_modelName_mols)
    return results


def apply_molgraph_stats(config, mols, smiles_modelName_mols=None):
    logger.info("Calculating Molecular Graph statistics...")
    severities = list(range(1, 12))
    config_structFilters = load_config(config["config_structFilters"])
    n_jobs = resolve_n_jobs(config_structFilters, config)
    scheduler = _resolve_scheduler(config_structFilters, "molgraph_scheduler")
    logger.info("MolGraph workers: %d (scheduler=%s)", n_jobs, scheduler)

    results = {"mol": mols}
    for s in severities:
        out = mc.functional.molecular_graph_filter(
            mols=mols,
            max_severity=s,
            n_jobs=n_jobs,
            scheduler=scheduler,
            progress=False,
            return_idx=False,
        )
        results[f"pass_{s}"] = out
    results = pd.DataFrame(results)

    if smiles_modelName_mols is not None:
        results = add_model_name_col(results, smiles_modelName_mols)
    return results


def apply_molcomplexity_filters(config, mols, smiles_modelName_mols=None):
    logger.info("Calculating Complexity filters...")
    config_path = config.get("config_structFilters")
    config_structFilters = load_config(config_path) if config_path else {}
    n_jobs = resolve_n_jobs(config_structFilters, config)
    logger.info("MolComplexity workers: %d (scheduler=processes)", n_jobs)

    alert_names = mc.complexity.ComplexityFilter.list_default_available_filters()
    items = [(i, mol) for i, mol in enumerate(mols)]
    rows = parallel_map(
        _compute_molcomplexity_one,
        items,
        n_jobs,
        initializer=_init_molcomplexity_worker,
        initargs=(alert_names,),
    )
    for row in rows:
        mol_idx = row.pop("_mol_idx")
        row["mol"] = mols[mol_idx]
    final_result = pd.DataFrame(rows)

    if smiles_modelName_mols is not None:
        final_result = add_model_name_col(final_result, smiles_modelName_mols)
    return final_result


def apply_bredt_filter(config, mols, smiles_modelName_mols=None):
    logger.info("Calculating Bredt filter...")
    config_structFilters = load_config(config["config_structFilters"])
    n_jobs = resolve_n_jobs(config_structFilters, config)
    scheduler = _resolve_scheduler(config_structFilters, "bredt_scheduler")
    logger.info("Bredt workers: %d (scheduler=%s)", n_jobs, scheduler)
    out = mc.functional.bredt_filter(
        mols=mols,
        n_jobs=n_jobs,
        scheduler=scheduler,
        progress=False,
        return_idx=False,
    )
    results = pd.DataFrame({"mol": mols, "pass": out})
    if smiles_modelName_mols is not None:
        results = add_model_name_col(results, smiles_modelName_mols)
    return results


def apply_protecting_groups(config, mols, smiles_modelName_mols=None):
    logger.info("Calculating Protecting Groups filter...")
    config_structFilters = load_config(config["config_structFilters"])
    n_jobs = resolve_n_jobs(config_structFilters, config)
    scheduler = _resolve_scheduler(config_structFilters, "protecting_groups_scheduler")
    logger.info("Protecting Groups workers: %d (scheduler=%s)", n_jobs, scheduler)
    out = mc.functional.protecting_groups_filter(
        mols=mols,
        n_jobs=n_jobs,
        scheduler=scheduler,
        progress=False,
        return_idx=False,
    )
    results = pd.DataFrame({"mol": mols, "pass": out})
    if smiles_modelName_mols is not None:
        results = add_model_name_col(results, smiles_modelName_mols)
    return results


def apply_ring_infraction(config, mols, smiles_modelName_mols=None):
    logger.info("Calculating Ring Infraction filter...")
    config_sf = load_config(config["config_structFilters"])
    n_jobs = resolve_n_jobs(config_sf, config)
    scheduler = _resolve_scheduler(config_sf, "ring_infraction_scheduler")
    logger.info("Ring Infraction workers: %d (scheduler=%s)", n_jobs, scheduler)
    min_size = config_sf.get("ring_infraction_hetcycle_min_size", 4)
    out = mc.functional.ring_infraction_filter(
        mols=mols,
        hetcycle_min_size=min_size,
        n_jobs=n_jobs,
        scheduler=scheduler,
        progress=False,
        return_idx=False,
    )
    results = pd.DataFrame({"mol": mols, "pass": out})
    if smiles_modelName_mols is not None:
        results = add_model_name_col(results, smiles_modelName_mols)
    return results


def apply_stereo_center(config, mols, smiles_modelName_mols=None):
    logger.info("Calculating Stereo Center filter...")
    config_sf = load_config(config["config_structFilters"])
    n_jobs = resolve_n_jobs(config_sf, config)
    scheduler = _resolve_scheduler(config_sf, "stereo_center_scheduler")
    logger.info("Stereo Center workers: %d (scheduler=%s)", n_jobs, scheduler)
    max_centers = config_sf.get("stereo_max_centers", 4)
    max_undefined = config_sf.get("stereo_max_undefined", 2)
    out = mc.functional.num_stereo_center_filter(
        mols=mols,
        max_stereo_centers=max_centers,
        max_undefined_stereo_centers=max_undefined,
        n_jobs=n_jobs,
        scheduler=scheduler,
        progress=False,
        return_idx=False,
    )
    results = pd.DataFrame({"mol": mols, "pass": out})
    if smiles_modelName_mols is not None:
        results = add_model_name_col(results, smiles_modelName_mols)
    return results


def apply_halogenicity(config, mols, smiles_modelName_mols=None):
    logger.info("Calculating Halogenicity filter...")
    config_sf = load_config(config["config_structFilters"])
    n_jobs = resolve_n_jobs(config_sf, config)
    scheduler = _resolve_scheduler(config_sf, "halogenicity_scheduler")
    logger.info("Halogenicity workers: %d (scheduler=%s)", n_jobs, scheduler)
    out = mc.functional.halogenicity_filter(
        mols=mols,
        thresh_F=config_sf.get("halogenicity_thresh_F", 6),
        thresh_Br=config_sf.get("halogenicity_thresh_Br", 3),
        thresh_Cl=config_sf.get("halogenicity_thresh_Cl", 3),
        n_jobs=n_jobs,
        scheduler=scheduler,
        progress=False,
        return_idx=False,
    )
    results = pd.DataFrame({"mol": mols, "pass": out})
    if smiles_modelName_mols is not None:
        results = add_model_name_col(results, smiles_modelName_mols)
    return results


def apply_symmetry(config, mols, smiles_modelName_mols=None):
    logger.info("Calculating Symmetry filter...")
    config_sf = load_config(config["config_structFilters"])
    n_jobs = resolve_n_jobs(config_sf, config)
    scheduler = _resolve_scheduler(config_sf, "symmetry_scheduler")
    logger.info("Symmetry workers: %d (scheduler=%s)", n_jobs, scheduler)
    threshold = config_sf.get("symmetry_threshold", 0.8)
    out = mc.functional.symmetry_filter(
        mols=mols,
        symmetry_threshold=threshold,
        n_jobs=n_jobs,
        scheduler=scheduler,
        progress=False,
        return_idx=False,
    )
    results = pd.DataFrame({"mol": mols, "pass": out})
    if smiles_modelName_mols is not None:
        results = add_model_name_col(results, smiles_modelName_mols)
    return results


def apply_nibr_filter(config, mols, smiles_modelName_mols=None):
    logger.info("Calculating NIBR filter...")
    config_structFilters = load_config(config["config_structFilters"])
    n_jobs = resolve_n_jobs(config_structFilters, config)
    scheduler = _resolve_scheduler(config_structFilters, "nibr_scheduler")
    logger.info("NIBR workers: %d (scheduler=%s)", n_jobs, scheduler)

    indexed_mols = list(enumerate(mols))
    n_workers = max(1, min(n_jobs, len(indexed_mols)))
    chunked_mols = _split_indexed_mols(indexed_mols, n_workers)
    chunk_payloads = [(chunk, True) for chunk in chunked_mols]
    chunk_results = parallel_map(
        _process_nibr_chunk,
        chunk_payloads,
        n_workers,
        chunksize=1,
    )

    rows = []
    for chunk_rows in chunk_results:
        rows.extend(chunk_rows)
    results = pd.DataFrame(rows)

    expected_indices = set(range(len(mols)))
    observed_indices = set(results["_mol_idx"].tolist()) if "_mol_idx" in results else set()
    if len(results) != len(mols) or observed_indices != expected_indices:
        logger.warning(
            "NIBR chunked processing mismatch (got %d/%d rows). Falling back to medchem native parallel call.",
            len(results),
            len(mols),
        )
        nibr_filters = mc.structural.NIBRFilters()
        results = nibr_filters(
            mols=mols,
            n_jobs=n_jobs,
            scheduler=scheduler,
            keep_details=True,
        )
    else:
        results = results.sort_values("_mol_idx").reset_index(drop=True)
        results["mol"] = [mols[i] for i in results["_mol_idx"]]
        results = results.drop(columns=["_mol_idx"])

    if smiles_modelName_mols is not None:
        results = add_model_name_col(results, smiles_modelName_mols)
    return results


def _process_lilly_batch(dfilter, batch, n_jobs, scheduler):
    """Process a single batch of molecules through Lilly filter."""
    try:
        batch_result = dfilter(mols=batch, n_jobs=n_jobs, scheduler=scheduler)
        if len(batch_result) != len(batch):
            logger.warning(
                "Lilly returned %d results for %d molecules.",
                len(batch_result),
                len(batch),
            )
            template = (
                batch_result.iloc[-1].to_dict() if len(batch_result) > 0 else None
            )
            batch_result = _ensure_dataframe_length(batch_result, len(batch), template)
        return batch_result
    except Exception as batch_error:
        if "Length of values" in str(
            batch_error
        ) or "does not match length of index" in str(batch_error):
            return _process_lilly_one_by_one(dfilter, batch, scheduler)
        smiles_list = [dm.to_smiles(m) if m else None for m in batch]
        return _create_failed_dataframe(
            len(batch), smiles_list, "batch_processing_failed"
        )


def _process_lilly_one_by_one(dfilter, batch, scheduler):
    """Process molecules one by one as fallback."""
    one_by_one_results = []
    for mol in batch:
        try:
            single_result = dfilter(mols=[mol], n_jobs=1, scheduler=scheduler)
            if len(single_result) > 0:
                one_by_one_results.append(single_result.iloc[0].to_dict())
            else:
                smi = dm.to_smiles(mol) if mol else None
                one_by_one_results.append(
                    _create_failed_row(smi, "unsupported_or_missing")
                )
        except Exception:
            smi = dm.to_smiles(mol) if mol else None
            one_by_one_results.append(_create_failed_row(smi, "processing_failed"))
    return pd.DataFrame(one_by_one_results)


def _run_lilly_in_batches(dfilter, valid_mols, n_jobs, scheduler, batch_size=500):
    """Run Lilly filter in batches with fallback to one-by-one processing."""
    batch_results = []
    for i in range(0, len(valid_mols), batch_size):
        batch = valid_mols[i : i + batch_size]
        batch_result = _process_lilly_batch(dfilter, batch, n_jobs, scheduler)
        batch_results.append(batch_result)

    if not batch_results:
        return None
    return pd.concat(batch_results, ignore_index=True)


def _process_lilly_chunk(args):
    """Process one Lilly chunk in an isolated worker process."""
    indexed_chunk, scheduler = args
    mol_indices = [item[0] for item in indexed_chunk]
    mol_chunk = [item[1] for item in indexed_chunk]

    dfilter = LillyDemeritsFilters()
    try:
        chunk_result = dfilter(mols=mol_chunk, n_jobs=1, scheduler=scheduler)
    except ValueError as error:
        if "Length of values" in str(error) or "does not match length of index" in str(
            error
        ):
            chunk_result = _process_lilly_one_by_one(dfilter, mol_chunk, scheduler)
        else:
            raise

    if len(chunk_result) != len(indexed_chunk):
        template = chunk_result.iloc[-1].to_dict() if len(chunk_result) > 0 else None
        chunk_result = _ensure_dataframe_length(chunk_result, len(indexed_chunk), template)

    chunk_result = chunk_result.reset_index(drop=True)
    if "mol" in chunk_result.columns:
        chunk_result = chunk_result.drop(columns=["mol"])
    chunk_result["_mol_idx"] = mol_indices
    return chunk_result.to_dict("records")


def _reconstruct_full_results(results, valid_indices, expected_len, input_smiles):
    """Reconstruct full results DataFrame including invalid molecules."""
    complete_results = []
    valid_idx = 0
    valid_indices_set = set(valid_indices)

    for orig_idx in range(expected_len):
        smi = input_smiles[orig_idx] if orig_idx < len(input_smiles) else None
        if orig_idx in valid_indices_set:
            if valid_idx < len(results):
                complete_results.append(results.iloc[valid_idx].to_dict())
            else:
                complete_results.append(
                    _create_failed_row(smi, "unsupported_or_missing")
                )
            valid_idx += 1
        else:
            complete_results.append(_create_failed_row(smi, "invalid_molecule"))

    return pd.DataFrame(complete_results)


def apply_lilly_filter(config, mols, smiles_modelName_mols=None):
    """Apply Lilly demerits filter to molecules."""
    if not LILLY_AVAILABLE:
        raise ImportError(
            "Lilly demerits filter is not available. "
            "This filter requires conda/mamba-installed binaries. "
            "Install with: conda install lilly-medchem-rules\n"
            "Or disable this filter by setting 'calculate_lilly: False' in config_structFilters.yml"
        )

    logger.info("Calculating Lilly filter...")
    config_structFilters = load_config(config["config_structFilters"])
    n_jobs = resolve_n_jobs(config_structFilters, config)
    scheduler = _resolve_scheduler(config_structFilters, "lilly_scheduler")
    if scheduler != "threads":
        logger.warning(
            "Lilly supports only threads scheduler. Falling back to 'threads' (got '%s').",
            scheduler,
        )
        scheduler = "threads"
    logger.info("Lilly workers: %d (scheduler=%s)", n_jobs, scheduler)

    if smiles_modelName_mols is not None:
        expected_len = len(smiles_modelName_mols)
        input_smiles = [
            item[0] if isinstance(item, tuple) and len(item) > 0 else None
            for item in smiles_modelName_mols
        ]
    else:
        expected_len = len(mols)
        input_smiles = [dm.to_smiles(mol) if mol is not None else None for mol in mols]

    # Collect valid molecules with original indices
    valid_mols = []
    valid_indexed = []
    valid_indices = []
    for idx, mol in enumerate(mols):
        if mol is not None:
            try:
                smi = dm.to_smiles(mol)
                if smi:
                    valid_mols.append(mol)
                    valid_indexed.append((idx, mol))
                    valid_indices.append(idx)
            except Exception:
                pass

    if not valid_mols:
        results = _create_failed_dataframe(
            expected_len, input_smiles, "invalid_molecule"
        )
        if smiles_modelName_mols is not None:
            results = add_model_name_col(results, smiles_modelName_mols)
        return results

    # Run Lilly across worker processes, one thread-backed call per worker.
    n_workers = max(1, min(n_jobs, len(valid_indexed)))
    chunked_mols = _split_indexed_mols(valid_indexed, n_workers)
    chunk_payloads = [(chunk, scheduler) for chunk in chunked_mols]

    try:
        chunk_results = parallel_map(
            _process_lilly_chunk,
            chunk_payloads,
            n_workers,
            chunksize=1,
        )
        rows = []
        for chunk_rows in chunk_results:
            rows.extend(chunk_rows)
        results = pd.DataFrame(rows)
    except Exception as error:
        logger.warning(
            "Parallel Lilly execution failed (%s). Falling back to batched native mode.",
            error,
        )
        dfilter = LillyDemeritsFilters()
        results = _run_lilly_in_batches(dfilter, valid_mols, n_jobs, scheduler)
        if results is None:
            raise ValueError("All Lilly batches failed in fallback mode.") from error
        template = results.iloc[-1].to_dict() if len(results) > 0 else None
        results = _ensure_dataframe_length(results, len(valid_mols), template)
    else:
        expected_indices = set(valid_indices)
        observed_indices = (
            set(results["_mol_idx"].tolist()) if "_mol_idx" in results.columns else set()
        )
        if len(results) != len(valid_mols) or observed_indices != expected_indices:
            logger.warning(
                "Lilly chunked processing mismatch (got %d/%d rows). Falling back to batched native mode.",
                len(results),
                len(valid_mols),
            )
            dfilter = LillyDemeritsFilters()
            results = _run_lilly_in_batches(dfilter, valid_mols, n_jobs, scheduler)
            if results is None:
                raise ValueError("All Lilly batches failed in fallback mode.")
            template = results.iloc[-1].to_dict() if len(results) > 0 else None
            results = _ensure_dataframe_length(results, len(valid_mols), template)
        else:
            results = results.sort_values("_mol_idx").reset_index(drop=True)
            results = results.drop(columns=["_mol_idx"])

    # Reconstruct full results including invalid molecules
    results = _reconstruct_full_results(
        results, valid_indices, expected_len, input_smiles
    )
    results["mol"] = [mols[i] if i < len(mols) else None for i in range(expected_len)]

    # Final length check
    results = _ensure_dataframe_length(results, expected_len)
    if len(results) != expected_len:
        raise ValueError(
            f"CRITICAL: Results length ({len(results)}) doesn't match expected ({expected_len})"
        )

    if smiles_modelName_mols is not None:
        results = add_model_name_col(results, smiles_modelName_mols)
    return results


def _get_pass_column(df, fallback_column=None, fallback_condition=None):
    """Determine the pass column name, creating it if necessary.

    Args:
        df: DataFrame to check
        fallback_column: Column name to use for creating pass column if not found
        fallback_condition: Lambda function to apply to fallback_column

    Returns:
        Name of the pass column
    """
    if "pass" in df.columns:
        return "pass"
    if "pass_filter" in df.columns:
        return "pass_filter"
    if fallback_column and fallback_column in df.columns and fallback_condition:
        df["pass"] = fallback_condition(df[fallback_column])
        return "pass"
    return None


def _ensure_pass_column_in_extended(
    filter_extended, pass_col, filter_results, fallback_column=None
):
    """Ensure the 'pass' column exists in filter_extended DataFrame."""
    if "pass" in filter_extended.columns:
        return filter_extended

    if pass_col and pass_col in filter_extended.columns and pass_col != "pass":
        filter_extended.rename(columns={pass_col: "pass"}, inplace=True)
    elif "pass_filter" in filter_extended.columns:
        filter_extended.rename(columns={"pass_filter": "pass"}, inplace=True)
    elif fallback_column and fallback_column in filter_extended.columns:
        filter_extended["pass"] = filter_extended[fallback_column] == 0
    elif "pass" in filter_results.columns and len(filter_extended) == len(
        filter_results
    ):
        filter_extended["pass"] = filter_results["pass"].values

    return filter_extended


def _create_base_stats_df(model_name, num_mol, **extra_columns):
    """Create base statistics DataFrame with model_name and num_mol."""
    data = {"model_name": [model_name], "num_mol": [num_mol]}
    data.update({k: [v] for k, v in extra_columns.items()})
    return pd.DataFrame(data)


def get_basic_stats(
    config, filter_results, model_name, filter_name, stat=None, extend=None
):
    """Calculate basic statistics for filter results."""
    model_col = (
        filter_results["model_name"] if "model_name" in filter_results.columns else None
    )
    is_multi = (not isinstance(model_name, str)) or (
        model_col is not None and model_col.nunique(dropna=True) > 1
    )
    if is_multi:
        if model_col is None:
            raise ValueError(
                "Multi-model statistics requested but 'model_name' column is missing"
            )
        all_res = []
        all_extended = []
        for model, group in filter_results.groupby("model_name"):
            res_df, filter_extended = get_basic_stats(
                config, group.copy(), model, filter_name, stat, extend
            )
            all_res.append(res_df)
            all_extended.append(filter_extended)
        return pd.concat(all_res, ignore_index=True), pd.concat(
            all_extended, ignore_index=True
        )

    num_mol = len(filter_results)
    filter_results.dropna(subset="mol", inplace=True)
    if isinstance(model_name, str):
        filter_results["model_name"] = model_name

    if filter_name == "common_alerts":
        res_df = _create_base_stats_df(
            model_name,
            num_mol,
            all_banned_ratio=(~filter_results["pass"]).mean(),
            any_banned_ratio=(~filter_results["pass_any"]).mean(),
        )
        for name in config["include_rulesets"]:
            res_df[f"{name}_banned_ratio"] = 1 - filter_results[f"pass_{name}"].mean()
        return common_postprocessing_statistics(filter_results, res_df, stat, extend)

    if filter_name == "molgraph_stats":
        res_df = _create_base_stats_df(model_name, num_mol)
        for i in range(1, 12):
            res_df[f"banned_ratio_s_{i}"] = 1 - filter_results[f"pass_{i}"].mean()
        res_df, filter_extended = common_postprocessing_statistics(
            filter_results, res_df, stat, extend
        )
        pass_cols = [
            col
            for col in filter_extended.columns
            if col.startswith("pass_") and col != "pass_any"
        ]
        filter_extended["pass"] = filter_extended[pass_cols].all(axis=1)
        return res_df, filter_extended

    if filter_name == "molcomplexity":
        res_df = _create_base_stats_df(
            model_name,
            num_mol,
            all_banned_ratio=1 - filter_results["pass"].mean(),
            any_banned_ratio=1 - filter_results["pass_any"].mean(),
        )
        for name in mc.complexity.ComplexityFilter.list_default_available_filters():
            res_df[f"{name}_banned_ratio"] = 1 - filter_results[f"pass_{name}"].mean()
        return common_postprocessing_statistics(filter_results, res_df, stat, extend)

    if filter_name in (
        "bredt",
        "protecting_groups",
        "ring_infraction",
        "stereo_center",
        "halogenicity",
        "symmetry",
    ):
        res_df = _create_base_stats_df(model_name, num_mol)
        res_df["banned_ratio"] = 1 - filter_results["pass"].mean()
        return common_postprocessing_statistics(filter_results, res_df, stat, extend)

    if filter_name == "NIBR":
        res_df = _create_base_stats_df(
            model_name,
            num_mol,
            mean_severity=filter_results.severity.mean(),
            max_severity=filter_results.severity.max(),
            mean_n_covalent_motif=filter_results.n_covalent_motif.mean(),
            mean_nonzero_special_mol=(filter_results.special_mol > 0).mean(),
        )
        pass_col = _get_pass_column(filter_results, "severity", lambda x: x == 0)
        res_df["banned_ratio"] = 1 - filter_results[pass_col].mean()
        res_df, filter_extended = common_postprocessing_statistics(
            filter_results, res_df, stat, extend
        )
        filter_extended = _ensure_pass_column_in_extended(
            filter_extended, pass_col, filter_results, "severity"
        )
        return res_df, filter_extended

    if filter_name == "lilly":
        res_df = _create_base_stats_df(
            model_name,
            num_mol,
            mean_noNA_demerit_score=filter_results.demerit_score.dropna().mean(),
        )
        pass_col = _get_pass_column(filter_results, "demerit_score", lambda x: x == 0)
        res_df["banned_ratio"] = 1 - filter_results[pass_col].mean()
        res_df, filter_extended = common_postprocessing_statistics(
            filter_results, res_df, stat, extend
        )
        filter_extended = _ensure_pass_column_in_extended(
            filter_extended, pass_col, filter_results, "demerit_score"
        )
        return res_df, filter_extended

    raise ValueError(f"Filter {filter_name} not found")


def check_paths(config, paths):
    all_filters = {}
    for k, v in config.items():
        if "calculate_" in k:
            k = k.replace("calculate_", "")
            all_filters[k] = v

    path_folders = []
    for path in paths:
        parts = path.split("/")
        if len(parts) >= 2:
            folder_name = parts[-2]
            path_folders.append(folder_name.lower())

    missing_filters = []
    for filter_name, enabled in all_filters.items():
        if enabled:
            filter_name_lower = filter_name.lower()
            filter_name_no_underscore = filter_name_lower.replace("_", "")
            found = any(
                filter_name_lower == folder
                or filter_name_no_underscore == folder.replace("_", "")
                or filter_name_lower in folder
                or filter_name_no_underscore in folder.replace("_", "")
                for folder in path_folders
            )
            if not found:
                missing_filters.append(filter_name_no_underscore)

    if len(missing_filters) > 0:
        raise AssertionError(
            f"Invalid filter name(s) missing: {', '.join(missing_filters)}"
        )
    return True


def plot_calculated_stats(config, stage_dir):
    """Plot calculated statistics for structural filters."""
    folder_to_save = Path(process_path(config["folder_to_save"]))
    config_structFilters = load_config(config["config_structFilters"])

    struct_folder = folder_to_save / stage_dir
    paths = list(struct_folder.glob("*/metrics.csv"))
    if not paths:
        paths = list(struct_folder.glob("*metrics.csv"))
    paths = [str(p) for p in paths]
    check_paths(config_structFilters, paths)

    datas = []
    filter_names = []
    all_model_names = set()

    for path in paths:
        data = pd.read_csv(path)

        all_model_names.update(data["model_name"].dropna().unique())
        data.set_index("model_name", inplace=True)

        banned_cols = [col for col in data.columns if "banned_ratio" in col]
        data_filtered = data[banned_cols + ["num_mol"]].copy()
        for banned_col in banned_cols:
            data_filtered.loc[:, f"num_banned_{banned_col}"] = (
                data_filtered[banned_col] * data_filtered["num_mol"]
            )
        datas.append(data_filtered)

        filter_name = path.split("/")[-1].replace("_metrics.csv", "")
        filter_names.append(filter_name)

    model_name_set = sorted(list(all_model_names))

    filter_results = {}
    filters_to_find = list(struct_folder.glob("*/filtered_molecules.csv"))
    if not filters_to_find:
        filters_to_find = list(struct_folder.glob("*filteredMols.csv"))
    filters_to_find = [str(p) for p in filters_to_find]

    for path in filters_to_find:
        try:
            filter_data = pd.read_csv(path)
            filter_name = path.split("/")[-1].split("filteredMols.csv")[0].strip("_")

            num_passed_by_model = None
            if "pass" in filter_data.columns:
                passed = filter_data[filter_data["pass"]]
                if len(passed) > 0:
                    num_passed_by_model = passed.groupby("model_name").size().to_dict()

            if num_passed_by_model is not None:
                filter_results[filter_name] = num_passed_by_model
            else:
                default_models = filter_data["model_name"].unique().tolist()
                filter_results[filter_name] = {m: 0 for m in default_models}

        except (IndexError, FileNotFoundError) as e:
            logger.warning("Could not process %s: %s", path, e)
            filter_results[filter_name] = {}

    all_models = model_name_set

    for filter_name, values in filter_results.items():
        if len(values) != len(all_models):
            for model in all_models:
                if model not in values:
                    filter_results[filter_name][model] = 0
        filter_results[filter_name] = dict(sorted(filter_results[filter_name].items()))

    n_plots = len(datas)
    n_cols = 2
    n_rows = (n_plots + n_cols - 1) // n_cols

    plt.figure(figsize=(40, 5 * n_rows))
    for idx, (data, filter_name) in enumerate(zip(datas, filter_names, strict=False)):
        ax = plt.subplot(n_rows, n_cols, idx + 1)
        models = data.index
        x = np.arange(len(models))
        width = 0.8
        total_mols = data["num_mol"].sum()
        total = ax.barh(
            x,
            data.loc[models, "num_mol"],
            width,
            label=f"Total Molecules ({format_number(total_mols)})",
            color="#E5E5E5",
            alpha=0.5,
        )

        clean_filter_name = filter_name.split("/")[-1].lower()
        for known_filter in filter_results.keys():
            if known_filter.lower() in clean_filter_name:
                for i, (model, passed) in enumerate(
                    filter_results[known_filter].items()
                ):
                    bar_center_x = data.loc[models, "num_mol"].values[0] / 2
                    bar_center_y = x[i]
                    model_total = (
                        data.loc[model, "num_mol"]
                        if model in data.index
                        else data["num_mol"].iloc[0]
                    )
                    if model_total != 0:
                        text = f"Passed molecules: {passed} ({(passed / model_total * 100):.1f}%)"
                    else:
                        text = f"Passed molecules: {passed} (0%)"
                    ax.annotate(
                        text,
                        (bar_center_x, bar_center_y),
                        ha="center",
                        va="center",
                        fontsize=12,
                        color="black",
                        fontweight="bold",
                        bbox=dict(
                            facecolor="white", alpha=0.7, edgecolor="none", pad=3
                        ),
                        zorder=1000,
                    )

        for i, model in enumerate(models):
            count = data.loc[model, "num_mol"]
            max_bar_width = data["num_mol"].max()
            text_x_position = max_bar_width
            ax.text(
                text_x_position,
                i,
                int(count),
                va="center",
                ha="left",
                fontsize=12,
                color="black",
                fontweight="bold",
            )

        banned_bars = []
        banned_percentages = []
        ratio_cols = [
            col
            for col in data.columns
            if "banned_ratio" in col and "num_banned" not in col
        ]
        colors = get_model_colors(model_names=ratio_cols, cmap="Paired")
        for col, color in zip(ratio_cols, colors.values(), strict=False):
            num_banned_col = f"num_banned_{col}"
            ratio_name = col.replace("banned_ratio", "").strip("_")
            ratio_name = clean_name(ratio_name)

            banned_count = data[num_banned_col]
            total_banned = banned_count.sum()
            banned_percent = (total_banned / total_mols) * 100 if total_mols > 0 else 0
            banned_percentages.append(banned_percent)

            if banned_percent == 0.0:
                label = f"{ratio_name} (0%)"
            else:
                label = f"{ratio_name} ({format_number(total_banned)}, {banned_percent:.1f}%)"

            bar = ax.barh(x, banned_count, width, label=label, color=color, alpha=0.8)
            banned_bars.append(bar)

        clean_filter_name = filter_name.split("/")[-1]
        ax.set_title(
            clean_name(clean_filter_name), fontsize=14, pad=20, fontweight="bold"
        )
        ax.set_yticks(x)
        ax.set_yticklabels(models, fontsize=12)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(format_number))
        ax.set_xlim(left=0)

    ax.set_xlabel(
        f"Number of Molecules (Total: {format_number(total_mols)})",
        fontsize=12,
        labelpad=10,
    )
    ax.set_ylabel("Models", fontsize=12, labelpad=10)

    ax.grid(True, axis="x", alpha=0.2, linestyle="--")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    sorted_indices = np.argsort(banned_percentages)[::-1]
    sorted_handles = [total] + [banned_bars[i] for i in sorted_indices]

    legend = ax.legend(
        handles=sorted_handles,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=11,
        ncol=1,
    )
    legend.get_frame().set_alpha(0.9)
    legend.get_frame().set_edgecolor("lightgray")

    plt.subplots_adjust(right=0.85, hspace=0.6, wspace=0.5)

    plots_dir = struct_folder / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(
        plots_dir / "molecule_counts_comparison.png",
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
        edgecolor="none",
    )
    plt.close()


def plot_restriction_ratios(config, stage_dir):
    """Plot restriction ratios for structural filters."""
    folder_to_save = Path(process_path(config["folder_to_save"]))
    folder_name = config["folder_to_save"].split("/")[-1]

    config_structFilters = load_config(config["config_structFilters"])

    struct_folder = folder_to_save / stage_dir
    paths = list(struct_folder.glob("*/metrics.csv"))
    if not paths:
        paths = list(struct_folder.glob("*metrics.csv"))
    paths = [str(p) for p in paths]
    check_paths(config_structFilters, paths)

    if not paths:
        logger.error("No data files found in %s", str(folder_to_save))
        return

    filter_data = {}
    model_names_filters = {}
    for path in paths:
        filter_name = path.split(f"{folder_name}/")[-1].split("_metrics.csv")[0]
        data = pd.read_csv(path)

        ratio_cols = [col for col in data.columns if "banned_ratio" in col]
        model_names_filters[filter_name] = dict(
            zip(data["model_name"].tolist(), data["num_mol"].tolist(), strict=False)
        )

        if not ratio_cols:
            continue

        clean_cols = {
            col: col.replace("_banned_ratio", "")
            .replace("banned_ratio", "")
            .replace("_s", "s")
            for col in ratio_cols
        }
        ratios = data[ratio_cols].rename(columns=clean_cols)
        actual_model_names = data["model_name"].tolist()
        ratios.index = actual_model_names

        row = ratios.iloc[0]
        if row.isna().all():
            continue

        all_value = None
        if "all" in row.index:
            all_value = row["all"]
            row = row.drop("all")

        sorted_values = row.sort_values(ascending=False)

        if all_value is not None:
            if all_value >= sorted_values.iloc[0]:
                sorted_index = pd.Index(["all"]).append(sorted_values.index)
            else:
                sorted_index = sorted_values.index.append(pd.Index(["all"]))
            ratios = ratios[sorted_index]
        else:
            ratios = ratios[sorted_values.index]

        filter_data[filter_name] = ratios
    if not filter_data:
        logger.error("No valid data to plot")
        return

    model_names_filters = pd.DataFrame(model_names_filters).reset_index()
    model_names_filters = model_names_filters.rename(columns={"index": "model_name"})

    plt.style.use("default")
    sns.set_style("white")
    sns.set_context("talk")

    n_filters = len(filter_data)
    n_cols = min(2, n_filters)
    n_rows = (n_filters + n_cols - 1) // n_cols

    fig = plt.figure(figsize=(16, 7 * n_rows))
    fig.suptitle(
        "Comparison of Restriction Ratios Across Different Filters",
        fontsize=16,
        y=0.98,
        fontweight="bold",
    )

    for idx, (filter_name, data) in enumerate(filter_data.items()):
        number_of_mols = np.array(model_names_filters[filter_name].tolist())
        ax = plt.subplot(n_rows, n_cols, idx + 1)
        for col in data.columns:
            if col not in ["any", "all", "model_name", "num_mol"]:
                data[col] = number_of_mols * (1 - np.array(data[col].tolist()))
        if not data.empty and data.notna().any().any():
            if "all" in data.columns:
                data.drop(columns=["all"], inplace=True)
            if "any" in data.columns:
                data.drop(columns=["any"], inplace=True)

            custom_cmap = LinearSegmentedColormap.from_list(
                "custom", ["white", "#B29EEE"]
            )
            show_cbar = idx == 1
            sns.heatmap(
                data.T,
                cmap=custom_cmap,
                ax=ax,
                cbar_kws={"label": "Passed Molecules", "format": "%d"},
                vmin=0,
                vmax=max(data.max()),
                fmt=".0f",
                annot=True,
                annot_kws={"size": 12, "rotation": 0, "color": "black"},
                cbar=show_cbar,
            )

            ax.set_title(
                f"{clean_name(filter_name)} Filter", fontsize=12, fontweight="bold"
            )
            plt.setp(ax.get_yticklabels(), rotation=0, ha="right", fontsize=12)
            plt.setp(ax.get_xticklabels(), rotation=0, ha="right", fontsize=12)
            ax.set_xlabel("Model")

            actual_model_names = data.index.tolist()
            if len(actual_model_names) == len(ax.get_xticklabels()):
                ax.set_xticklabels(actual_model_names)

        else:
            ax.text(
                0.5,
                0.5,
                "No data available",
                horizontalalignment="center",
                verticalalignment="center",
                transform=ax.transAxes,
            )
            ax.set_title(
                f"{clean_name(filter_name)} Filter",
                pad=10,
                fontsize=11,
                fontweight="bold",
            )

    plt.tight_layout()
    plots_dir = struct_folder / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(
        plots_dir / "restriction_ratios_comparison.png",
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
        edgecolor="none",
    )
    plt.close()


def combine_filter_results_in_memory(output_dir, input_df, pass_mask_by_filter):
    """Combine per-filter pass masks in-memory and persist stage outputs."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    identity_cols = [c for c in ["smiles", "model_name", "mol_idx"] if c in input_df.columns]
    combined = input_df[identity_cols].copy()

    pass_columns = []
    for filter_name, pass_df in pass_mask_by_filter.items():
        col_name = f"pass_{filter_name}"
        pass_columns.append(col_name)

        if pass_df is None or len(pass_df) == 0:
            combined[col_name] = False
            continue

        local_df = pass_df.copy()
        pass_col = None
        if "pass" in local_df.columns:
            pass_col = "pass"
        elif "pass_filter" in local_df.columns:
            pass_col = "pass_filter"

        if pass_col is None:
            combined[col_name] = False
            continue

        merge_cols = [c for c in identity_cols if c in local_df.columns]
        if not merge_cols:
            combined[col_name] = False
            continue

        cols_to_take = merge_cols + [pass_col]
        local_df = local_df[cols_to_take].drop_duplicates(subset=merge_cols, keep="last")
        local_df = local_df.rename(columns={pass_col: col_name})
        combined = combined.merge(local_df, on=merge_cols, how="left")
        combined[col_name] = combined[col_name].fillna(False).astype(bool)

    if pass_columns:
        combined["_pass_all"] = combined[pass_columns].all(axis=1)
    else:
        combined["_pass_all"] = False

    filtered_df = combined[combined["_pass_all"]][identity_cols].copy()
    filtered_df.to_csv(output_path / "filtered_molecules.csv", index=False)

    failed_cols = identity_cols + pass_columns
    failed_df = combined[~combined["_pass_all"]][failed_cols].copy()
    failed_df.to_csv(output_path / "failed_molecules.csv", index=False)
    return filtered_df, failed_df


def filter_data(config, stage_dir):
    """Filter and combine data from all structural filters.

    Args:
        config: Configuration dictionary
        stage_dir: Stage directory path
    """
    base_folder = Path(process_path(config["folder_to_save"]))
    folder_to_save = base_folder / stage_dir

    paths = list(folder_to_save.glob("*/filtered_molecules.csv"))
    if not paths:
        paths = list(folder_to_save.glob("*filteredMols.csv"))
    paths = [str(p) for p in paths]

    columns_to_drop = ["pass", "any_pass", "name", "pass_any"]
    datas = []
    for path in paths:
        data = pd.read_csv(path)
        for col in columns_to_drop:
            if col in data.columns:
                data.drop(columns=[col], inplace=True)
        datas.append(data)

    if len(datas) > 0:
        filtered_data = datas[0].copy()

        for df in datas[1:]:
            merge_cols = ["smiles", "model_name"]
            existing_cols = set(filtered_data.columns) - set(merge_cols)
            new_cols = [
                col
                for col in df.columns
                if col not in existing_cols and col not in merge_cols
            ]

            if new_cols:
                cols_to_merge = merge_cols + new_cols
                filtered_data = filtered_data.merge(
                    df[cols_to_merge], on=merge_cols, how="inner"
                )
            else:
                filtered_data = filtered_data.merge(
                    df[merge_cols], on=merge_cols, how="inner"
                )
    else:
        filtered_data = pd.DataFrame(columns=["smiles", "model_name", "mol_idx"])

    if "mol_idx" not in filtered_data.columns:
        filtered_data["mol_idx"] = None

    cols = ["smiles", "model_name", "mol_idx"]
    out_df = filtered_data[cols].copy()
    out_df.to_csv(folder_to_save / "filtered_molecules.csv", index=False)

    is_post_descriptors = (
        "03_structural_filters_post" in stage_dir or stage_dir == "StructFilters"
    )
    if is_post_descriptors:
        descriptors_path = base_folder / "Descriptors" / "passDescriptorsSMILES.csv"
        if descriptors_path.exists():
            input_path = str(descriptors_path)
        else:
            sampled_path = base_folder / "sampled_molecules.csv"
            if sampled_path.exists():
                input_path = str(sampled_path)
            else:
                try:
                    from hedgehog.stages.structFilters.main import _get_input_path

                    input_path = _get_input_path(config, stage_dir, str(base_folder))
                except Exception:
                    input_path = None
    else:
        sampled_path = base_folder / "sampled_molecules.csv"
        if sampled_path.exists():
            input_path = str(sampled_path)
        else:
            try:
                from hedgehog.stages.structFilters.main import _get_input_path

                input_path = _get_input_path(config, stage_dir, str(base_folder))
            except Exception:
                input_path = None

    if input_path and Path(input_path).exists():
        try:
            all_input = pd.read_csv(input_path)
            if len(out_df) > 0:
                merge_cols = ["smiles", "model_name"]
                merged = all_input.merge(
                    out_df[merge_cols], on=merge_cols, how="left", indicator=True
                )
                fail_molecules = merged[merged["_merge"] == "left_only"].drop(
                    columns=["_merge"]
                )
            else:
                fail_molecules = all_input.copy()

            if len(fail_molecules) > 0:
                extended_paths = [str(p) for p in folder_to_save.glob("*extended.csv")]
                all_extended = None
                for ext_path in extended_paths:
                    try:
                        ext_df = pd.read_csv(ext_path)
                        if all_extended is None:
                            all_extended = ext_df.copy()
                        else:
                            merge_cols = ["smiles", "model_name"]
                            pass_cols = [
                                col
                                for col in ext_df.columns
                                if col.startswith("pass_") or col == "pass"
                            ]
                            if pass_cols:
                                cols_to_merge = merge_cols + pass_cols
                                all_extended = all_extended.merge(
                                    ext_df[cols_to_merge],
                                    on=merge_cols,
                                    how="outer",
                                    suffixes=("", "_dup"),
                                )
                                for col in pass_cols:
                                    if f"{col}_dup" in all_extended.columns:
                                        all_extended[col] = all_extended[col].fillna(
                                            all_extended[f"{col}_dup"]
                                        )
                                        all_extended = all_extended.drop(
                                            columns=[f"{col}_dup"]
                                        )
                    except Exception:
                        continue

                if all_extended is not None:
                    merge_cols = ["smiles", "model_name"]
                    pass_cols = [
                        col
                        for col in all_extended.columns
                        if col.startswith("pass_") or col == "pass"
                    ]
                    if pass_cols:
                        cols_to_merge = merge_cols + pass_cols
                        fail_molecules = fail_molecules.merge(
                            all_extended[cols_to_merge], on=merge_cols, how="left"
                        )
                        for col in pass_cols:
                            if col in fail_molecules.columns:
                                fail_molecules[col] = fail_molecules[col].fillna(False)

                id_cols = ["smiles", "model_name", "mol_idx"]
                pass_cols_final = [
                    col
                    for col in fail_molecules.columns
                    if col.startswith("pass_") or col == "pass"
                ]
                fail_cols = [
                    c for c in id_cols if c in fail_molecules.columns
                ] + pass_cols_final
                fail_molecules[fail_cols].to_csv(
                    folder_to_save / "failed_molecules.csv", index=False
                )
        except Exception as e:
            logger.warning("Could not create failStructFiltersSMILES.csv: %s", e)

    return filtered_data


def inject_identity_columns_to_all_csvs(config, stage_dir):
    """Ensure identity columns are ordered consistently in all CSVs."""
    base_folder = Path(process_path(config["folder_to_save"]))
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


def _get_breakdown_folder(file_path):
    """Get the CommonAlertsBreakdown folder path, creating it if necessary."""
    path_to_save = Path(file_path).parent / "CommonAlertsBreakdown"
    path_to_save.mkdir(parents=True, exist_ok=True)
    return str(path_to_save)


def _save_plot(file_path, filename, dpi=600):
    """Save plot to CommonAlertsBreakdown folder and close it."""
    output_path = Path(_get_breakdown_folder(file_path)) / filename
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close()


def _aggregate_reasons(filter_reasons):
    """Aggregate all reasons across filters into a single dictionary."""
    all_reasons = {}
    for reasons in filter_reasons.values():
        for reason, count in reasons:
            all_reasons[reason] = all_reasons.get(reason, 0) + count
    return all_reasons


def _truncate_string(s, max_length, suffix="..."):
    """Truncate string if longer than max_length."""
    if len(s) > max_length:
        return s[: max_length - len(suffix)] + suffix
    return s


def _calculate_grid_layout(num_items):
    """Calculate rows and columns for grid layout."""
    if num_items <= 3:
        return 1, num_items
    if num_items <= 6:
        return 2, 3
    if num_items <= 9:
        return 3, 3
    if num_items <= 12:
        return 3, 4
    cols = 4
    rows = (num_items + cols - 1) // cols
    return rows, cols


def analyze_filter_failures(file_path):
    """Analyze filter failures from extended CSV file and generate visualizations."""
    logger.debug("Analyzing filter failures from: %s", file_path)
    df = pd.read_csv(file_path, low_memory=False)

    filter_columns = [
        col
        for col in df.columns
        if col.startswith("pass_") and col not in ("pass", "pass_any")
    ]

    if not filter_columns:
        return None, None, None

    filter_failures = {}
    filter_reasons = {}
    all_detailed_reasons = {}

    for col in filter_columns:
        filter_name = col.replace("pass_", "")
        failures = (~df[col]).sum()
        total = len(df)

        filter_failures[filter_name] = {
            "failures": failures,
            "total": total,
            "percentage": (failures / total) * 100,
        }

        reasons_col = f"reasons_{filter_name}"
        if reasons_col in df.columns:
            failed_molecules = df[~df[col]]
            reasons_data = failed_molecules[reasons_col].dropna()

            reason_counts = {}
            for reasons_str in reasons_data:
                if pd.notna(reasons_str) and str(reasons_str).strip():
                    for reason in str(reasons_str).split(";"):
                        reason = reason.strip()
                        if reason:
                            reason_counts[reason] = reason_counts.get(reason, 0) + 1

            filter_reasons[filter_name] = sorted(
                reason_counts.items(), key=lambda x: x[1], reverse=True
            )
            all_detailed_reasons[filter_name] = reason_counts

    _create_main_filter_plot(filter_failures, file_path)
    _create_individual_filter_plots(filter_failures, filter_reasons, file_path)
    _create_multi_panel_filter_plot(filter_failures, filter_reasons, file_path)

    all_reasons = _aggregate_reasons(filter_reasons)
    top_reasons = sorted(all_reasons.items(), key=lambda x: x[1], reverse=True)[:5]

    if top_reasons:
        logger.info(
            "Top 5 most common filter failure reasons (molecules may have multiple reasons):"
        )
        for i, (reason, count) in enumerate(top_reasons, 1):
            logger.info("  %d. %s: %d failures", i, reason, count)

    _create_complete_reasons_breakdown(all_detailed_reasons, filter_failures, file_path)
    _create_comprehensive_overview(filter_reasons, file_path)
    _create_summary_table(filter_failures, filter_reasons, file_path)

    return filter_failures, filter_reasons, all_detailed_reasons


def _create_main_filter_plot(filter_failures, file_path):
    """Create main filter failures bar chart."""
    plot_data = [
        {
            "filter": name,
            "failures": stats["failures"],
            "percentage": stats["percentage"],
        }
        for name, stats in filter_failures.items()
    ]
    plot_df = pd.DataFrame(plot_data).sort_values("failures", ascending=False)

    plt.figure(figsize=(max(16, len(plot_df) * 0.6), 16))
    plt.bar(
        range(len(plot_df)),
        plot_df["failures"],
        color="steelblue",
        alpha=0.8,
        width=0.3,
    )

    plt.xlabel("Filters", fontsize=20)
    plt.ylabel("Number of Molecules Failed", fontsize=20)
    plt.title(
        "Number of Molecules Failed by Each Filter", fontsize=26, fontweight="bold"
    )
    plt.xticks(
        range(len(plot_df)), plot_df["filter"], rotation=45, ha="right", fontsize=16
    )

    max_failures = max(plot_df["failures"])
    for i, (_, row) in enumerate(plot_df.iterrows()):
        plt.text(
            i,
            row["failures"] + max_failures * 0.01,
            f"{row['failures']}\n({row['percentage']:.1f}%)",
            ha="center",
            va="bottom",
            fontsize=14,
        )

    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    _save_plot(file_path, "filter_failures_plot.png")


def _create_individual_filter_plots(filter_failures, filter_reasons, file_path):
    """Create individual plots for each filter showing failure reasons."""
    for filter_name, stats in filter_failures.items():
        if stats["failures"] == 0:
            continue

        reasons_data = filter_reasons.get(filter_name, [])
        if not reasons_data:
            continue

        plot_data = [
            {
                "Reason": reason,
                "Count": count,
                "Percentage": (count / stats["failures"]) * 100,
            }
            for reason, count in reasons_data
        ]
        plot_df = pd.DataFrame(plot_data).sort_values("Count", ascending=False)

        plt.figure(figsize=(max(16, len(plot_df) * 0.6), 20))
        plt.bar(
            range(len(plot_df)),
            plot_df["Count"],
            color="steelblue",
            alpha=0.8,
            width=0.3,
        )

        plt.xlabel("Failure Reasons", fontsize=20)
        plt.ylabel("Number of Molecules Failed", fontsize=20)
        plt.title(
            f"{filter_name.upper()} - Failure Reasons ({len(plot_df)} reasons, {stats['failures']} total failures)",
            fontsize=26,
            fontweight="bold",
        )
        plt.xticks(
            range(len(plot_df)),
            plot_df["Reason"],
            rotation=45,
            ha="right",
            fontsize=max(10, min(16, 300 // len(plot_df))),
        )

        max_count = max(plot_df["Count"])
        for i, (_, row) in enumerate(plot_df.iterrows()):
            plt.text(
                i,
                row["Count"] + max_count * 0.01,
                f"{row['Count']}\n({row['Percentage']:.1f}%)",
                ha="center",
                va="bottom",
                fontsize=12,
            )

        plt.grid(axis="y", alpha=0.3)
        plt.tight_layout()
        _save_plot(file_path, f"{filter_name}_reasons_plot.png")


def _create_multi_panel_filter_plot(filter_failures, filter_reasons, file_path):
    """Create multi-panel plot showing all filters with reasons."""
    sorted_filters = [
        (name, stats)
        for name, stats in sorted(
            filter_failures.items(), key=lambda x: x[1]["failures"], reverse=True
        )
        if stats["failures"] > 0
    ]

    if not sorted_filters:
        return

    rows, cols = _calculate_grid_layout(len(sorted_filters))
    plt.figure(figsize=(cols * 6, rows * 6))

    for i, (filter_name, _stats) in enumerate(sorted_filters):
        all_reasons_data = filter_reasons.get(filter_name, [])
        reason_names = [r[0] for r in all_reasons_data]
        reason_counts = [r[1] for r in all_reasons_data]

        if len(reason_names) > 10:
            title_suffix = f"(Top 10 of {len(all_reasons_data)} reasons)"
            reason_names = reason_names[:10]
            reason_counts = reason_counts[:10]
        else:
            title_suffix = f"({len(all_reasons_data)} reasons)"

        plt.subplot(rows, cols, i + 1)
        plt.bar(
            range(len(reason_names)),
            reason_counts,
            color="steelblue",
            alpha=0.8,
            width=0.3,
        )

        plt.xlabel("Reasons", fontsize=14)
        plt.ylabel("Molecules Failed", fontsize=14)
        plt.title(
            f"{filter_name.upper()}\n{title_suffix}", fontsize=16, fontweight="bold"
        )

        truncated_names = [_truncate_string(name, 15, "...") for name in reason_names]
        plt.xticks(
            range(len(truncated_names)),
            truncated_names,
            rotation=45,
            ha="right",
            fontsize=12,
        )

        max_count = max(reason_counts) if reason_counts else 0
        for j, count in enumerate(reason_counts):
            plt.text(
                j,
                count + max_count * 0.01,
                f"{count}",
                ha="center",
                va="bottom",
                fontsize=11,
            )

        plt.grid(axis="y", alpha=0.3)

    plt.tight_layout(h_pad=1.5, w_pad=1.0)
    _save_plot(file_path, "all_filters_reasons_plot.png")


def _create_complete_reasons_breakdown(
    all_detailed_reasons, filter_failures, file_path
):
    """Create complete CSV breakdown of all reasons."""
    breakdown_data = []
    for filter_name, reasons_dict in all_detailed_reasons.items():
        total_failures = filter_failures[filter_name]["failures"]
        for reason, count in sorted(
            reasons_dict.items(), key=lambda x: x[1], reverse=True
        ):
            breakdown_data.append(
                {
                    "Ruleset": filter_name,
                    "Reason": reason,
                    "Count": count,
                    "Percentage_of_Filter_Failures": (count / total_failures) * 100
                    if total_failures > 0
                    else 0,
                    "Total_Filter_Failures": total_failures,
                }
            )

    breakdown_df = pd.DataFrame(breakdown_data)
    output_path = (
        Path(_get_breakdown_folder(file_path)) / "complete_reasons_breakdown.csv"
    )
    breakdown_df.to_csv(output_path, index=False)
    return breakdown_df


def _create_comprehensive_overview(filter_reasons, file_path):
    """Create comprehensive overview plot of most common failure reasons."""
    all_reasons = _aggregate_reasons(filter_reasons)
    top_reasons = sorted(all_reasons.items(), key=lambda x: x[1], reverse=True)

    if not top_reasons:
        return

    display_count = min(30, len(top_reasons))
    reason_names = [
        _truncate_string(r[0], 30, "...") for r in top_reasons[:display_count]
    ]
    reason_counts = [r[1] for r in top_reasons[:display_count]]

    plt.figure(figsize=(max(16, display_count * 0.6), 16))
    plt.bar(
        range(len(reason_names)), reason_counts, color="darkgreen", alpha=0.7, width=0.3
    )

    plt.xlabel("Failure Reasons", fontsize=20)
    plt.ylabel("Total Number of Molecules Failed", fontsize=20)
    plt.title(
        f"Most Common Molecular Filter Failure Reasons (Top {display_count} of {len(top_reasons)})",
        fontsize=26,
        fontweight="bold",
    )
    plt.xticks(
        range(len(reason_names)), reason_names, rotation=45, ha="right", fontsize=16
    )

    max_count = max(reason_counts)
    for i, count in enumerate(reason_counts):
        plt.text(
            i,
            count + max_count * 0.01,
            f"{count}",
            ha="center",
            va="bottom",
            fontsize=14,
        )

    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    _save_plot(file_path, "comprehensive_reasons_overview.png")

    # Save CSV summary
    all_reasons_df = pd.DataFrame(top_reasons, columns=["Reason", "Total_Count"])
    output_path = Path(_get_breakdown_folder(file_path)) / "all_reasons_summary.csv"
    all_reasons_df.to_csv(output_path, index=False)


def _create_summary_table(filter_failures, filter_reasons, file_path):
    """Create summary table CSV with filter statistics."""
    summary_data = []
    for filter_name, stats in filter_failures.items():
        row = {
            "Ruleset": filter_name,
            "Total_Failures": stats["failures"],
            "Failure_Percentage": stats["percentage"],
            "Total_Molecules": stats["total"],
            "Unique_Reasons_Count": len(filter_reasons.get(filter_name, [])),
        }

        reasons = filter_reasons.get(filter_name, [])
        for i, (reason, count) in enumerate(reasons[:5], 1):
            row[f"Top_Reason_{i}"] = reason
            row[f"Top_Reason_{i}_Count"] = count
            row[f"Top_Reason_{i}_Percentage"] = (
                (count / stats["failures"]) * 100 if stats["failures"] > 0 else 0
            )

        summary_data.append(row)

    summary_df = pd.DataFrame(summary_data).sort_values(
        "Total_Failures", ascending=False
    )
    output_path = Path(_get_breakdown_folder(file_path)) / "filter_summary_table.csv"
    summary_df.to_csv(output_path, index=False)
    logger.debug("Summary table saved to: %s", output_path)
    return summary_df


def plot_filter_failures_analysis(config, stage_dir):
    """Analyze and plot filter failures for extended CSV files."""
    folder_to_save = Path(process_path(config["folder_to_save"]))
    struct_folder = folder_to_save / stage_dir

    if not struct_folder.exists():
        return

    extended_files = [str(p) for p in struct_folder.glob("*_extended.csv")]

    if not extended_files:
        logger.debug("No extended CSV files found for failure analysis")
        return

    for file_path in extended_files:
        try:
            analyze_filter_failures(file_path)
        except Exception as e:
            logger.debug("Error analyzing filter failures for %s: %s", file_path, e)
