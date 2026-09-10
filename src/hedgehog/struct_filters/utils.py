import contextlib
import importlib
import io
import json
import os
import sys
import threading
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
_ALERT_COMPILED_SMARTS: list[dict] | None = None
_ALERT_RULESET_NAMES: list[str] | None = None
_MOLGRAPH_CATALOG = None
_MOLGRAPH_SEVERITY_BY_ENTRY: list[int] | None = None
_MOLGRAPH_RULE_DETAILS_BY_ENTRY: list[dict] | None = None
_PROTECTING_GROUPS_CATALOG = None
_SPECIAL_PROTECTING_GROUP_QUERIES = None


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


def _silence_worker_stdio() -> None:
    """Redirect worker stdout/stderr to /dev/null."""
    try:
        devnull_fd = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull_fd, 1)
        os.dup2(devnull_fd, 2)
        os.close(devnull_fd)
    except OSError:
        pass


def _resolve_scheduler(config_struct_filters, scheduler_key, default="processes"):
    """Resolve scheduler name for medchem parallel APIs."""
    scheduler_raw = config_struct_filters.get(scheduler_key)
    if scheduler_raw is None:
        scheduler_raw = config_struct_filters.get("parallel_scheduler", default)
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
        complexity_filter = mc.complexity.ComplexityFilter(
            complexity_metric=metric_name
        )
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


def _init_molcomplexity_worker_quiet(alert_names: list[str]) -> None:
    _silence_worker_stdio()
    _init_molcomplexity_worker(alert_names)


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


mc = _import_medchem_quietly()

try:
    from medchem.structural.lilly_demerits import LillyDemeritsFilters

    LILLY_AVAILABLE = True
except ImportError:
    LILLY_AVAILABLE = False
    LillyDemeritsFilters = None

from hedgehog._constants import CFG_STRUCT_FILTERS, KEY_FOLDER_TO_SAVE
from hedgehog.configs.logger import load_config, logger
from hedgehog.struct_filters._helpers import resolve_include_rulesets
from hedgehog.struct_filters.common_alert_diagnostics import (
    HITS_JSON_COLUMN,
    hit_records_to_json,
    make_compiled_alert_rule,
    make_hit_record,
)
from hedgehog.utils.datamol_import import import_datamol_quietly
from hedgehog.utils.parallel import parallel_map, resolve_n_jobs
from hedgehog.utils.paths import process_path as _shared_process_path

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


def process_path(folder_to_save, key_word=None):
    """Backward-compatible wrapper around shared process_path helper."""
    return _shared_process_path(folder_to_save, key_word)


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
    """Select alerts by whole ruleset, optionally dropping exact SMARTS."""
    df = pd.read_csv(config["alerts_data_path"])
    available = df["rule_set_name"].dropna().astype(str).tolist()
    available_unique = list(dict.fromkeys(available))
    include_rulesets = set(
        resolve_include_rulesets(config.get("include_rulesets"), available_unique)
    )
    mask = df["rule_set_name"].astype(str).isin(include_rulesets)

    excluded_smarts = {
        str(value)
        for value in (config.get("exclude_smarts") or [])
        if value is not None
    }
    if excluded_smarts:
        mask &= ~df["smarts"].astype(str).isin(excluded_smarts)

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


def prepare_structfilters_input(df, subsample, n_jobs, progress_cb=None):
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
    parse_payload = [(row_idx, item[0], item[1], item[2]) for row_idx, item in items]

    worker_count = max(1, int(n_jobs))
    parsed = parallel_map(
        _parse_smiles_item,
        parse_payload,
        worker_count,
        progress=progress_cb,
        initializer=_silence_worker_stdio
        if worker_count > 1 and len(parse_payload) > 1
        else None,
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


def _filter_valid_molecule_tuples(smiles):
    """Filter out invalid molecules from a list of (smiles, model, mol, ...) tuples."""
    cleaned = []
    for item in smiles:
        if len(item) < 3:
            continue
        smi, model, mol = item[0], item[1], item[2]
        if mol is None:
            continue
        if len(item) >= 4:
            cleaned.append((smi, model, mol, item[3]))
        else:
            cleaned.append((smi, model, mol))
    return cleaned


def _invoke_filter(apply_filter, config, mols, smiles, progress_cb=None):
    """Invoke a filter function with optional progress callback.

    Tries passing progress_cb as a keyword argument first; falls back to
    calling without it if the filter does not accept that parameter.
    """
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

    smiles = _filter_valid_molecule_tuples(smiles)
    mols = [it[2] for it in smiles]

    if len(mols) == 0:
        return None

    return _invoke_filter(apply_filter, config, mols, smiles, progress_cb)


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

    is_tuple_format = isinstance(smiles[0], tuple)
    if is_tuple_format:
        smiles = _filter_valid_molecule_tuples(smiles)
        mols = [it[2] for it in smiles]
    else:
        mols, smiles = dropna(mols, smiles)

    assert len(mols) == len(smiles), f"{len(mols)}, {len(smiles)}"
    if not is_tuple_format:
        assert len(mols) <= subsample

    for mol, smi in zip(mols, smiles, strict=False):
        smi_val = smi[0] if isinstance(smi, tuple) else smi
        assert mol is not None, f"{smi_val}"

    if len(mols) == 0:
        return None

    if is_tuple_format:
        return _invoke_filter(apply_filter, config, mols, smiles, progress_cb)

    if progress_cb is not None:
        try:
            return apply_filter(config, mols, progress_cb=progress_cb)
        except TypeError:
            pass
    return apply_filter(config, mols)


def _align_result_length(final_result, expected_len, smiles_with_model):
    """Pad or trim final_result to match expected_len.

    Returns the adjusted DataFrame.
    """
    actual_len = len(final_result)
    if actual_len == expected_len:
        return final_result

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

    return final_result


def _realign_result_by_temp_mol_idx(final_result):
    """Restore original molecule order when worker code returns ``_mol_idx``."""
    if "_mol_idx" not in final_result.columns:
        return final_result
    return final_result.sort_values("_mol_idx").reset_index(drop=True)


def add_model_name_col(final_result, smiles_with_model):
    """Add smiles, model_name, and mol_idx columns from input data."""
    expected_len = len(smiles_with_model)
    final_result = _realign_result_by_temp_mol_idx(final_result)
    final_result = _align_result_length(final_result, expected_len, smiles_with_model)

    smiles_vals = [item[0] for item in smiles_with_model]
    model_vals = [
        item[1] if item[1] is not None else "single" for item in smiles_with_model
    ]
    mol_idx_vals = [item[3] if len(item) >= 4 else None for item in smiles_with_model]

    final_result["smiles"] = smiles_vals
    final_result["model_name"] = model_vals
    final_result["mol_idx"] = mol_idx_vals
    if "_mol_idx" in final_result.columns:
        final_result = final_result.drop(columns=["_mol_idx"])

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

        hit_records = []
        for rule in compiled_smarts:
            ruleset = rule["ruleset"]
            patt = rule["patt"]
            desc = rule["description"]
            matches = mol.GetSubstructMatches(patt, uniquify=True)
            if matches:
                row_data[f"pass_{ruleset}"] = False
                if desc and desc not in reasons_map[ruleset]:
                    reasons_map[ruleset].append(desc)
                for match_id, atom_ids in enumerate(matches):
                    hit_records.append(
                        make_hit_record(
                            mol=mol,
                            rule=rule,
                            atom_ids=atom_ids,
                            match_id=match_id,
                            n_matches_for_rule=len(matches),
                        )
                    )

        for name in rule_set_names:
            row_data[f"reasons_{name}"] = ";".join(reasons_map[name])
        row_data[HITS_JSON_COLUMN] = hit_records_to_json(hit_records)
    else:
        for name in rule_set_names:
            row_data[f"pass_{name}"] = False
            row_data[f"reasons_{name}"] = "invalid_molecule"
        row_data[HITS_JSON_COLUMN] = hit_records_to_json([])

    return row_data


def _init_alert_worker(
    compiled_smarts: list[dict],
    rule_set_names: list[str],
) -> None:
    global _ALERT_COMPILED_SMARTS, _ALERT_RULESET_NAMES
    _ALERT_COMPILED_SMARTS = compiled_smarts
    _ALERT_RULESET_NAMES = rule_set_names


def _init_alert_worker_quiet(
    compiled_smarts: list[dict],
    rule_set_names: list[str],
) -> None:
    _silence_worker_stdio()
    _init_alert_worker(compiled_smarts, rule_set_names)


def _compile_alert_smarts(alert_data) -> list[dict]:
    """Pre-compile SMARTS patterns from alert data to avoid per-molecule recompilation."""
    compiled: list[dict] = []
    for _, row in alert_data.iterrows():
        smarts = row.get("smarts")
        ruleset = row.get("rule_set_name")
        desc = row.get("description") or ""
        if not isinstance(smarts, str) or not isinstance(ruleset, str):
            continue
        patt = Chem.MolFromSmarts(smarts)
        if patt is not None:
            compiled.append(
                make_compiled_alert_rule(
                    row,
                    ruleset=ruleset,
                    patt=patt,
                    smarts=smarts,
                    description=str(desc),
                )
            )
    return compiled


def _setup_alerts_progress(total_items, progress_cb):
    """Set up progress tracking and heartbeat thread for structural alerts.

    Returns (progress_wrapper, heartbeat_stop, heartbeat_thread).
    """
    progress_log_step = max(1, min(500, total_items // 20)) if total_items else 1
    heartbeat_interval_seconds = 30.0
    logger.info(
        "Common Alerts progress logging: step=%d molecules, heartbeat=%.0fs",
        progress_log_step,
        heartbeat_interval_seconds,
    )

    done_count = 0
    next_progress_log = progress_log_step
    progress_lock = threading.Lock()
    heartbeat_stop = threading.Event()

    def _progress_wrapper(done: int, total: int) -> None:
        nonlocal done_count, next_progress_log
        with progress_lock:
            done_count = done

        should_log = done == total or done >= next_progress_log
        if should_log:
            logger.info(
                "Common Alerts progress: %d/%d molecules (%.1f%%)",
                done,
                total,
                100.0 * done / max(1, total),
            )
            while next_progress_log <= done:
                next_progress_log += progress_log_step

        if progress_cb is not None:
            progress_cb(done, total)

    def _heartbeat_logger() -> None:
        while not heartbeat_stop.wait(heartbeat_interval_seconds):
            with progress_lock:
                current_done = done_count
            if total_items == 0:
                break
            if current_done == 0:
                logger.info(
                    "Common Alerts still running: 0/%d molecules complete; waiting for first worker batch",
                    total_items,
                )
                continue
            logger.info(
                "Common Alerts still running: %d/%d molecules (%.1f%%)",
                current_done,
                total_items,
                100.0 * current_done / max(1, total_items),
            )

    logger.info(
        "Common Alerts progress: 0/%d molecules (0.0%%)",
        total_items,
    )

    heartbeat_thread = None
    if total_items > 1:
        heartbeat_thread = threading.Thread(
            target=_heartbeat_logger,
            name="common-alerts-heartbeat",
            daemon=True,
        )
        heartbeat_thread.start()

    return _progress_wrapper, heartbeat_stop, heartbeat_thread


def _build_alerts_results(results_list, mols):
    """Reconstruct alerts results DataFrame from parallel worker output."""
    # Restore mol column (RDKit Mol objects are not picklable across processes).
    for row_data in results_list:
        idx = row_data.get("_mol_idx")
        row_data["mol"] = mols[idx]

    results = pd.DataFrame(results_list)
    if "_mol_idx" in results.columns:
        # apply_structural_alerts runs workers unordered for throughput; sort back to
        # the original molecule order before identity columns are attached.
        results = results.sort_values("_mol_idx").reset_index(drop=True)

    # Calculate pass and pass_any columns.
    pass_cols = [
        c
        for c in results.columns
        if c.startswith("pass_") and c not in {"pass", "pass_any"}
    ]
    results["pass"] = results[pass_cols].all(axis=1)
    results["pass_any"] = results[pass_cols].any(axis=1)
    if "_mol_idx" in results.columns:
        results = results.drop(columns=["_mol_idx"])
    return results


def _resolve_common_alerts_start_method(config_sf: dict) -> str | None:
    """Resolve multiprocessing start method for Common Alerts workers.

    Priority:
    1. config_structFilters.yml: ``common_alerts_start_method``
    2. environment: ``HEDGEHOG_COMMON_ALERTS_START_METHOD``
    3. platform default: ``fork`` on Linux to reduce worker warm-up latency
    """
    raw = config_sf.get("common_alerts_start_method")
    if raw is None:
        raw = os.environ.get("HEDGEHOG_COMMON_ALERTS_START_METHOD")

    if raw is None:
        return "fork" if sys.platform.startswith("linux") else None

    method = str(raw).strip().lower()
    if method in {"", "auto", "default", "none"}:
        return None
    return method


def apply_structural_alerts(
    config, mols, smiles_model_name_mols=None, progress_cb=None
):
    logger.info("Calculating Common Alerts...")

    # Load config and filter alerts ONCE outside
    config_sf = load_config(config[CFG_STRUCT_FILTERS])
    alert_data = filter_alerts(config_sf)

    rule_set_names = alert_data["rule_set_name"].unique().tolist()

    logger.info(
        "Processing %d filtered molecules with %d alert rule sets",
        len(mols),
        len(rule_set_names),
    )

    compiled_smarts = _compile_alert_smarts(alert_data)

    items = [(i, mol) for i, mol in enumerate(mols)]
    total_items = len(items)
    n_jobs = resolve_n_jobs(config_sf, config)
    start_method = _resolve_common_alerts_start_method(config_sf)
    logger.info("Common Alerts workers: %d", n_jobs)
    logger.info(
        "Common Alerts start method: %s",
        start_method if start_method else "auto",
    )

    progress_wrapper, heartbeat_stop, heartbeat_thread = _setup_alerts_progress(
        total_items, progress_cb
    )

    worker_initializer = (
        _init_alert_worker_quiet
        if n_jobs > 1 and len(items) > 1
        else _init_alert_worker
    )
    try:
        results_list = parallel_map(
            _check_alerts_single_mol,
            items,
            n_jobs,
            progress=progress_wrapper,
            initializer=worker_initializer,
            initargs=(compiled_smarts, rule_set_names),
            preserve_order=False,
            start_method=start_method,
        )
    finally:
        heartbeat_stop.set()
        if heartbeat_thread is not None:
            heartbeat_thread.join(timeout=0.1)

    logger.info(
        "Common Alerts completed: %d/%d molecules (%.1f%%)",
        total_items,
        total_items,
        100.0 if total_items else 0.0,
    )

    results = _build_alerts_results(results_list, mols)

    if smiles_model_name_mols is not None:
        results = add_model_name_col(results, smiles_model_name_mols)
    return results


def _get_molgraph_catalog_and_severity():
    global _MOLGRAPH_CATALOG
    global _MOLGRAPH_SEVERITY_BY_ENTRY
    global _MOLGRAPH_RULE_DETAILS_BY_ENTRY

    if (
        _MOLGRAPH_CATALOG is None
        or _MOLGRAPH_SEVERITY_BY_ENTRY is None
        or _MOLGRAPH_RULE_DETAILS_BY_ENTRY is None
    ):
        graph_path = mc.utils.loader.get_data_path("graph.csv")
        graph_df = pd.read_csv(graph_path)
        _MOLGRAPH_SEVERITY_BY_ENTRY = graph_df["severity"].astype(int).tolist()
        _MOLGRAPH_RULE_DETAILS_BY_ENTRY = [
            {
                "pattern_index": int(index),
                "label": row.get("labels", ""),
                "description": row.get("description", ""),
                "smarts": row.get("smarts", ""),
                "severity": int(row["severity"]),
            }
            for index, row in graph_df.iterrows()
        ]
        _MOLGRAPH_CATALOG = mc.catalogs.NamedCatalogs.unstable_graph(
            severity_threshold=1
        )
    return _MOLGRAPH_CATALOG, _MOLGRAPH_SEVERITY_BY_ENTRY


def _molgraph_diagnostics_for_mol(mol):
    catalog, severity_by_entry = _get_molgraph_catalog_and_severity()
    matched_entries = catalog.GetMatches(mol)
    max_severity = 0
    matched_rules = []
    for entry in matched_entries:
        try:
            entry_idx = int(entry.GetDescription())
        except (TypeError, ValueError):
            continue
        if 0 <= entry_idx < len(severity_by_entry):
            severity = int(severity_by_entry[entry_idx])
            max_severity = max(max_severity, severity)
            if _MOLGRAPH_RULE_DETAILS_BY_ENTRY is not None:
                matched_rules.append(_MOLGRAPH_RULE_DETAILS_BY_ENTRY[entry_idx])
    return {
        "molgraph_max_severity": max_severity,
        "molgraph_matched_patterns_json": json.dumps(
            matched_rules, separators=(",", ":")
        ),
        "molgraph_matched_pattern_count": len(matched_rules),
    }


def _compute_molgraph_diagnostics(mols, n_jobs, scheduler):
    if len(mols) == 0:
        return []
    return dm.parallelized(
        _molgraph_diagnostics_for_mol,
        mols,
        n_jobs=n_jobs,
        scheduler=scheduler,
        progress=False,
    )


def apply_molgraph_stats(config, mols, smiles_model_name_mols=None, progress_cb=None):
    logger.info("Calculating Molecular Graph statistics...")
    total_mols = len(mols)
    config_sf = load_config(config[CFG_STRUCT_FILTERS])
    n_jobs = resolve_n_jobs(config_sf, config)
    scheduler = _resolve_scheduler(config_sf, "molgraph_scheduler")
    logger.info("MolGraph workers: %d", n_jobs)

    diagnostics = _compute_molgraph_diagnostics(mols, n_jobs, scheduler)
    max_severities = [row["molgraph_max_severity"] for row in diagnostics]
    results = {
        "mol": mols,
        "molgraph_max_severity": max_severities,
        "molgraph_matched_patterns_json": [
            row["molgraph_matched_patterns_json"] for row in diagnostics
        ],
        "molgraph_matched_pattern_count": [
            row["molgraph_matched_pattern_count"] for row in diagnostics
        ],
    }
    if progress_cb is not None:
        progress_cb(total_mols, max(1, total_mols))

    for s in range(1, 12):
        results[f"pass_{s}"] = [sev < s for sev in max_severities]
    logger.info("MolGraph derived pass_1..pass_11 from one severity catalog pass")
    results = pd.DataFrame(results)

    if smiles_model_name_mols is not None:
        results = add_model_name_col(results, smiles_model_name_mols)
    return results


def apply_molcomplexity_filters(config, mols, smiles_model_name_mols=None):
    logger.info("Calculating Complexity filters...")
    config_path = config.get(CFG_STRUCT_FILTERS)
    config_sf = load_config(config_path) if config_path else {}
    n_jobs = resolve_n_jobs(config_sf, config)
    logger.info("MolComplexity workers: %d", n_jobs)

    alert_names = mc.complexity.ComplexityFilter.list_default_available_filters()
    items = [(i, mol) for i, mol in enumerate(mols)]
    rows = parallel_map(
        _compute_molcomplexity_one,
        items,
        n_jobs,
        initializer=(
            _init_molcomplexity_worker_quiet
            if n_jobs > 1 and len(items) > 1
            else _init_molcomplexity_worker
        ),
        initargs=(alert_names,),
    )
    for row in rows:
        mol_idx = row.pop("_mol_idx")
        row["mol"] = mols[mol_idx]
    final_result = pd.DataFrame(rows)

    if smiles_model_name_mols is not None:
        final_result = add_model_name_col(final_result, smiles_model_name_mols)
    return final_result


def _apply_simple_medchem_filter(
    config,
    mols,
    smiles_model_name_mols,
    filter_name,
    mc_func,
    scheduler_key,
    extra_kwargs=None,
):
    """Generic helper for simple medchem filters that share the same skeleton.

    Args:
        config: Pipeline configuration dictionary.
        mols: List of RDKit molecule objects.
        smiles_model_name_mols: Optional list of (smiles, model_name, mol, mol_idx) tuples.
        filter_name: Human-readable filter name for logging.
        mc_func: The medchem functional filter callable.
        scheduler_key: Config key used to resolve the parallel scheduler.
        extra_kwargs: Optional callable ``(config_sf) -> dict`` that returns
            additional keyword arguments for *mc_func*.
    """
    logger.info("Calculating %s filter...", filter_name)
    config_sf = load_config(config[CFG_STRUCT_FILTERS])
    n_jobs = resolve_n_jobs(config_sf, config)
    scheduler = _resolve_scheduler(config_sf, scheduler_key)
    logger.info("%s workers: %d", filter_name, n_jobs)
    kwargs = dict(
        mols=mols,
        n_jobs=n_jobs,
        scheduler=scheduler,
        progress=False,
        return_idx=False,
    )
    if extra_kwargs:
        kwargs.update(extra_kwargs(config_sf))
    out = mc_func(**kwargs)
    results = pd.DataFrame({"mol": mols, "pass": out})
    results["status"] = np.where(results["pass"], "ok", "warning")
    results["reason"] = np.where(
        results["pass"],
        "",
        f"{filter_name} diagnostic criterion triggered",
    )
    if smiles_model_name_mols is not None:
        results = add_model_name_col(results, smiles_model_name_mols)
    return results


def apply_bredt_filter(config, mols, smiles_model_name_mols=None):
    return _apply_simple_medchem_filter(
        config,
        mols,
        smiles_model_name_mols,
        "Bredt",
        mc.functional.bredt_filter,
        "bredt_scheduler",
    )


def _get_protecting_groups_catalog():
    """Load Hedgehog's curated protecting-group catalog once per worker."""
    global _PROTECTING_GROUPS_CATALOG
    if _PROTECTING_GROUPS_CATALOG is None:
        catalog_path = Path(__file__).parent / "data" / "protecting_groups.csv"
        _PROTECTING_GROUPS_CATALOG = mc.groups.ChemicalGroup(
            "protecting_groups",
            groups_db=catalog_path,
        )
    return _PROTECTING_GROUPS_CATALOG


def _get_special_protecting_group_queries():
    """Compile motifs that cannot use MedChem's exact terminal matcher."""
    global _SPECIAL_PROTECTING_GROUP_QUERIES
    if _SPECIAL_PROTECTING_GROUP_QUERIES is None:
        catalog = _get_protecting_groups_catalog()
        row = catalog.data.loc[catalog.data["name"] == "n-tert-butoxymethyl"].iloc[0]
        query = Chem.MolFromSmarts(str(row["smarts"]))
        if query is None:
            raise ValueError("Invalid n-tert-butoxymethyl SMARTS")
        _SPECIAL_PROTECTING_GROUP_QUERIES = {"n-tert-butoxymethyl": query}
    return _SPECIAL_PROTECTING_GROUP_QUERIES


def _match_protecting_groups(mol):
    """Return curated terminal protecting-group names matched by one molecule."""
    if mol is None:
        return []
    matches = _get_protecting_groups_catalog().get_matches(
        mol,
        exact_match=True,
        terminal_only=True,
    )
    names = [] if matches is None else matches["name"].astype(str).tolist()
    # MedChem's source entry includes the protected imidazole itself, so its
    # exact+terminal combination can never match a substituted molecule. This
    # curated SMARTS describes the actual N-CH2-O-tBu protecting-group motif.
    for name, query in _get_special_protecting_group_queries().items():
        if name not in names and mol.HasSubstructMatch(query):
            names.append(name)
    return names


def apply_protecting_groups(config, mols, smiles_model_name_mols=None):
    """Apply the curated terminal protecting-group catalog."""
    logger.info("Calculating Protecting Groups filter...")
    config_sf = load_config(config[CFG_STRUCT_FILTERS])
    n_jobs = resolve_n_jobs(config_sf, config)
    scheduler = _resolve_scheduler(config_sf, "protecting_groups_scheduler")
    logger.info("Protecting Groups workers: %d", n_jobs)

    matched_groups = dm.parallelized(
        _match_protecting_groups,
        mols,
        n_jobs=n_jobs,
        scheduler=scheduler,
        progress=False,
    )
    matched_groups = [groups or [] for groups in matched_groups]
    passed = np.asarray([len(groups) == 0 for groups in matched_groups], dtype=bool)
    matched_text = [";".join(groups) for groups in matched_groups]

    results = pd.DataFrame(
        {
            "mol": mols,
            "pass": passed,
            "matched_protecting_groups": matched_text,
            "n_protecting_groups": [len(groups) for groups in matched_groups],
        }
    )
    results["status"] = np.where(results["pass"], "ok", "warning")
    results["reason"] = np.where(
        results["pass"],
        "",
        "Protecting groups detected: " + results["matched_protecting_groups"],
    )
    if smiles_model_name_mols is not None:
        results = add_model_name_col(results, smiles_model_name_mols)
    return results


def apply_ring_infraction(config, mols, smiles_model_name_mols=None):
    return _apply_simple_medchem_filter(
        config,
        mols,
        smiles_model_name_mols,
        "Ring Infraction",
        mc.functional.ring_infraction_filter,
        "ring_infraction_scheduler",
        extra_kwargs=lambda cfg: {
            "hetcycle_min_size": cfg.get("ring_infraction_hetcycle_min_size", 4),
        },
    )


def _compute_stereo_center_row(args):
    mol_idx, mol, max_stereo_centers, max_undefined_stereo_centers = args
    if mol is None:
        return {
            "_mol_idx": mol_idx,
            "pass": False,
            "n_stereo_centers": np.nan,
            "n_undefined_stereo_centers": np.nan,
            "stereo_max_centers": max_stereo_centers,
            "stereo_max_undefined": max_undefined_stereo_centers,
            "undefined_stereo_pass": False,
            "undefined_stereo_reason": "invalid molecule",
            "status": "warning",
            "reason": "invalid molecule",
        }

    prepared_mol = Chem.Mol(mol)
    Chem.AssignStereochemistry(prepared_mol, cleanIt=True, force=True)
    try:
        stereo_centers = Chem.FindMolChiralCenters(
            prepared_mol,
            includeUnassigned=True,
            useLegacyImplementation=False,
        )
    except RuntimeError:
        # RDKit's modern CIP labeler can reject otherwise parseable unusual
        # tetrahedral centres (for example [P@@H2]) with a carriers.size()
        # post-condition violation.  The legacy implementation still returns
        # the centre and its assignment, so one such molecule must not abort
        # the entire structural-filter stage.
        stereo_centers = Chem.FindMolChiralCenters(
            prepared_mol,
            includeUnassigned=True,
            useLegacyImplementation=True,
        )
    n_stereo_centers = len(stereo_centers)
    n_undefined_stereo_centers = sum(1 for _, label in stereo_centers if label == "?")
    passed = n_stereo_centers < max_stereo_centers
    undefined_stereo_pass = n_undefined_stereo_centers <= max_undefined_stereo_centers
    reasons = []
    if n_stereo_centers >= max_stereo_centers:
        reasons.append(
            f"stereocenters={n_stereo_centers} >= cutoff={max_stereo_centers}"
        )
    undefined_reason = ""
    if not undefined_stereo_pass:
        undefined_reason = (
            f"undefined_stereocenters={n_undefined_stereo_centers} "
            f"> maximum={max_undefined_stereo_centers}"
        )
    return {
        "_mol_idx": mol_idx,
        "pass": passed,
        "n_stereo_centers": n_stereo_centers,
        "n_undefined_stereo_centers": n_undefined_stereo_centers,
        "stereo_max_centers": max_stereo_centers,
        "stereo_max_undefined": max_undefined_stereo_centers,
        "undefined_stereo_pass": undefined_stereo_pass,
        "undefined_stereo_reason": undefined_reason,
        "status": "ok" if passed else "warning",
        "reason": "; ".join(reasons),
    }


def apply_stereo_center(config, mols, smiles_model_name_mols=None):
    logger.info("Calculating Stereo Center filter...")
    config_sf = load_config(config[CFG_STRUCT_FILTERS])
    n_jobs = resolve_n_jobs(config_sf, config)
    logger.info("Stereo Center workers: %d", n_jobs)

    max_stereo_centers = int(config_sf.get("stereo_max_centers", 4))
    max_undefined_stereo_centers = int(config_sf.get("stereo_max_undefined", 2))
    items = [
        (
            mol_idx,
            mol,
            max_stereo_centers,
            max_undefined_stereo_centers,
        )
        for mol_idx, mol in enumerate(mols)
    ]
    rows = parallel_map(_compute_stereo_center_row, items, n_jobs)
    result = pd.DataFrame(
        [
            {
                "mol": mols[row["_mol_idx"]],
                **{k: v for k, v in row.items() if k != "_mol_idx"},
            }
            for row in rows
        ]
    )
    if smiles_model_name_mols is not None:
        result = add_model_name_col(result, smiles_model_name_mols)
    return result


def _halogen_count(mol, atomic_number):
    if mol is None:
        return np.nan
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == atomic_number)


def apply_halogenicity(config, mols, smiles_model_name_mols=None):
    config_sf = load_config(config[CFG_STRUCT_FILTERS])
    thresholds = {
        "F": int(config_sf.get("halogenicity_thresh_F", 6)),
        "Br": int(config_sf.get("halogenicity_thresh_Br", 3)),
        "Cl": int(config_sf.get("halogenicity_thresh_Cl", 3)),
    }
    results = _apply_simple_medchem_filter(
        config,
        mols,
        None,
        "Halogenicity",
        mc.functional.halogenicity_filter,
        "halogenicity_scheduler",
        extra_kwargs=lambda _cfg: {
            "thresh_F": thresholds["F"],
            "thresh_Br": thresholds["Br"],
            "thresh_Cl": thresholds["Cl"],
        },
    )
    for symbol, atomic_number in (("F", 9), ("Br", 35), ("Cl", 17)):
        results[f"n_{symbol}"] = [_halogen_count(mol, atomic_number) for mol in mols]
        results[f"threshold_{symbol}"] = thresholds[symbol]
    results["reason"] = results.apply(
        lambda row: "; ".join(
            f"{symbol}={int(row[f'n_{symbol}'])} > cutoff={thresholds[symbol]}"
            for symbol in ("F", "Br", "Cl")
            if pd.notna(row[f"n_{symbol}"]) and row[f"n_{symbol}"] > thresholds[symbol]
        ),
        axis=1,
    )
    if smiles_model_name_mols is not None:
        results = add_model_name_col(results, smiles_model_name_mols)
    return results


def _compute_symmetry_row(args):
    mol_idx, mol, threshold = args
    if mol is None:
        return {
            "_mol_idx": mol_idx,
            "pass": False,
            "symmetry_score": np.nan,
            "symmetry_threshold": threshold,
            "status": "warning",
            "reason": "invalid molecule",
        }

    try:
        from medchem.utils.graph import score_symmetry

        score = float(score_symmetry(mol))
    except Exception as exc:
        return {
            "_mol_idx": mol_idx,
            # A calculation error is not evidence that the molecule violates
            # the symmetry threshold. Preserve it as an explicit diagnostic.
            "pass": True,
            "symmetry_score": np.nan,
            "symmetry_threshold": threshold,
            "status": "warning",
            "reason": f"symmetry calculation failed: {type(exc).__name__}: {exc}",
        }

    passed = score <= threshold
    reasons = []
    if not passed:
        reasons.append(f"symmetry={score:.6g} > cutoff={threshold}")
    return {
        "_mol_idx": mol_idx,
        "pass": passed,
        "symmetry_score": score,
        "symmetry_threshold": threshold,
        "status": "ok" if passed else "warning",
        "reason": "; ".join(reasons),
    }


def apply_symmetry(config, mols, smiles_model_name_mols=None):
    logger.info("Calculating Symmetry diagnostic...")
    config_sf = load_config(config[CFG_STRUCT_FILTERS])
    n_jobs = resolve_n_jobs(config_sf, config)
    threshold = float(config_sf.get("symmetry_threshold", 0.8))
    logger.info("Symmetry workers: %d; threshold: %g", n_jobs, threshold)
    rows = parallel_map(
        _compute_symmetry_row,
        [(mol_idx, mol, threshold) for mol_idx, mol in enumerate(mols)],
        n_jobs,
        chunksize=1,
    )
    results = pd.DataFrame(
        [
            {
                "mol": mols[row["_mol_idx"]],
                **{key: value for key, value in row.items() if key != "_mol_idx"},
            }
            for row in rows
        ]
    )
    if smiles_model_name_mols is not None:
        results = add_model_name_col(results, smiles_model_name_mols)
    return results


def apply_nibr_filter(config, mols, smiles_model_name_mols=None):
    logger.info("Calculating NIBR filter...")
    config_sf = load_config(config[CFG_STRUCT_FILTERS])
    n_jobs = resolve_n_jobs(config_sf, config)
    scheduler = _resolve_scheduler(config_sf, "nibr_scheduler")
    logger.info("NIBR workers: %d", n_jobs)

    indexed_mols = list(enumerate(mols))
    n_workers = max(1, min(n_jobs, len(indexed_mols)))
    chunked_mols = _split_indexed_mols(indexed_mols, n_workers)
    chunk_payloads = [(chunk, True) for chunk in chunked_mols]
    chunk_results = parallel_map(
        _process_nibr_chunk,
        chunk_payloads,
        n_workers,
        chunksize=1,
        initializer=_silence_worker_stdio if n_workers > 1 else None,
    )

    rows = []
    for chunk_rows in chunk_results:
        rows.extend(chunk_rows)
    results = pd.DataFrame(rows)

    expected_indices = set(range(len(mols)))
    observed_indices = (
        set(results["_mol_idx"].tolist()) if "_mol_idx" in results else set()
    )
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

    if smiles_model_name_mols is not None:
        results = add_model_name_col(results, smiles_model_name_mols)
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
        logger.warning(
            "Lilly batch failed (%s). Retrying %d molecules individually.",
            batch_error,
            len(batch),
        )
        return _process_lilly_one_by_one(dfilter, batch, scheduler)


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


def _resolve_lilly_filter_options(config_sf):
    """Build validated Lilly scorer options from the structural-filter config."""
    raw_cutoff = config_sf.get("lilly_demerit_cutoff")
    if raw_cutoff is None:
        return {}

    try:
        cutoff = int(raw_cutoff)
    except (TypeError, ValueError) as error:
        raise ValueError("lilly_demerit_cutoff must be a positive integer") from error
    if isinstance(raw_cutoff, float) and not raw_cutoff.is_integer():
        raise ValueError("lilly_demerit_cutoff must be a positive integer")
    if cutoff < 1:
        raise ValueError("lilly_demerit_cutoff must be a positive integer")
    return {"dthresh": cutoff}


def _new_lilly_filter(lilly_filter_options):
    """Create a Lilly scorer with the configured native demerit cutoff."""
    return LillyDemeritsFilters(**lilly_filter_options)


def _process_lilly_chunk(args):
    """Process one Lilly chunk in an isolated worker process."""
    indexed_chunk, scheduler, lilly_filter_options = args
    mol_indices = [item[0] for item in indexed_chunk]
    mol_chunk = [item[1] for item in indexed_chunk]

    dfilter = _new_lilly_filter(lilly_filter_options)
    try:
        chunk_result = dfilter(mols=mol_chunk, n_jobs=1, scheduler=scheduler)
    except Exception:
        chunk_result = _process_lilly_one_by_one(dfilter, mol_chunk, scheduler)

    if len(chunk_result) != len(indexed_chunk):
        template = chunk_result.iloc[-1].to_dict() if len(chunk_result) > 0 else None
        chunk_result = _ensure_dataframe_length(
            chunk_result, len(indexed_chunk), template
        )

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


def apply_lilly_filter(config, mols, smiles_model_name_mols=None):
    """Apply Lilly demerits filter to molecules."""
    if not LILLY_AVAILABLE:
        raise ImportError(
            "Lilly demerits filter is not available. "
            "This filter requires conda/mamba-installed binaries. "
            "Install with: conda install lilly-medchem-rules\n"
            "Or disable this filter by setting 'calculate_lilly: False' in config_sf.yml"
        )

    logger.info("Calculating Lilly filter...")
    config_sf = load_config(config[CFG_STRUCT_FILTERS])
    n_jobs = resolve_n_jobs(config_sf, config)
    scheduler = _resolve_scheduler(config_sf, "lilly_scheduler", default="threads")
    lilly_filter_options = _resolve_lilly_filter_options(config_sf)
    if scheduler != "threads":
        logger.warning(
            "Lilly supports only threads scheduler. Falling back to 'threads' (got '%s').",
            scheduler,
        )
        scheduler = "threads"
    logger.info("Lilly workers: %d", n_jobs)
    logger.info(
        "Lilly demerit cutoff: %s",
        lilly_filter_options.get("dthresh", "native default"),
    )

    if smiles_model_name_mols is not None:
        expected_len = len(smiles_model_name_mols)
        input_smiles = [
            item[0] if isinstance(item, tuple) and len(item) > 0 else None
            for item in smiles_model_name_mols
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
        if smiles_model_name_mols is not None:
            results = add_model_name_col(results, smiles_model_name_mols)
        return results

    # Run Lilly across worker processes, one thread-backed call per worker.
    n_workers = max(1, min(n_jobs, len(valid_indexed)))
    chunked_mols = _split_indexed_mols(valid_indexed, n_workers)
    chunk_payloads = [
        (chunk, scheduler, lilly_filter_options) for chunk in chunked_mols
    ]

    try:
        chunk_results = parallel_map(
            _process_lilly_chunk,
            chunk_payloads,
            n_workers,
            chunksize=1,
            initializer=_silence_worker_stdio if n_workers > 1 else None,
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
        dfilter = _new_lilly_filter(lilly_filter_options)
        results = _run_lilly_in_batches(dfilter, valid_mols, n_jobs, scheduler)
        if results is None:
            raise ValueError("All Lilly batches failed in fallback mode.") from error
        template = results.iloc[-1].to_dict() if len(results) > 0 else None
        results = _ensure_dataframe_length(results, len(valid_mols), template)
    else:
        expected_indices = set(valid_indices)
        observed_indices = (
            set(results["_mol_idx"].tolist())
            if "_mol_idx" in results.columns
            else set()
        )
        if len(results) != len(valid_mols) or observed_indices != expected_indices:
            logger.warning(
                "Lilly chunked processing mismatch (got %d/%d rows). Falling back to batched native mode.",
                len(results),
                len(valid_mols),
            )
            dfilter = _new_lilly_filter(lilly_filter_options)
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

    if smiles_model_name_mols is not None:
        results = add_model_name_col(results, smiles_model_name_mols)
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


def _stats_common_alerts(config, filter_results, model_name, num_mol, stat, extend):
    """Compute statistics for common_alerts filter."""
    res_df = _create_base_stats_df(
        model_name,
        num_mol,
        all_banned_ratio=(~filter_results["pass"]).mean(),
        any_banned_ratio=(~filter_results["pass_any"]).mean(),
    )
    ruleset_names = [
        column[len("pass_") :]
        for column in filter_results.columns
        if column.startswith("pass_") and column not in {"pass", "pass_any"}
    ]
    for name in ruleset_names:
        res_df[f"{name}_banned_ratio"] = 1 - filter_results[f"pass_{name}"].mean()
    return common_postprocessing_statistics(filter_results, res_df, stat, extend)


def _stats_molgraph(config, filter_results, model_name, num_mol, stat, extend):
    """Compute statistics for molgraph_stats filter."""
    res_df = _create_base_stats_df(model_name, num_mol)
    for i in range(1, 12):
        res_df[f"banned_ratio_s_{i}"] = 1 - filter_results[f"pass_{i}"].mean()
    res_df, filter_extended = common_postprocessing_statistics(
        filter_results, res_df, stat, extend
    )
    raw_threshold = config.get("molgraph_max_severity", 5)
    try:
        threshold = int(raw_threshold)
    except (TypeError, ValueError) as error:
        raise ValueError(
            "molgraph_max_severity must be an integer from 1 to 11"
        ) from error
    if isinstance(raw_threshold, float) and not raw_threshold.is_integer():
        raise ValueError("molgraph_max_severity must be an integer from 1 to 11")
    if not 1 <= threshold <= 11:
        raise ValueError("molgraph_max_severity must be an integer from 1 to 11")
    filter_extended["pass"] = filter_extended["molgraph_max_severity"] < threshold
    return res_df, filter_extended


def _stats_molcomplexity(config, filter_results, model_name, num_mol, stat, extend):
    """Compute statistics for molcomplexity filter."""
    res_df = _create_base_stats_df(
        model_name,
        num_mol,
        all_banned_ratio=1 - filter_results["pass"].mean(),
        any_banned_ratio=1 - filter_results["pass_any"].mean(),
    )
    for name in mc.complexity.ComplexityFilter.list_default_available_filters():
        res_df[f"{name}_banned_ratio"] = 1 - filter_results[f"pass_{name}"].mean()
    return common_postprocessing_statistics(filter_results, res_df, stat, extend)


def _stats_simple_banned(config, filter_results, model_name, num_mol, stat, extend):
    """Compute statistics for simple banned-ratio filters (bredt, protecting_groups, etc.)."""
    res_df = _create_base_stats_df(model_name, num_mol)
    res_df["banned_ratio"] = 1 - filter_results["pass"].mean()
    return common_postprocessing_statistics(filter_results, res_df, stat, extend)


def _stats_nibr(config, filter_results, model_name, num_mol, stat, extend):
    """Compute statistics using the published NIBR severity policy."""
    raw_threshold = config.get("nibr_max_severity", 10)
    try:
        threshold = int(raw_threshold)
    except (TypeError, ValueError) as error:
        raise ValueError("nibr_max_severity must be a positive integer") from error
    if isinstance(raw_threshold, float) and not raw_threshold.is_integer():
        raise ValueError("nibr_max_severity must be a positive integer")
    if threshold < 1:
        raise ValueError("nibr_max_severity must be a positive integer")

    severity = pd.to_numeric(filter_results["severity"], errors="coerce")
    severity_pass = severity.lt(threshold).fillna(False)
    if "pass_filter" in filter_results.columns:
        native_pass = filter_results["pass_filter"].fillna(False).astype(bool)
    else:
        native_pass = pd.Series(True, index=filter_results.index, dtype=bool)
    publication_pass = native_pass & severity_pass
    res_df = _create_base_stats_df(
        model_name,
        num_mol,
        mean_severity=filter_results.severity.mean(),
        max_severity=filter_results.severity.max(),
        mean_n_covalent_motif=filter_results.n_covalent_motif.mean(),
        mean_nonzero_special_mol=(filter_results.special_mol > 0).mean(),
    )
    res_df["banned_ratio"] = 1 - publication_pass.mean()
    res_df["severity_threshold"] = threshold
    if "pass_filter" in filter_results.columns:
        res_df["native_pass_filter_banned_ratio"] = (
            1 - filter_results["pass_filter"].fillna(False).astype(bool).mean()
        )
    res_df, filter_extended = common_postprocessing_statistics(
        filter_results, res_df, stat, extend
    )
    if "pass_filter" in filter_extended.columns:
        filter_extended["native_pass_filter"] = filter_extended["pass_filter"]
    extended_severity = pd.to_numeric(filter_extended["severity"], errors="coerce")
    extended_severity_pass = extended_severity.lt(threshold).fillna(False)
    if "pass_filter" in filter_extended.columns:
        extended_native_pass = filter_extended["pass_filter"].fillna(False).astype(bool)
    else:
        extended_native_pass = pd.Series(True, index=filter_extended.index, dtype=bool)
    filter_extended["pass"] = extended_native_pass & extended_severity_pass
    return res_df, filter_extended


def _stats_lilly(config, filter_results, model_name, num_mol, stat, extend):
    """Compute statistics for lilly filter."""
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


_STATS_REGISTRY = {
    "common_alerts": _stats_common_alerts,
    "molgraph_stats": _stats_molgraph,
    "molcomplexity": _stats_molcomplexity,
    "bredt": _stats_simple_banned,
    "protecting_groups": _stats_simple_banned,
    "ring_infraction": _stats_simple_banned,
    "stereo_center": _stats_simple_banned,
    "halogenicity": _stats_simple_banned,
    "symmetry": _stats_simple_banned,
    "NIBR": _stats_nibr,
    "lilly": _stats_lilly,
}


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

    handler = _STATS_REGISTRY.get(filter_name)
    if handler is None:
        raise ValueError(f"Filter {filter_name} not found")
    return handler(config, filter_results, model_name, num_mol, stat, extend)


def _pass_mask_from_column(
    filter_extended: pd.DataFrame, pass_column: str
) -> pd.DataFrame:
    """Build an identity-keyed pass mask from one calculated decision column."""
    identity_cols = [
        col
        for col in ("smiles", "model_name", "mol_idx")
        if col in filter_extended.columns
    ]
    if pass_column not in filter_extended.columns:
        return pd.DataFrame(columns=[*identity_cols, "pass"])
    result = filter_extended[identity_cols + [pass_column]].copy()
    if pass_column != "pass":
        result = result.rename(columns={pass_column: "pass"})
    result["pass"] = result["pass"].fillna(False).astype(bool)
    return result.drop_duplicates(subset=identity_cols, keep="last")


def build_structural_policy_pass_masks(
    config: dict,
    filter_name: str,
    filter_extended: pd.DataFrame,
    default_mask: pd.DataFrame | None = None,
) -> dict[str, pd.DataFrame]:
    """Build independently selectable hard-policy masks from one calculation."""
    masks = {
        filter_name: (
            default_mask
            if default_mask is not None
            else build_structural_enforcement_pass_mask(
                config, filter_name, filter_extended
            )
        )
    }
    undefined_policy_configured = (
        "filter_undefined_stereo_center" in config
        or "undefined_stereo_center" in (config.get("enforced_filters") or [])
    )
    if filter_name == "stereo_center" and undefined_policy_configured:
        masks["undefined_stereo_center"] = _pass_mask_from_column(
            filter_extended, "undefined_stereo_pass"
        )
    return masks


def merge_structural_policy_aliases(
    profile: pd.DataFrame, filter_name: str
) -> pd.DataFrame:
    """Expose policy-specific stereo decisions in the liability profile."""
    if filter_name != "stereo_center":
        return profile
    out = profile.copy()
    column_aliases = {
        "stereo_center__undefined_stereo_pass": ("undefined_stereo_center__pass"),
        "stereo_center__undefined_stereo_reason": ("undefined_stereo_center__reason"),
        "stereo_center__n_undefined_stereo_centers": ("undefined_stereo_center__count"),
        "stereo_center__stereo_max_undefined": ("undefined_stereo_center__maximum"),
    }
    for source, target in column_aliases.items():
        if source in out.columns:
            out[target] = out[source]
    pass_column = "undefined_stereo_center__pass"
    if pass_column in out.columns:
        out["undefined_stereo_center__enforcement_pass"] = (
            out[pass_column].fillna(False).astype(bool)
        )
    return out


def build_structural_enforcement_pass_mask(
    config: dict,
    filter_name: str,
    filter_extended: pd.DataFrame,
) -> pd.DataFrame:
    """Build the survival mask independently from diagnostic calculation."""
    identity_cols = [
        col
        for col in ("smiles", "model_name", "mol_idx")
        if col in filter_extended.columns
    ]
    if filter_name != "common_alerts":
        if "pass" in filter_extended.columns:
            pass_column = "pass"
        elif "pass_filter" in filter_extended.columns:
            pass_column = "pass_filter"
        else:
            return pd.DataFrame(columns=[*identity_cols, "pass"])
        result = filter_extended[identity_cols + [pass_column]].copy()
        if pass_column != "pass":
            result = result.rename(columns={pass_column: "pass"})
        result["pass"] = result["pass"].fillna(False).astype(bool)
        return result.drop_duplicates(subset=identity_cols, keep="last")

    available_rulesets = {
        column.removeprefix("pass_")
        for column in filter_extended.columns
        if column.startswith("pass_") and column != "pass_any"
    }
    include_rulesets = {
        str(value)
        for value in config.get("common_alerts_filter_include_rulesets", []) or []
    }
    exclude_rulesets = {
        str(value)
        for value in config.get("common_alerts_filter_exclude_rulesets", []) or []
    }
    unknown_rulesets = (include_rulesets | exclude_rulesets) - available_rulesets
    if unknown_rulesets:
        names = ", ".join(sorted(unknown_rulesets))
        raise ValueError(f"Unknown Common Alerts filter ruleset(s): {names}")

    if include_rulesets:
        selected_rulesets = include_rulesets - exclude_rulesets
    else:
        selected_rulesets = available_rulesets - exclude_rulesets

    def _passes_selected_alerts(value) -> bool:
        records = _common_alert_profile_values(value)[0]
        for record in records:
            if not isinstance(record, dict):
                continue
            if record.get("ruleset") in selected_rulesets:
                return False
        return True

    result = filter_extended[identity_cols].copy()
    alert_json = filter_extended.get(
        "alert_hits_json", pd.Series("[]", index=filter_extended.index)
    )
    result["pass"] = alert_json.map(_passes_selected_alerts).astype(bool)
    return result.drop_duplicates(subset=identity_cols, keep="last")


def attach_structural_enforcement_pass(
    config: dict,
    filter_name: str,
    filter_extended: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Attach the configured survival decision without replacing native pass."""
    mask = build_structural_enforcement_pass_mask(config, filter_name, filter_extended)
    identity_cols = [
        col
        for col in ("smiles", "model_name", "mol_idx")
        if col in filter_extended.columns and col in mask.columns
    ]
    if not identity_cols or mask.empty:
        result = filter_extended.copy()
        result["enforcement_pass"] = False
        return result, mask
    enforcement = mask.rename(columns={"pass": "enforcement_pass"})
    result = filter_extended.drop(columns=["enforcement_pass"], errors="ignore").merge(
        enforcement, on=identity_cols, how="left"
    )
    result["enforcement_pass"] = result["enforcement_pass"].fillna(False).astype(bool)
    return result, mask


def initialize_structural_liability_profile(input_df: pd.DataFrame) -> pd.DataFrame:
    """Create the identity spine for one molecule-level structural profile."""
    identity_cols = [
        col for col in ("smiles", "model_name", "mol_idx") if col in input_df
    ]
    return input_df[identity_cols].drop_duplicates().copy()


def merge_structural_liability_profile(
    profile: pd.DataFrame,
    filter_name: str,
    filter_extended: pd.DataFrame,
) -> pd.DataFrame:
    """Add every row-level output from one structural filter to the profile."""
    identity_cols = [
        col
        for col in ("smiles", "model_name", "mol_idx")
        if col in profile.columns and col in filter_extended.columns
    ]
    if not identity_cols:
        return profile
    details = filter_extended.drop(columns=["mol"], errors="ignore").copy()
    details = details.drop_duplicates(subset=identity_cols, keep="last")
    details = details.rename(
        columns={
            col: f"{filter_name}__{col}"
            for col in details.columns
            if col not in identity_cols
        }
    )
    return profile.merge(details, on=identity_cols, how="left")


def _profile_pass(profile: pd.DataFrame, filter_name: str) -> pd.Series:
    column = f"{filter_name}__pass"
    if column not in profile.columns:
        return pd.Series(False, index=profile.index, dtype=bool)
    return profile[column].fillna(False).astype(bool)


def _profile_enforcement_pass(profile: pd.DataFrame, filter_name: str) -> pd.Series:
    column = f"{filter_name}__enforcement_pass"
    if column in profile.columns:
        return profile[column].fillna(False).astype(bool)
    return _profile_pass(profile, filter_name)


def _common_alert_profile_values(value):
    """Summarize atom-level Common Alerts JSON without discarding raw hits."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        records = []
    else:
        try:
            records = json.loads(str(value))
        except (TypeError, ValueError, json.JSONDecodeError):
            records = []
    if not isinstance(records, list):
        records = []
    rule_ids = sorted(
        {
            str(record.get("rule_id"))
            for record in records
            if isinstance(record, dict) and record.get("rule_id") is not None
        }
    )
    return records, len(records), len(rule_ids), ";".join(rule_ids)


def finalize_structural_liability_profile(
    profile: pd.DataFrame,
    enforced_filters,
    calculated_filters,
) -> pd.DataFrame:
    """Add hard-pass and warning summaries while retaining all raw diagnostics."""
    out = profile.copy()
    calculated = list(calculated_filters)
    hard_filters = (
        calculated
        if enforced_filters is None
        else [name for name in calculated if name in set(enforced_filters)]
    )
    diagnostic_filters = [name for name in calculated if name not in hard_filters]
    out["stage3_hard_filters"] = ";".join(hard_filters)
    hard_filter_set = set(hard_filters)
    for name in calculated:
        out[f"{name}__filter_enabled"] = name in hard_filter_set

    if hard_filters:
        hard_matrix = pd.concat(
            [
                _profile_enforcement_pass(out, name).rename(name)
                for name in hard_filters
            ],
            axis=1,
        )
        out["stage3_hard_pass"] = hard_matrix.all(axis=1)
        out["hard_failed_filters"] = hard_matrix.apply(
            lambda row: ";".join(name for name, passed in row.items() if not passed),
            axis=1,
        )
    else:
        out["stage3_hard_pass"] = True
        out["hard_failed_filters"] = ""

    if diagnostic_filters:
        diagnostic_matrix = pd.concat(
            [_profile_pass(out, name).rename(name) for name in diagnostic_filters],
            axis=1,
        )
        out["structural_warning_count"] = (~diagnostic_matrix).sum(axis=1)
        out["diagnostic_failed_filters"] = diagnostic_matrix.apply(
            lambda row: ";".join(name for name, passed in row.items() if not passed),
            axis=1,
        )
    else:
        out["structural_warning_count"] = 0
        out["diagnostic_failed_filters"] = ""

    ruleset_columns = [
        col
        for col in out.columns
        if col.startswith("common_alerts__pass_")
        and col not in {"common_alerts__pass_any"}
    ]
    out["common_alert_ruleset_hit_count"] = (
        (~out[ruleset_columns].fillna(False).astype(bool)).sum(axis=1)
        if ruleset_columns
        else 0
    )
    alert_json_column = "common_alerts__alert_hits_json"
    if alert_json_column in out.columns:
        alert_values = out[alert_json_column].map(_common_alert_profile_values)
        out["common_alert_match_count"] = alert_values.map(lambda value: value[1])
        out["common_alert_unique_rule_count"] = alert_values.map(lambda value: value[2])
        out["common_alert_rule_ids"] = alert_values.map(lambda value: value[3])
        for pass_column in ruleset_columns:
            ruleset = pass_column.removeprefix("common_alerts__pass_")
            prefix = f"common_alert_ruleset__{ruleset}"
            out[f"{prefix}__match_count"] = alert_values.map(
                lambda value, name=ruleset: sum(
                    1
                    for record in value[0]
                    if isinstance(record, dict) and record.get("ruleset") == name
                )
            )
            out[f"{prefix}__unique_rule_count"] = alert_values.map(
                lambda value, name=ruleset: len(
                    {
                        str(record.get("rule_id"))
                        for record in value[0]
                        if isinstance(record, dict)
                        and record.get("ruleset") == name
                        and record.get("rule_id") is not None
                    }
                )
            )
    else:
        out["common_alert_match_count"] = 0
        out["common_alert_unique_rule_count"] = 0
        out["common_alert_rule_ids"] = ""

    if {"NIBR", "molgraph_stats"}.issubset(calculated):
        out["nibr_molgraph_policy_pass"] = _profile_pass(out, "NIBR") & _profile_pass(
            out, "molgraph_stats"
        )
    if {"lilly", "molgraph_stats"}.issubset(calculated):
        out["lilly160_molgraph_policy_pass"] = _profile_pass(
            out, "lilly"
        ) & _profile_pass(out, "molgraph_stats")
    return out


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
    folder_to_save = Path(process_path(config[KEY_FOLDER_TO_SAVE]))
    config_sf = load_config(config[CFG_STRUCT_FILTERS])

    struct_folder = folder_to_save / stage_dir
    paths = list(struct_folder.glob("*/metrics.csv"))
    if not paths:
        paths = list(struct_folder.glob("*metrics.csv"))
    paths = [str(p) for p in paths]
    check_paths(config_sf, paths)

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
    folder_to_save = Path(process_path(config[KEY_FOLDER_TO_SAVE]))
    folder_name = config[KEY_FOLDER_TO_SAVE].split("/")[-1]

    config_sf = load_config(config[CFG_STRUCT_FILTERS])

    struct_folder = folder_to_save / stage_dir
    paths = list(struct_folder.glob("*/metrics.csv"))
    if not paths:
        paths = list(struct_folder.glob("*metrics.csv"))
    paths = [str(p) for p in paths]
    check_paths(config_sf, paths)

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

    identity_cols = [
        c for c in ["smiles", "model_name", "mol_idx"] if c in input_df.columns
    ]
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
        local_df = local_df[cols_to_take].drop_duplicates(
            subset=merge_cols, keep="last"
        )
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
    base_folder = Path(process_path(config[KEY_FOLDER_TO_SAVE]))
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
        descriptor_candidates = [
            base_folder
            / "stages"
            / "02_descriptors_initial"
            / "filtered_molecules.csv",
            base_folder
            / "stages"
            / "02_descriptors_initial"
            / "filtered"
            / "filtered_molecules.csv",
            base_folder / "Descriptors" / "passDescriptorsSMILES.csv",
        ]
        descriptors_path = next(
            (path for path in descriptor_candidates if path.exists()), None
        )
        if descriptors_path is not None:
            input_path = str(descriptors_path)
        else:
            sampled_path = base_folder / "sampled_molecules.csv"
            if sampled_path.exists():
                input_path = str(sampled_path)
            else:
                try:
                    from hedgehog.struct_filters.main import _get_input_path

                    input_path = _get_input_path(config, stage_dir, str(base_folder))
                except Exception:
                    input_path = None
    else:
        sampled_path = base_folder / "sampled_molecules.csv"
        if sampled_path.exists():
            input_path = str(sampled_path)
        else:
            try:
                from hedgehog.struct_filters.main import _get_input_path

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


def _parse_failure_reasons(df, col, filter_name):
    """Parse semicolon-delimited reasons from failed molecules for one filter column.

    Returns (sorted_reasons_list, reason_counts_dict).
    """
    reasons_col = f"reasons_{filter_name}"
    if reasons_col not in df.columns:
        return [], {}

    failed_molecules = df[~df[col]]
    reasons_data = failed_molecules[reasons_col].dropna()

    reason_counts: dict[str, int] = {}
    for reasons_str in reasons_data:
        if not (pd.notna(reasons_str) and str(reasons_str).strip()):
            continue
        for reason in str(reasons_str).split(";"):
            reason = reason.strip()
            if reason:
                reason_counts[reason] = reason_counts.get(reason, 0) + 1

    sorted_reasons = sorted(reason_counts.items(), key=lambda x: x[1], reverse=True)
    return sorted_reasons, reason_counts


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

        sorted_reasons, reason_counts = _parse_failure_reasons(df, col, filter_name)
        if sorted_reasons:
            filter_reasons[filter_name] = sorted_reasons
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
    folder_to_save = Path(process_path(config[KEY_FOLDER_TO_SAVE]))
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
