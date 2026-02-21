import contextlib
import importlib
import io
import os
import threading
import warnings
from pathlib import Path

import pandas as pd
from rdkit import Chem

from hedgehog.configs.logger import load_config, logger
from hedgehog.utils.datamol_import import import_datamol_quietly
from hedgehog.utils.parallel import parallel_map, resolve_n_jobs

warnings.filterwarnings("ignore", category=FutureWarning)

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
    except (ValueError, RuntimeError):
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
                except (ValueError, TypeError, RuntimeError):
                    passed = False
        else:
            try:
                passed = bool(complexity_filter(prepared_mol))
            except (ValueError, TypeError, RuntimeError):
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


def _init_alert_worker_quiet(
    compiled_smarts: list[tuple[str, object, str]],
    rule_set_names: list[str],
) -> None:
    _silence_worker_stdio()
    _init_alert_worker(compiled_smarts, rule_set_names)


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
    logger.info("Common Alerts workers: %d", n_jobs)
    items = [(i, mol) for i, mol in enumerate(mols)]
    total_items = len(items)

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
            progress=_progress_wrapper,
            initializer=worker_initializer,
            initargs=(compiled_smarts, rule_set_names),
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
    logger.info("MolGraph workers: %d", n_jobs)

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

    if smiles_modelName_mols is not None:
        final_result = add_model_name_col(final_result, smiles_modelName_mols)
    return final_result


def apply_bredt_filter(config, mols, smiles_modelName_mols=None):
    logger.info("Calculating Bredt filter...")
    config_structFilters = load_config(config["config_structFilters"])
    n_jobs = resolve_n_jobs(config_structFilters, config)
    scheduler = _resolve_scheduler(config_structFilters, "bredt_scheduler")
    logger.info("Bredt workers: %d", n_jobs)
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
    logger.info("Protecting Groups workers: %d", n_jobs)
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
    logger.info("Ring Infraction workers: %d", n_jobs)
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
    logger.info("Stereo Center workers: %d", n_jobs)
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
    logger.info("Halogenicity workers: %d", n_jobs)
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
    logger.info("Symmetry workers: %d", n_jobs)
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
    except Exception as batch_error:  # noqa: BLE001 — intentional: Lilly filter raises various internal errors
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
        except Exception:  # noqa: BLE001 — intentional: Lilly filter raises various internal errors
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
    logger.info("Lilly workers: %d", n_jobs)

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
            except (ValueError, TypeError):
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
            initializer=_silence_worker_stdio if n_workers > 1 else None,
        )
        rows = []
        for chunk_rows in chunk_results:
            rows.extend(chunk_rows)
        results = pd.DataFrame(rows)
    except Exception as error:  # noqa: BLE001 — intentional: parallel Lilly may fail in various ways
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
