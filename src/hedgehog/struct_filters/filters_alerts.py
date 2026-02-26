"""Common structural alerts filter (SMARTS-based screening)."""

import threading

import pandas as pd
from rdkit import Chem

from hedgehog._constants import CFG_STRUCT_FILTERS
from hedgehog.configs.logger import load_config, logger
from hedgehog.struct_filters._helpers import (
    add_model_name_col,
    dm,
    filter_alerts,
    _silence_worker_stdio,
)
from hedgehog.utils.parallel import parallel_map, resolve_n_jobs

# Module-level globals for worker processes
_ALERT_COMPILED_SMARTS: list[tuple[str, object, str]] | None = None
_ALERT_RULESET_NAMES: list[str] | None = None


def _check_alerts_single_mol(args):
    """Check all alert SMARTS patterns against a single molecule."""
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


def _compile_alert_smarts(alert_data) -> list[tuple[str, object, str]]:
    """Pre-compile SMARTS patterns from alert data to avoid per-molecule recompilation."""
    compiled: list[tuple[str, object, str]] = []
    for _, row in alert_data.iterrows():
        smarts = row.get("smarts")
        ruleset = row.get("rule_set_name")
        desc = row.get("description") or ""
        if not isinstance(smarts, str) or not isinstance(ruleset, str):
            continue
        patt = Chem.MolFromSmarts(smarts)
        if patt is not None:
            compiled.append((ruleset, patt, str(desc)))
    return compiled


def _setup_alerts_progress(total_items, progress_cb):
    """Set up progress tracking and heartbeat thread for structural alerts."""
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
    for row_data in results_list:
        idx = row_data.pop("_mol_idx")
        row_data["mol"] = mols[idx]

    results = pd.DataFrame(results_list)

    pass_cols = [
        c
        for c in results.columns
        if c.startswith("pass_") and c not in {"pass", "pass_any"}
    ]
    results["pass"] = results[pass_cols].all(axis=1)
    results["pass_any"] = results[pass_cols].any(axis=1)
    return results


def _resolve_common_alerts_n_jobs(config_sf: dict, config: dict, total_items: int) -> int:
    """Resolve worker count for Common Alerts using size-aware defaults.

    Policy:
    - total_items < 1_000  -> 1 worker
    - total_items < 10_000 -> 12 workers
    - otherwise            -> all available workers
    """
    auto_n_jobs = bool(config_sf.get("common_alerts_auto_n_jobs", True))
    if not auto_n_jobs:
        return resolve_n_jobs(config_sf, config)

    max_workers = resolve_n_jobs({"n_jobs": -1}, config)

    small_threshold_raw = config_sf.get("common_alerts_small_input_threshold", 1000)
    small_n_jobs_raw = config_sf.get("common_alerts_small_input_n_jobs", 1)
    medium_n_jobs_raw = config_sf.get("common_alerts_large_input_n_jobs", 12)

    try:
        small_threshold = int(small_threshold_raw)
    except (TypeError, ValueError):
        logger.warning(
            "Invalid common_alerts_small_input_threshold=%r; using 1000.",
            small_threshold_raw,
        )
        small_threshold = 1000

    try:
        small_n_jobs = int(small_n_jobs_raw)
    except (TypeError, ValueError):
        logger.warning(
            "Invalid common_alerts_small_input_n_jobs=%r; using 1.",
            small_n_jobs_raw,
        )
        small_n_jobs = 1

    try:
        medium_n_jobs = int(medium_n_jobs_raw)
    except (TypeError, ValueError):
        logger.warning(
            "Invalid common_alerts_large_input_n_jobs=%r; using 12.",
            medium_n_jobs_raw,
        )
        medium_n_jobs = 12

    if total_items < small_threshold:
        n_jobs = small_n_jobs
    elif total_items < 10_000:
        n_jobs = medium_n_jobs
    else:
        n_jobs = max_workers

    return max(1, min(n_jobs, max_workers))


def apply_structural_alerts(
    config, mols, smiles_model_name_mols=None, progress_cb=None
):
    logger.info("Calculating Common Alerts...")

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
    n_jobs = _resolve_common_alerts_n_jobs(config_sf, config, total_items)
    logger.info("Common Alerts workers: %d", n_jobs)

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
