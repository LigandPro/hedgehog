"""MolEval generative metrics integration for HEDGEHOG reports.

Computes intrinsic distribution quality metrics (MOSES-style) at multiple
pipeline stages using the vendored moleval package from MolScore v1.9.5.
"""

import logging
import random
import statistics
import warnings
from collections.abc import Callable
from typing import Any

from hedgehog.utils.datamol_import import import_datamol_quietly
from hedgehog.utils.mce18 import compute_mce18
from hedgehog.utils.parallel import resolve_n_jobs

dm = import_datamol_quietly()

logger = logging.getLogger(__name__)

# Metric key -> config flag mapping
_METRIC_CONFIG_MAP: dict[str, str] = {
    "Validity": "validity",
    "Uniqueness": "uniqueness",
    "IntDiv1": "internal_diversity",
    "IntDiv2": "internal_diversity",
    "SEDiv": "se_diversity",
    "ScaffDiv": "scaffold_diversity",
    "ScaffUniqueness": "scaffold_diversity",
    "FG": "functional_groups",
    "RS": "ring_systems",
    "Filters": "filters",
    "MCE18": "mce18",
}

# Intrinsic metrics to extract (excludes count-based keys like '#', '# valid', etc.)
_INTRINSIC_KEYS = list(_METRIC_CONFIG_MAP.keys())


def _close_getmetrics_pool(gm: Any) -> None:
    """Close a GetMetrics pool robustly across vendor variants.

    Some vendored versions expose ``close_pool`` as an instance method,
    while others shadow it with a boolean flag on the instance. This helper
    supports both layouts and avoids leaking worker pools.
    """
    close_attr = getattr(gm, "close_pool", None)
    if callable(close_attr):
        close_attr()
        return

    close_method = getattr(type(gm), "close_pool", None)
    if callable(close_method):
        close_method(gm)


def is_available() -> bool:
    """Check if the vendored moleval package is importable."""
    try:
        from hedgehog.vendor.moleval.metrics.metrics import GetMetrics  # noqa: F401

        return True
    except ImportError:
        return False


def compute_stage_metrics(
    stage_smiles: dict[str, list[str]],
    config: dict[str, Any],
    seed: int = 42,
    progress_cb: Callable[[int, int, str | None], None] | None = None,
) -> dict[str, Any]:
    """Compute intrinsic generative metrics for each pipeline stage.

    Args:
        stage_smiles: Mapping of stage name to list of SMILES strings.
        config: MolEval configuration dict (from config_moleval.yml).

    Returns:
        Dictionary with keys 'by_stage', 'stages', and 'metrics'.
        Empty dict if run=false or no data available.
    """
    if not config.get("run", True):
        return {}

    if not stage_smiles:
        return {}

    if not is_available():
        logger.warning("MolEval not available, skipping generative metrics")
        return {}

    rng = random.Random(seed)

    from hedgehog.vendor.moleval.metrics.metrics import GetMetrics

    n_jobs = resolve_n_jobs(stage_config=config, default=-1)
    device = config.get("device", "cpu")
    max_molecules = config.get("max_molecules", 2000)
    # Avoid spawning unnecessary processes for tiny datasets.
    # This keeps the report phase lightweight and avoids exhausting file descriptors
    # on systems with conservative per-process limits.
    n_jobs = max(1, n_jobs)

    stage_items = [(name, smi) for name, smi in stage_smiles.items() if smi]
    if not stage_items:
        return {}

    total_stages = len(stage_items)
    if progress_cb is not None:
        progress_cb(0, total_stages, "MolEval")

    by_stage: dict[str, dict[str, float]] = {}
    all_metrics: set[str] = set()
    gm_by_jobs: dict[int, Any] = {}
    mce18_cache: dict[str, float | None] = {}

    try:
        for idx, (stage_name, smiles_list) in enumerate(stage_items, start=1):
            # Subsample if exceeding max_molecules
            if len(smiles_list) > max_molecules:
                smiles_list = rng.sample(smiles_list, max_molecules)

            stage_jobs = max(1, min(n_jobs, len(smiles_list)))
            gm = gm_by_jobs.get(stage_jobs)
            if gm is None:
                gm = GetMetrics(n_jobs=stage_jobs, device=device, run_fcd=False)
                gm_by_jobs[stage_jobs] = gm

            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    raw_metrics = gm.calculate(
                        gen=smiles_list,
                        calc_valid=True,
                        calc_unique=True,
                    )
            except Exception as e:
                logger.warning("MolEval failed for stage %s: %s", stage_name, e)
                if progress_cb is not None:
                    progress_cb(idx, total_stages, f"MolEval failed: {stage_name}")
                continue

            filtered = _filter_metrics(raw_metrics, config)

            # Add MCE-18 if enabled (computed separately, not from GetMetrics)
            if config.get("mce18", True):
                mce18_values = []
                for smi in smiles_list:
                    if smi not in mce18_cache:
                        mol = dm.to_mol(smi)
                        mce18_cache[smi] = (
                            compute_mce18(mol) if mol is not None else None
                        )
                    cached_val = mce18_cache[smi]
                    if cached_val is not None:
                        mce18_values.append(cached_val)
                if mce18_values:
                    filtered["MCE18"] = round(statistics.mean(mce18_values), 4)

            if filtered:
                by_stage[stage_name] = filtered
                all_metrics.update(filtered.keys())

            if progress_cb is not None:
                progress_cb(idx, total_stages, f"MolEval: {stage_name}")
    finally:
        for gm in gm_by_jobs.values():
            try:
                _close_getmetrics_pool(gm)
            except Exception as e:
                logger.debug("Failed to close MolEval process pool: %s", e)

    if not by_stage:
        return {}

    stages = list(by_stage.keys())
    metrics = sorted(all_metrics)

    return {
        "by_stage": by_stage,
        "stages": stages,
        "metrics": metrics,
    }


def _filter_metrics(
    raw: dict[str, Any],
    config: dict[str, Any],
) -> dict[str, float]:
    """Filter raw metrics to only those enabled in config.

    Normalizes variant keys like 'SEDiv@1k' -> 'SEDiv' for consistency.

    Args:
        raw: Raw metrics dict from GetMetrics.calculate().
        config: MolEval configuration dict.

    Returns:
        Filtered dict of metric_key -> float value.
    """
    # Normalize variant keys (SEDiv@1k -> SEDiv, SPDiv@1k -> SPDiv)
    normalized = dict(raw)
    for key in list(normalized.keys()):
        if key.startswith("SEDiv@"):
            normalized["SEDiv"] = normalized.pop(key)
        elif key.startswith("SPDiv@"):
            normalized["SPDiv"] = normalized.pop(key)

    result: dict[str, float] = {}
    for key in _INTRINSIC_KEYS:
        if key not in normalized:
            continue
        config_flag = _METRIC_CONFIG_MAP.get(key, "")
        if not config.get(config_flag, True):
            continue
        value = normalized[key]
        if isinstance(value, (int, float)):
            result[key] = round(float(value), 4)
    return result
