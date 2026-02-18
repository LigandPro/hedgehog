"""MolEval generative metrics integration for HEDGEHOG reports.

Computes intrinsic distribution quality metrics (MOSES-style) at multiple
pipeline stages using the vendored moleval package from MolScore v1.9.5.
"""

import logging
import random
import statistics
import warnings
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

    gm = GetMetrics(n_jobs=n_jobs, device=device, run_fcd=False)

    by_stage: dict[str, dict[str, float]] = {}
    all_metrics: set[str] = set()

    for stage_name, smiles_list in stage_smiles.items():
        if not smiles_list:
            continue

        # Subsample if exceeding max_molecules
        if len(smiles_list) > max_molecules:
            smiles_list = rng.sample(smiles_list, max_molecules)

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
            continue

        filtered = _filter_metrics(raw_metrics, config)

        # Add MCE-18 if enabled (computed separately, not from GetMetrics)
        if config.get("mce18", True):
            mce18_values = []
            for smi in smiles_list:
                mol = dm.to_mol(smi)
                if mol is not None:
                    val = compute_mce18(mol)
                    if val is not None:
                        mce18_values.append(val)
            if mce18_values:
                filtered["MCE18"] = round(statistics.mean(mce18_values), 4)

        if filtered:
            by_stage[stage_name] = filtered
            all_metrics.update(filtered.keys())

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
