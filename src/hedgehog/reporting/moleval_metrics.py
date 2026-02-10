"""MolEval generative metrics integration for HEDGEHOG reports.

Computes intrinsic distribution quality metrics (MOSES-style) at multiple
pipeline stages using the vendored moleval package from MolScore v1.9.5.
"""

import logging
import random
import warnings
from typing import Any

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

    n_jobs = config.get("n_jobs", 1)
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


# =========================================================================
# Drug-likeness Compliance (rules_filter)
# =========================================================================

# Default rule sets for compliance reporting
_DEFAULT_RULE_SETS = ["rule_of_five", "rule_of_veber", "rule_of_generative_design"]

# Display names for rules
_RULE_DISPLAY_NAMES: dict[str, str] = {
    "rule_of_five": "Lipinski",
    "rule_of_veber": "Veber",
    "rule_of_ghose": "Ghose",
    "rule_of_egan": "Egan",
    "rule_of_generative_design": "GenDesign",
    "rule_of_generative_design_strict": "GenDesign(strict)",
    "rule_of_three": "Ro3",
    "rule_of_reos": "REOS",
    "rule_of_cns": "CNS",
    "rule_of_oprea": "Oprea",
    "rule_of_xu": "Xu",
    "rule_of_druglike_soft": "DrugLike",
    "rule_of_leadlike_soft": "LeadLike",
}

# Default chemical groups for monitoring
_DEFAULT_CHEMICAL_GROUPS = ["rings_in_drugs", "privileged_scaffolds"]


def compute_rules_compliance(
    stage_smiles: dict[str, list[str]],
    config: dict[str, Any],
    seed: int = 42,
) -> dict[str, Any]:
    """Compute drug-likeness rule compliance rates across pipeline stages.

    For each stage and each rule set, calculates the fraction of molecules
    that pass the rule.  Shows how compliance evolves through the pipeline.

    Args:
        stage_smiles: Mapping of stage name to list of SMILES strings.
        config: MolEval configuration dict (from config_moleval.yml).

    Returns:
        Dictionary with keys 'by_stage', 'stages', 'rules', or empty dict.
    """
    if not config.get("rules_compliance", False):
        return {}

    if not stage_smiles:
        return {}

    try:
        import datamol as dm
        import medchem as mc
    except ImportError:
        logger.warning("medchem/datamol not available, skipping rules compliance")
        return {}

    rng = random.Random(seed)
    rule_names = config.get("rules_sets", _DEFAULT_RULE_SETS)
    max_molecules = config.get("max_molecules", 2000)

    by_stage: dict[str, dict[str, float]] = {}

    for stage_name, smiles_list in stage_smiles.items():
        if not smiles_list:
            continue

        if len(smiles_list) > max_molecules:
            smiles_list = rng.sample(smiles_list, max_molecules)

        mols = [dm.to_mol(s) for s in smiles_list]
        mols = [m for m in mols if m is not None]
        if not mols:
            continue

        stage_results: dict[str, float] = {}
        for rule_name in rule_names:
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    passed = mc.functional.rules_filter(
                        mols=mols,
                        rules=[rule_name],
                        return_idx=False,
                        n_jobs=1,
                        progress=False,
                    )
                display = _RULE_DISPLAY_NAMES.get(rule_name, rule_name)
                stage_results[display] = round(float(passed.mean()), 4)
            except Exception as e:
                logger.debug(
                    "rules_filter failed for %s/%s: %s", stage_name, rule_name, e
                )

        if stage_results:
            by_stage[stage_name] = stage_results

    if not by_stage:
        return {}

    stages = list(by_stage.keys())
    rules = sorted({r for vals in by_stage.values() for r in vals})

    return {
        "by_stage": by_stage,
        "stages": stages,
        "rules": rules,
    }


def compute_group_match_rates(
    stage_smiles: dict[str, list[str]],
    config: dict[str, Any],
    seed: int = 42,
) -> dict[str, Any]:
    """Compute chemical group match rates across pipeline stages.

    For each stage and each chemical group, calculates the fraction of
    molecules that contain the group (drug-relevant scaffolds).

    Args:
        stage_smiles: Mapping of stage name to list of SMILES strings.
        config: MolEval configuration dict (from config_moleval.yml).

    Returns:
        Dictionary with keys 'by_stage', 'stages', 'groups', or empty dict.
    """
    if not config.get("chemical_groups", False):
        return {}

    if not stage_smiles:
        return {}

    try:
        import datamol as dm
        import medchem as mc
        from medchem.groups import ChemicalGroup
    except ImportError:
        logger.warning("medchem/datamol not available, skipping chemical groups")
        return {}

    rng = random.Random(seed)
    group_names = config.get("chemical_group_names", _DEFAULT_CHEMICAL_GROUPS)
    max_molecules = config.get("max_molecules", 2000)

    by_stage: dict[str, dict[str, float]] = {}

    for stage_name, smiles_list in stage_smiles.items():
        if not smiles_list:
            continue

        if len(smiles_list) > max_molecules:
            smiles_list = rng.sample(smiles_list, max_molecules)

        mols = [dm.to_mol(s) for s in smiles_list]
        mols = [m for m in mols if m is not None]
        if not mols:
            continue

        stage_results: dict[str, float] = {}
        for group_name in group_names:
            try:
                cg = ChemicalGroup(group_name)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    # chemical_group_filter returns True = does NOT contain group
                    # We want match rate = fraction that DO contain the group
                    passed = mc.functional.chemical_group_filter(
                        mols=mols,
                        chemical_group=cg,
                        return_idx=False,
                        n_jobs=1,
                        progress=False,
                    )
                match_rate = 1.0 - float(passed.mean())
                stage_results[group_name] = round(match_rate, 4)
            except Exception as e:
                logger.debug(
                    "chemical_group_filter failed for %s/%s: %s",
                    stage_name,
                    group_name,
                    e,
                )

        if stage_results:
            by_stage[stage_name] = stage_results

    if not by_stage:
        return {}

    stages = list(by_stage.keys())
    groups = sorted({g for vals in by_stage.values() for g in vals})

    return {
        "by_stage": by_stage,
        "stages": stages,
        "groups": groups,
    }
