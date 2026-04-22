"""Generator reality scorecard for HEDGEHOG reports."""

from __future__ import annotations

import math
import statistics
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pandas as pd

COMPONENT_ORDER = [
    "yield",
    "physchem",
    "structural",
    "synthesis",
    "docking_pose",
    "diversity",
]

STRUCTURAL_PASS_COLUMNS = [
    "pass_common_alerts",
    "pass_molgraph_stats",
    "pass_molcomplexity",
    "pass_NIBR",
    "pass_bredt",
    "pass_lilly",
    "pass_protecting_groups",
    "pass_ring_infraction",
    "pass_stereo_center",
    "pass_halogenicity",
]

DESCRIPTOR_SUMMARY_COLUMNS = {
    "QED": "qed",
    "MolWt": "molWt",
    "LogP": "logP",
    "TPSA": "tpsa",
    "Fsp3": "fsp3",
}

DOCKING_AFFINITY_COLUMNS = [
    "gnina_affinity",
    "gnina_minimizedAffinity",
    "minimizedAffinity",
    "affinity",
]
DOCKING_CNNSCORE_COLUMNS = ["gnina_cnnscore", "gnina_CNNscore", "CNNscore"]
DOCKING_CNNAFFINITY_COLUMNS = [
    "gnina_cnnaffinity",
    "gnina_CNNaffinity",
    "CNNaffinity",
]
DOCKING_PASS_COLUMNS = [
    "pass_search_box",
    "pass_pose_quality",
    "pass_interactions",
    "pass_conformer_deviation",
    "pass",
]

DIVERSITY_METRICS = [
    "IntDiv1",
    "IntDiv2",
    "ScaffDiv",
    "ScaffUniqueness",
    "SEDiv",
]


def clip01(x: float) -> float:
    """Clamp a numeric value to the 0..1 interval."""
    if math.isnan(x):
        return 0.0
    return max(0.0, min(1.0, x))


def linear_score(
    value: float,
    bad: float,
    good: float,
    higher_is_better: bool,
) -> float:
    """Normalize a value between bad and good thresholds on a 0..1 scale."""
    if good == bad:
        return 1.0 if value >= good else 0.0
    if higher_is_better:
        return clip01((value - bad) / (good - bad))
    return clip01((bad - value) / (bad - good))


def sigmoid_score(value: float, midpoint: float, scale: float) -> float:
    """Map a value to 0..1 with a logistic curve."""
    if scale == 0:
        return 1.0 if value >= midpoint else 0.0
    exponent = -(value - midpoint) / scale
    return 1.0 / (1.0 + math.exp(exponent))


def safe_median(values: Sequence[float]) -> float | None:
    """Return the median for numeric non-null values, or None when empty."""
    clean = [float(v) for v in values if pd.notna(v)]
    if not clean:
        return None
    return float(statistics.median(clean))


def first_existing(base_path: Path, candidates: list[str]) -> Path | None:
    """Return the first existing path relative to base_path."""
    for candidate in candidates:
        path = base_path / candidate
        if path.exists():
            return path
    return None


def compute_weighted_scores(
    base_path: Path,
    available_models: list[str],
    config: dict[str, Any],
    moleval: dict[str, Any] | None = None,
    initial_count: int | None = None,
    final_count: int | None = None,
) -> dict[str, Any]:
    """Compute Generator Reality Assessment data for report_data.json."""
    if not config or config.get("run", True) is False:
        return {}

    base_path = Path(base_path)
    models = _resolve_models(base_path, available_models)
    if not models:
        return {}

    weights = _normalized_config_weights(config)
    candidate_pool_weights = _normalized_config_weights(
        config, key="candidate_pool_weights"
    )
    result = {
        "config": {
            "version": config.get("version", "v1"),
            "mode": config.get("mode", "generator_reality"),
            "weights": weights,
            "candidate_pool_weights": candidate_pool_weights,
        },
        "models": {},
        "ranking": [],
    }

    for model_name in models:
        evidence = collect_model_evidence(
            base_path=base_path,
            model_name=model_name,
            moleval=moleval or {},
            initial_count=initial_count,
            final_count=final_count,
        )
        scored = score_model_evidence(
            model_name,
            evidence,
            config,
            weights,
            candidate_pool_weights,
        )
        result["models"][model_name] = scored

    ranking = [
        {
            "model_name": model_name,
            "overall": score["overall"],
            "grade": score["grade"],
        }
        for model_name, score in result["models"].items()
        if score.get("overall") is not None
    ]
    result["ranking"] = sorted(ranking, key=lambda item: item["overall"], reverse=True)
    return result


def collect_model_evidence(
    base_path: Path,
    model_name: str,
    moleval: dict[str, Any] | None = None,
    initial_count: int | None = None,
    final_count: int | None = None,
) -> dict[str, Any]:
    """Collect all scorecard evidence for a model without scoring it."""
    return {
        "yield": collect_yield_evidence(
            base_path, model_name, initial_count, final_count
        ),
        "physchem": collect_physchem_evidence(base_path, model_name),
        "structural": collect_structural_evidence(base_path, model_name),
        "synthesis": collect_synthesis_evidence(base_path, model_name),
        "docking_pose": collect_docking_pose_evidence(base_path, model_name),
        "diversity": collect_diversity_evidence(moleval or {}),
    }


def collect_yield_evidence(
    base_path: Path,
    model_name: str,
    initial_count: int | None = None,
    final_count: int | None = None,
) -> dict[str, Any]:
    """Collect initial/final molecule counts for yield scoring."""
    warnings = []
    input_df = _read_first_csv(base_path, ["input/sampled_molecules.csv"], warnings)
    final_df = _read_first_csv(
        base_path,
        [
            "output/final_molecules.csv",
            "final_molecules.csv",
            "stages/07_descriptors_final/filtered/filtered_molecules.csv",
            "stages/06_docking_filters/filtered_molecules.csv",
        ],
        warnings,
    )

    n_initial = _filtered_count(input_df, model_name)
    n_final = _filtered_count(final_df, model_name)

    if n_initial is None and model_name == "__all__" and initial_count is not None:
        n_initial = initial_count
    if n_final is None and model_name == "__all__" and final_count is not None:
        n_final = final_count

    metrics = {
        "initial_count": int(n_initial or 0),
        "final_count": int(n_final or 0),
    }
    if metrics["initial_count"] > 0:
        metrics["final_retention_rate"] = (
            metrics["final_count"] / metrics["initial_count"]
        )
    return {
        "available": metrics["initial_count"] > 0 or metrics["final_count"] > 0,
        "metrics": metrics,
        "warnings": warnings,
    }


def collect_physchem_evidence(base_path: Path, model_name: str) -> dict[str, Any]:
    """Collect descriptor pass rates and key descriptor summaries."""
    warnings = []
    flags = _read_first_csv(
        base_path,
        [
            "stages/01_descriptors_initial/filtered/pass_flags.csv",
            "stages/07_descriptors_final/filtered/pass_flags.csv",
        ],
        warnings,
    )
    descriptors = _read_first_csv(
        base_path,
        [
            "stages/01_descriptors_initial/metrics/descriptors_all.csv",
            "stages/01_descriptors_initial/filtered/descriptors_passed.csv",
            "stages/07_descriptors_final/metrics/descriptors_all.csv",
            "stages/07_descriptors_final/filtered/descriptors_passed.csv",
        ],
        warnings,
    )

    flags = _filter_model(flags, model_name)
    descriptors = _filter_model(descriptors, model_name)
    pass_cols = _pass_columns(flags, suffix="_pass") if flags is not None else []
    pass_rates = {
        col: _mean_bool(flags[col])
        for col in pass_cols
        if flags is not None and not flags.empty
    }

    metrics: dict[str, Any] = {"pass_rates": pass_rates}
    if flags is not None and not flags.empty and pass_cols:
        bool_flags = flags[pass_cols].map(_truthy)
        metrics["all_pass_rate"] = float(bool_flags.all(axis=1).mean())
        metrics["all_pass_count"] = int(bool_flags.all(axis=1).sum())
        metrics["total_count"] = int(len(bool_flags))
        metrics["mean_flag_pass_rate"] = statistics.mean(pass_rates.values())
        worst_flag, worst_rate = min(pass_rates.items(), key=lambda item: item[1])
        metrics["worst_flag"] = worst_flag
        metrics["worst_flag_pass_rate"] = worst_rate
    if descriptors is not None and not descriptors.empty:
        col_map = {col.lower(): col for col in descriptors.columns}
        for display, metric_name in DESCRIPTOR_SUMMARY_COLUMNS.items():
            col = col_map.get(display.lower())
            if col:
                values = pd.to_numeric(descriptors[col], errors="coerce").dropna()
                if not values.empty:
                    metrics[f"mean_{metric_name}"] = float(values.mean())
                    metrics[f"median_{metric_name}"] = float(values.median())

    return {
        "available": bool(pass_rates),
        "metrics": metrics,
        "warnings": warnings,
    }


def collect_structural_evidence(base_path: Path, model_name: str) -> dict[str, Any]:
    """Collect structural filter pass rates."""
    warnings = []
    filtered = _read_first_csv(
        base_path,
        ["stages/03_structural_filters_post/filtered_molecules.csv"],
        warnings,
    )
    failed = _read_first_csv(
        base_path,
        ["stages/03_structural_filters_post/failed_molecules.csv"],
        warnings,
    )

    filtered = _filter_model(filtered, model_name)
    failed = _filter_model(failed, model_name)
    if filtered is not None and not filtered.empty and failed is not None:
        for col in STRUCTURAL_PASS_COLUMNS:
            if col in failed.columns and col not in filtered.columns:
                filtered = filtered.copy()
                filtered[col] = True
    frames = [df for df in [filtered, failed] if df is not None and not df.empty]
    combined = pd.concat(frames, ignore_index=True) if frames else None

    pass_rates = {}
    if combined is not None:
        for col in STRUCTURAL_PASS_COLUMNS:
            if col in combined.columns:
                pass_rates[col] = _mean_bool(combined[col])

    if not pass_rates and failed is not None and filtered is not None:
        total = len(failed) + len(filtered)
        if total:
            pass_rates["overall"] = len(filtered) / total

    filtered_count = int(len(filtered)) if filtered is not None else 0
    failed_count = int(len(failed)) if failed is not None else 0
    total_count = filtered_count + failed_count
    metrics = {
        "pass_rates": pass_rates,
        "filtered_count": filtered_count,
        "failed_count": failed_count,
    }
    if total_count > 0:
        metrics["stage_pass_rate"] = filtered_count / total_count
    if pass_rates:
        worst_filter, worst_rate = min(pass_rates.items(), key=lambda item: item[1])
        metrics["mean_filter_pass_rate"] = statistics.mean(pass_rates.values())
        metrics["worst_filter"] = worst_filter
        metrics["worst_filter_pass_rate"] = worst_rate

    return {
        "available": bool(pass_rates),
        "metrics": metrics,
        "warnings": warnings,
    }


def collect_synthesis_evidence(base_path: Path, model_name: str) -> dict[str, Any]:
    """Collect synthesis score evidence."""
    warnings = []
    df = _read_first_csv(
        base_path,
        [
            "stages/04_synthesis/synthesis_extended.csv",
            "stages/04_synthesis/synthesis_scores.csv",
            "stages/04_synthesis/filtered_molecules.csv",
        ],
        warnings,
    )
    df = _filter_model(df, model_name)
    if df is None or df.empty:
        return {"available": False, "metrics": {}, "warnings": warnings}

    metrics: dict[str, Any] = {"count": int(len(df))}
    for col, metric in [
        ("sa_score", "median_sa"),
        ("syba_score", "median_syba"),
        ("ra_score", "median_ra"),
        ("search_time", "median_search_time"),
    ]:
        if col in df.columns:
            median = safe_median(pd.to_numeric(df[col], errors="coerce").tolist())
            if median is not None:
                metrics[metric] = median

    if "solved" in df.columns:
        metrics["solve_rate"] = _mean_bool(df["solved"])

    required_any = {"median_sa", "median_syba", "median_ra", "solve_rate"}
    return {
        "available": any(key in metrics for key in required_any),
        "metrics": metrics,
        "warnings": warnings,
    }


def collect_docking_pose_evidence(base_path: Path, model_name: str) -> dict[str, Any]:
    """Collect docking score and pose-filter evidence."""
    warnings = []
    candidates = [
        "stages/06_docking_filters/metrics.csv",
        "stages/06_docking_filters/filtered_poses.csv",
        "output/final_molecules.csv",
        "final_molecules.csv",
    ]
    df = _read_first_csv_with_any_columns(
        base_path,
        candidates,
        [
            *DOCKING_AFFINITY_COLUMNS,
            *DOCKING_CNNSCORE_COLUMNS,
            *DOCKING_CNNAFFINITY_COLUMNS,
            *DOCKING_PASS_COLUMNS,
        ],
        warnings,
    )
    df = _filter_model(df, model_name)
    if df is None or df.empty:
        return {"available": False, "metrics": {}, "warnings": warnings}

    metrics: dict[str, Any] = {"count": int(len(df))}
    for columns, metric in [
        (DOCKING_AFFINITY_COLUMNS, "median_affinity"),
        (DOCKING_CNNSCORE_COLUMNS, "median_cnnscore"),
        (DOCKING_CNNAFFINITY_COLUMNS, "median_cnnaffinity"),
    ]:
        col = _first_column(df, columns)
        if col:
            median = safe_median(pd.to_numeric(df[col], errors="coerce").tolist())
            if median is not None:
                metrics[metric] = median

    pose_rates = {
        col: _mean_bool(df[col]) for col in DOCKING_PASS_COLUMNS if col in df.columns
    }
    if pose_rates:
        metrics["pose_pass_rates"] = pose_rates
        metrics["pose_pass_rate"] = statistics.mean(pose_rates.values())

    score_keys = {
        "median_affinity",
        "median_cnnscore",
        "median_cnnaffinity",
        "pose_pass_rate",
    }
    return {
        "available": any(key in metrics for key in score_keys),
        "metrics": metrics,
        "warnings": warnings,
    }


def collect_diversity_evidence(moleval: dict[str, Any]) -> dict[str, Any]:
    """Collect diversity metrics from already-computed MolEval output."""
    by_stage = (moleval or {}).get("by_stage", {})
    source_stage = "Input" if by_stage.get("Input") else "DockingFilters"
    stage_metrics = by_stage.get(source_stage, {})
    metrics = {
        metric: float(stage_metrics[metric])
        for metric in DIVERSITY_METRICS
        if stage_metrics.get(metric) is not None
    }
    if metrics:
        metrics["source_stage"] = source_stage
    return {
        "available": bool(metrics),
        "metrics": metrics,
        "warnings": [] if metrics else ["MolEval diversity metrics not available"],
    }


def score_model_evidence(
    model_name: str,
    evidence: dict[str, Any],
    config: dict[str, Any],
    weights: dict[str, float] | None = None,
    candidate_pool_weights: dict[str, float] | None = None,
) -> dict[str, Any]:
    """Score one model from pre-collected evidence."""
    weights = weights or _normalized_config_weights(config)
    candidate_pool_weights = candidate_pool_weights or _normalized_config_weights(
        config, key="candidate_pool_weights"
    )
    scorers = {
        "yield": score_yield,
        "physchem": score_physchem,
        "structural": score_structural,
        "synthesis": score_synthesis,
        "docking_pose": score_docking_pose,
        "diversity": score_diversity,
    }

    components = {}
    candidate_pool_components = {}
    warnings = []
    for component in COMPONENT_ORDER:
        component_evidence = evidence.get(component, {})
        component_warnings = component_evidence.get("warnings", [])
        warnings.extend(component_warnings)
        score = None
        candidate_pool_score = None
        if component_evidence.get("available"):
            score = scorers[component](component_evidence, config)
            candidate_pool_score = score_candidate_pool_component(
                component,
                component_evidence,
                config,
            )
        components[component] = {
            "score": round(score, 1) if score is not None else None,
            "weight": weights.get(component, 0.0),
            "available": score is not None,
            "evidence": component_evidence.get("metrics", {}),
            "evidence_text": _evidence_text(
                component, component_evidence.get("metrics", {})
            ),
        }
        candidate_pool_components[component] = {
            "score": (
                round(candidate_pool_score, 1)
                if candidate_pool_score is not None
                else None
            ),
            "weight": candidate_pool_weights.get(component, 0.0),
            "available": candidate_pool_score is not None,
        }

    uncapped_overall = _weighted_overall(components, weights)
    overall = uncapped_overall
    hard_caps_triggered = _applicable_hard_caps(evidence, config)
    hard_caps_applied = []
    if overall is not None and hard_caps_triggered:
        cap_value = min(cap["cap"] for cap in hard_caps_triggered)
        if cap_value < overall:
            overall = cap_value
            hard_caps_applied = hard_caps_triggered
            warnings.extend(cap["warning"] for cap in hard_caps_applied)
    candidate_pool_overall = _weighted_overall(
        candidate_pool_components,
        candidate_pool_weights,
    )
    final_count = int(
        evidence.get("yield", {}).get("metrics", {}).get("final_count", 0)
    )
    available_count = sum(1 for item in components.values() if item["available"])
    confidence = assign_confidence(final_count, available_count, config)
    bottlenecks = generate_bottlenecks(components)

    return {
        "overall": round(overall, 1) if overall is not None else None,
        "uncapped_overall": (
            round(uncapped_overall, 1) if uncapped_overall is not None else None
        ),
        "grade": assign_grade(overall) if overall is not None else "Weak",
        "confidence": confidence,
        "candidate_pool_quality": {
            "overall": (
                round(candidate_pool_overall, 1)
                if candidate_pool_overall is not None
                else None
            ),
            "grade": (
                assign_grade(candidate_pool_overall)
                if candidate_pool_overall is not None
                else "Weak"
            ),
            "components": candidate_pool_components,
        },
        "hard_caps_triggered": hard_caps_triggered,
        "hard_caps_applied": hard_caps_applied,
        "components": components,
        "bottlenecks": bottlenecks,
        "warnings": warnings,
        "model_name": model_name,
    }


def score_yield(evidence: dict[str, Any], config: dict[str, Any]) -> float:
    """Score final retention against the configured target on a 0..100 scale."""
    metrics = evidence.get("metrics", {})
    yield_cfg = config.get("yield", {})
    mode = str(yield_cfg.get("mode", "retention")).lower()
    if mode == "absolute":
        return score_yield_absolute(evidence, config)

    n_final = float(metrics.get("final_count") or 0)
    n_initial = float(metrics.get("initial_count") or 0)
    target = float(
        yield_cfg.get(
            "target_final_retention",
            config.get("target_final_retention", 0.10),
        )
        or 0.10
    )
    if n_initial <= 0 or target <= 0:
        return 0.0
    retention = n_final / n_initial
    return 100 * clip01(retention / target)


def score_yield_absolute(evidence: dict[str, Any], config: dict[str, Any]) -> float:
    """Score survivor pool size with the original count-saturation formula."""
    metrics = evidence.get("metrics", {})
    n_final = float(metrics.get("final_count") or 0)
    n_initial = float(metrics.get("initial_count") or 0)
    target = float(config.get("target_final_count", 100) or 100)
    yield_cfg = config.get("yield", {})
    count_weight = float(yield_cfg.get("count_weight", 0.70))
    retention_weight = float(yield_cfg.get("retention_weight", 0.30))
    total_weight = count_weight + retention_weight
    if total_weight <= 0:
        count_weight, retention_weight, total_weight = 0.70, 0.30, 1.0

    count_score = 1 - math.exp(-n_final / target) if target > 0 else 0.0
    if n_initial > 0:
        retention_score = math.log1p(n_final) / math.log1p(n_initial)
    else:
        retention_score = 0.0
    score = (
        count_weight * clip01(count_score) + retention_weight * clip01(retention_score)
    ) / total_weight
    return 100 * clip01(score)


def score_physchem(evidence: dict[str, Any], _config: dict[str, Any]) -> float:
    """Score descriptor compatibility by all-pass gate survival."""
    metrics = evidence.get("metrics", {})
    if metrics.get("all_pass_rate") is not None:
        return 100 * clip01(float(metrics["all_pass_rate"]))
    return _mean_score(metrics.get("pass_rates", {}).values())


def score_structural(evidence: dict[str, Any], config: dict[str, Any]) -> float:
    """Score structural compatibility by stage survival and worst gate."""
    metrics = evidence.get("metrics", {})
    if metrics.get("stage_pass_rate") is None:
        return _mean_score(metrics.get("pass_rates", {}).values())

    structural_cfg = config.get("structural", {})
    stage_weight = float(structural_cfg.get("stage_pass_weight", 0.80))
    worst_weight = float(structural_cfg.get("worst_filter_weight", 0.20))
    stage_pass = clip01(float(metrics["stage_pass_rate"]))
    worst_pass = clip01(float(metrics.get("worst_filter_pass_rate", stage_pass)))
    return _weighted_parts_score(
        [
            (stage_weight, stage_pass),
            (worst_weight, worst_pass),
        ]
    )


def score_candidate_pool_component(
    component: str,
    evidence: dict[str, Any],
    config: dict[str, Any],
) -> float:
    """Score survivor-pool quality with pre-v2 partial-pass semantics."""
    metrics = evidence.get("metrics", {})
    if component == "yield":
        return score_yield_absolute(evidence, config)
    if component == "physchem":
        if metrics.get("mean_flag_pass_rate") is not None:
            return 100 * clip01(float(metrics["mean_flag_pass_rate"]))
        return _mean_score(metrics.get("pass_rates", {}).values())
    if component == "structural":
        if metrics.get("mean_filter_pass_rate") is not None:
            return 100 * clip01(float(metrics["mean_filter_pass_rate"]))
        return _mean_score(metrics.get("pass_rates", {}).values())
    scorers = {
        "synthesis": score_synthesis,
        "docking_pose": score_docking_pose,
        "diversity": score_diversity,
    }
    return scorers[component](evidence, config)


def score_synthesis(evidence: dict[str, Any], config: dict[str, Any]) -> float:
    """Score synthesis tractability on a 0..100 scale."""
    metrics = evidence.get("metrics", {})
    synth_cfg = config.get("synthesis", {})
    parts = []

    if metrics.get("solve_rate") is not None:
        parts.append((0.40, clip01(float(metrics["solve_rate"]))))
    if metrics.get("median_sa") is not None:
        parts.append(
            (
                0.25,
                linear_score(
                    float(metrics["median_sa"]),
                    float(synth_cfg.get("sa_max", 4.5)),
                    float(synth_cfg.get("sa_min", 1.0)),
                    higher_is_better=False,
                ),
            )
        )
    if metrics.get("median_ra") is not None:
        parts.append(
            (
                0.20,
                linear_score(
                    float(metrics["median_ra"]),
                    float(synth_cfg.get("ra_min", 0.5)),
                    float(synth_cfg.get("ra_max", 1.0)),
                    higher_is_better=True,
                ),
            )
        )
    if metrics.get("median_syba") is not None:
        parts.append(
            (
                0.10,
                sigmoid_score(
                    float(metrics["median_syba"]),
                    float(synth_cfg.get("syba_midpoint", 0.0)),
                    float(synth_cfg.get("syba_scale", 50.0)),
                ),
            )
        )
    if metrics.get("median_search_time") is not None:
        target = float(synth_cfg.get("target_search_time_sec", 30.0))
        scale = float(synth_cfg.get("search_time_scale_sec", 20.0))
        parts.append(
            (
                0.05,
                linear_score(
                    float(metrics["median_search_time"]),
                    target + scale,
                    target,
                    higher_is_better=False,
                ),
            )
        )
    return _weighted_parts_score(parts)


def score_docking_pose(evidence: dict[str, Any], config: dict[str, Any]) -> float:
    """Score docking affinity and pose quality on a 0..100 scale."""
    metrics = evidence.get("metrics", {})
    docking_cfg = config.get("docking", {})
    parts = []
    if metrics.get("median_affinity") is not None:
        parts.append(
            (
                0.35,
                linear_score(
                    float(metrics["median_affinity"]),
                    float(docking_cfg.get("bad_affinity", -6.0)),
                    float(docking_cfg.get("good_affinity", -9.0)),
                    higher_is_better=False,
                ),
            )
        )
    if metrics.get("median_cnnscore") is not None:
        parts.append(
            (
                0.20,
                linear_score(
                    float(metrics["median_cnnscore"]),
                    float(docking_cfg.get("bad_cnnscore", 0.35)),
                    float(docking_cfg.get("good_cnnscore", 0.85)),
                    higher_is_better=True,
                ),
            )
        )
    if metrics.get("median_cnnaffinity") is not None:
        parts.append(
            (
                0.15,
                linear_score(
                    float(metrics["median_cnnaffinity"]),
                    float(docking_cfg.get("bad_cnnaffinity", 4.5)),
                    float(docking_cfg.get("good_cnnaffinity", 6.5)),
                    higher_is_better=True,
                ),
            )
        )
    if metrics.get("pose_pass_rate") is not None:
        parts.append((0.30, clip01(float(metrics["pose_pass_rate"]))))
    return _weighted_parts_score(parts)


def score_diversity(evidence: dict[str, Any], _config: dict[str, Any]) -> float:
    """Score diversity from MolEval metrics on a 0..100 scale."""
    metrics = evidence.get("metrics", {})
    parts = [
        (0.25, metrics.get("IntDiv1")),
        (0.25, metrics.get("IntDiv2")),
        (0.25, metrics.get("ScaffDiv")),
        (0.15, metrics.get("ScaffUniqueness")),
        (0.10, metrics.get("SEDiv")),
    ]
    return _weighted_parts_score(
        [(weight, clip01(float(value))) for weight, value in parts if value is not None]
    )


def assign_grade(score: float) -> str:
    """Convert a weighted score to a report grade."""
    if score >= 90:
        return "Excellent"
    if score >= 80:
        return "Strong"
    if score >= 65:
        return "Moderate"
    return "Weak"


def assign_confidence(
    final_count: int,
    available_components: int,
    config: dict[str, Any],
) -> str:
    """Assign confidence from final count and component availability."""
    confidence_cfg = config.get("confidence", {})
    high_min = int(confidence_cfg.get("min_final_molecules_high", 100))
    medium_min = int(confidence_cfg.get("min_final_molecules_medium", 30))
    if final_count >= high_min and available_components >= 5:
        return "High"
    if final_count >= medium_min and available_components >= 4:
        return "Medium"
    return "Low"


def generate_bottlenecks(components: dict[str, dict[str, Any]]) -> list[str]:
    """Return up to three human-readable bottlenecks from weak components."""
    available = [
        (name, component["score"])
        for name, component in components.items()
        if component.get("available") and component.get("score") is not None
    ]
    available.sort(key=lambda item: item[1])
    return [_bottleneck_text(name, score) for name, score in available if score < 85][
        :3
    ]


def _resolve_models(base_path: Path, available_models: list[str]) -> list[str]:
    if available_models:
        return sorted(str(model) for model in available_models)

    models = set()
    for rel_path in [
        "input/sampled_molecules.csv",
        "output/final_molecules.csv",
        "final_molecules.csv",
    ]:
        df = _read_csv(base_path / rel_path)
        if df is not None and "model_name" in df.columns:
            models.update(str(model) for model in df["model_name"].dropna().unique())
    if models:
        return sorted(models)

    final_path = first_existing(
        base_path,
        ["output/final_molecules.csv", "final_molecules.csv"],
    )
    return ["__all__"] if final_path else []


def _normalized_config_weights(
    config: dict[str, Any],
    key: str = "weights",
) -> dict[str, float]:
    raw_weights = config.get(key, {})
    weights = {
        component: float(raw_weights.get(component, 0.0))
        for component in COMPONENT_ORDER
    }
    total = sum(value for value in weights.values() if value > 0)
    if total <= 0:
        default = 1.0 / len(COMPONENT_ORDER)
        return {component: default for component in COMPONENT_ORDER}
    return {component: max(0.0, value) / total for component, value in weights.items()}


def _weighted_overall(
    components: dict[str, dict[str, Any]],
    weights: dict[str, float],
) -> float | None:
    available = [
        (weights.get(name, 0.0), component["score"])
        for name, component in components.items()
        if component.get("available") and component.get("score") is not None
    ]
    weight_sum = sum(weight for weight, _score in available)
    if weight_sum <= 0:
        return None
    return sum(weight * score for weight, score in available) / weight_sum


def _weighted_parts_score(parts: list[tuple[float, float]]) -> float:
    weight_sum = sum(weight for weight, _value in parts)
    if weight_sum <= 0:
        return 0.0
    return 100 * clip01(
        sum(weight * clip01(value) for weight, value in parts) / weight_sum
    )


def _applicable_hard_caps(
    evidence: dict[str, Any],
    config: dict[str, Any],
) -> list[dict[str, Any]]:
    hard_caps_cfg = config.get("hard_caps", {})
    if not hard_caps_cfg:
        return []

    caps = []

    structural_rate = (
        evidence.get("structural", {}).get("metrics", {}).get("stage_pass_rate")
    )
    if structural_rate is not None and float(structural_rate) < float(
        hard_caps_cfg.get("structural_stage_pass_rate_below", 0.20)
    ):
        caps.append(
            {
                "metric": "structural_stage_pass_rate",
                "value": float(structural_rate),
                "cap": float(hard_caps_cfg.get("structural_stage_pass_rate_cap", 60.0)),
                "warning": "Structural stage pass rate triggers score cap",
            }
        )

    descriptor_rate = (
        evidence.get("physchem", {}).get("metrics", {}).get("all_pass_rate")
    )
    if descriptor_rate is not None and float(descriptor_rate) < float(
        hard_caps_cfg.get("descriptor_all_pass_rate_below", 0.50)
    ):
        caps.append(
            {
                "metric": "descriptor_all_pass_rate",
                "value": float(descriptor_rate),
                "cap": float(hard_caps_cfg.get("descriptor_all_pass_rate_cap", 70.0)),
                "warning": "Descriptor all-pass rate triggers score cap",
            }
        )

    retention_rate = (
        evidence.get("yield", {}).get("metrics", {}).get("final_retention_rate")
    )
    if retention_rate is not None and float(retention_rate) < float(
        hard_caps_cfg.get("final_retention_rate_below", 0.05)
    ):
        caps.append(
            {
                "metric": "final_retention_rate",
                "value": float(retention_rate),
                "cap": float(hard_caps_cfg.get("final_retention_rate_cap", 70.0)),
                "warning": "Final retention rate triggers score cap",
            }
        )

    return caps


def _mean_score(values: Sequence[float]) -> float:
    clean = [clip01(float(value)) for value in values if value is not None]
    if not clean:
        return 0.0
    return 100 * statistics.mean(clean)


def _read_first_csv(
    base_path: Path,
    candidates: list[str],
    warnings: list[str],
) -> pd.DataFrame | None:
    saw_existing = False
    for candidate in candidates:
        path = base_path / candidate
        if not path.exists():
            continue
        saw_existing = True
        df = _read_csv(path)
        if df is None:
            warnings.append(f"Could not read CSV: {path}")
            continue
        if df.empty:
            warnings.append(f"Empty CSV: {path}")
            continue
        return df
    if not saw_existing:
        warnings.append(f"Missing CSV: {candidates[0]}")
    return None


def _read_first_csv_with_any_columns(
    base_path: Path,
    candidates: list[str],
    columns: list[str],
    warnings: list[str],
) -> pd.DataFrame | None:
    seen_existing = False
    lower_columns = {column.lower() for column in columns}
    for candidate in candidates:
        path = base_path / candidate
        if not path.exists():
            continue
        seen_existing = True
        df = _read_csv(path)
        if df is None:
            warnings.append(f"Could not read CSV: {path}")
            continue
        if {column.lower() for column in df.columns} & lower_columns:
            return df
    if not seen_existing:
        warnings.append(f"Missing CSV: {candidates[0]}")
    return None


def _read_csv(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    try:
        return pd.read_csv(path)
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
        return None


def _filter_model(df: pd.DataFrame | None, model_name: str) -> pd.DataFrame | None:
    if df is None:
        return None
    if model_name != "__all__" and "model_name" in df.columns:
        return df[df["model_name"] == model_name]
    return df


def _filtered_count(df: pd.DataFrame | None, model_name: str) -> int | None:
    df = _filter_model(df, model_name)
    if df is None:
        return None
    return int(len(df))


def _pass_columns(df: pd.DataFrame, suffix: str) -> list[str]:
    return [col for col in df.columns if col.endswith(suffix)]


def _mean_bool(series: pd.Series) -> float:
    values = series.dropna().map(_truthy)
    if values.empty:
        return 0.0
    return float(values.mean())


def _truthy(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if pd.isna(value):
        return False
    if isinstance(value, (int, float)):
        return bool(value)
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def _first_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for col in candidates:
        if col in df.columns:
            return col
    lower_map = {col.lower(): col for col in df.columns}
    for col in candidates:
        if col.lower() in lower_map:
            return lower_map[col.lower()]
    return None


def _evidence_text(component: str, metrics: dict[str, Any]) -> str:
    if component == "yield":
        retention = metrics.get("final_retention_rate")
        retention_text = (
            f", retention {100 * retention:.1f}%" if retention is not None else ""
        )
        return (
            f"{metrics.get('final_count', 0)} final / "
            f"{metrics.get('initial_count', 0)} initial"
            f"{retention_text}"
        )
    if component == "synthesis":
        parts = []
        if metrics.get("solve_rate") is not None:
            parts.append(f"solve rate {100 * metrics['solve_rate']:.1f}%")
        if metrics.get("median_sa") is not None:
            parts.append(f"median SA {metrics['median_sa']:.2f}")
        return ", ".join(parts)
    if component == "docking_pose":
        parts = []
        if metrics.get("median_affinity") is not None:
            parts.append(f"median affinity {metrics['median_affinity']:.2f}")
        if metrics.get("pose_pass_rate") is not None:
            parts.append(f"pose pass {100 * metrics['pose_pass_rate']:.1f}%")
        return ", ".join(parts)
    if component in {"physchem", "structural"}:
        if component == "physchem" and metrics.get("all_pass_rate") is not None:
            text = f"descriptor all-pass {100 * metrics['all_pass_rate']:.1f}%"
            if metrics.get("mean_flag_pass_rate") is not None:
                text += f", mean flag pass {100 * metrics['mean_flag_pass_rate']:.1f}%"
            return text
        if component == "structural" and metrics.get("stage_pass_rate") is not None:
            text = f"stage pass {100 * metrics['stage_pass_rate']:.1f}%"
            if (
                metrics.get("worst_filter")
                and metrics.get("worst_filter_pass_rate") is not None
            ):
                text += (
                    f", worst {metrics['worst_filter']} "
                    f"{100 * metrics['worst_filter_pass_rate']:.1f}%"
                )
            return text
        pass_rates = metrics.get("pass_rates", {})
        if pass_rates:
            return f"mean pass rate {100 * statistics.mean(pass_rates.values()):.1f}%"
    if component == "diversity":
        return ", ".join(
            f"{metric} {metrics[metric]:.3f}"
            for metric in DIVERSITY_METRICS
            if metric in metrics
        )
    return ""


def _bottleneck_text(component: str, score: float) -> str:
    if component == "synthesis":
        return "Synthesis solve rate below ideal"
    if component == "docking_pose":
        return "Docking/CNN confidence is moderate"
    if component == "diversity":
        return "Sphere-exclusion diversity is moderate"
    if component == "yield":
        return "Final retention is below target"
    if component == "physchem":
        return "Descriptor all-pass rate is below ideal"
    if component == "structural":
        return "Structural gate survival is below ideal"
    return f"{component} score is {score:.1f}"
