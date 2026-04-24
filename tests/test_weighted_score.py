"""Tests for Generator Reality Assessment scoring."""

from __future__ import annotations

import math

import pandas as pd
import pytest

from hedgehog.reporting.weighted_score import (
    assign_confidence,
    assign_grade,
    clip01,
    collect_physchem_evidence,
    compute_weighted_scores,
    linear_score,
    score_diversity,
    score_docking_pose,
    score_model_evidence,
    score_physchem,
    score_structural,
    score_synthesis,
    score_yield,
    sigmoid_score,
)

BASE_CONFIG = {
    "run": True,
    "version": "v1",
    "mode": "generator_reality",
    "target_final_count": 100,
    "target_final_retention": 0.10,
    "weights": {
        "yield": 0.30,
        "physchem": 0.15,
        "structural": 0.25,
        "synthesis": 0.10,
        "docking_pose": 0.15,
        "diversity": 0.05,
    },
    "candidate_pool_weights": {
        "yield": 0.10,
        "physchem": 0.15,
        "structural": 0.15,
        "synthesis": 0.20,
        "docking_pose": 0.30,
        "diversity": 0.10,
    },
    "yield": {
        "mode": "retention",
        "target_final_retention": 0.10,
        "count_weight": 0.70,
        "retention_weight": 0.30,
    },
    "structural": {"stage_pass_weight": 0.80, "worst_filter_weight": 0.20},
    "docking": {
        "bad_affinity": -6.0,
        "good_affinity": -9.0,
        "bad_cnnscore": 0.35,
        "good_cnnscore": 0.85,
        "bad_cnnaffinity": 4.5,
        "good_cnnaffinity": 6.5,
    },
    "synthesis": {
        "sa_min": 1.0,
        "sa_max": 4.5,
        "ra_min": 0.5,
        "ra_max": 1.0,
        "syba_midpoint": 0.0,
        "syba_scale": 50.0,
        "target_search_time_sec": 30.0,
        "search_time_scale_sec": 20.0,
    },
    "confidence": {
        "min_final_molecules_high": 100,
        "min_final_molecules_medium": 30,
    },
    "hard_caps": {
        "structural_stage_pass_rate_below": 0.20,
        "structural_stage_pass_rate_cap": 60.0,
        "descriptor_all_pass_rate_below": 0.50,
        "descriptor_all_pass_rate_cap": 70.0,
        "final_retention_rate_below": 0.05,
        "final_retention_rate_cap": 70.0,
    },
}


def test_clip01_clamps_values():
    assert clip01(-1) == 0
    assert clip01(0.4) == 0.4
    assert clip01(2) == 1


def test_linear_score_handles_both_directions():
    assert linear_score(0.85, 0.35, 0.85, higher_is_better=True) == 1
    assert linear_score(0.35, 0.35, 0.85, higher_is_better=True) == 0
    assert linear_score(-9.0, -6.0, -9.0, higher_is_better=False) == 1
    assert linear_score(-6.0, -6.0, -9.0, higher_is_better=False) == 0


def test_sigmoid_score_midpoint_is_half():
    assert sigmoid_score(0.0, midpoint=0.0, scale=50.0) == pytest.approx(0.5)
    assert sigmoid_score(50.0, midpoint=0.0, scale=50.0) > 0.5


def test_yield_score_uses_final_retention_target():
    evidence = {
        "metrics": {
            "final_count": 50,
            "initial_count": 1000,
            "final_retention_rate": 0.05,
        }
    }
    assert score_yield(evidence, BASE_CONFIG) == pytest.approx(50.0)


def test_yield_absolute_score_keeps_candidate_pool_formula():
    evidence = {"metrics": {"final_count": 100, "initial_count": 1000}}
    config = {**BASE_CONFIG, "yield": {**BASE_CONFIG["yield"], "mode": "absolute"}}
    expected = 100 * (
        0.70 * (1 - math.exp(-1)) + 0.30 * (math.log1p(100) / math.log1p(1000))
    )
    assert score_yield(evidence, config) == pytest.approx(expected)


def test_synthesis_score_is_bounded():
    evidence = {
        "metrics": {
            "solve_rate": 0.8,
            "median_sa": 2.5,
            "median_ra": 0.8,
            "median_syba": 20.0,
            "median_search_time": 25.0,
        }
    }
    score = score_synthesis(evidence, BASE_CONFIG)
    assert 0 <= score <= 100
    assert score > 50


def test_docking_score_is_bounded():
    evidence = {
        "metrics": {
            "median_affinity": -8.5,
            "median_cnnscore": 0.75,
            "median_cnnaffinity": 6.0,
            "pose_pass_rate": 0.9,
        }
    }
    score = score_docking_pose(evidence, BASE_CONFIG)
    assert 0 <= score <= 100
    assert score > 70


def test_diversity_score_uses_available_metrics():
    evidence = {
        "metrics": {
            "IntDiv1": 0.9,
            "IntDiv2": 0.8,
            "ScaffDiv": 0.7,
            "ScaffUniqueness": 0.6,
            "SEDiv": 0.5,
        }
    }
    assert score_diversity(evidence, BASE_CONFIG) == pytest.approx(74.0)


def test_structural_score_uses_stage_pass_and_worst_filter():
    evidence = {
        "metrics": {
            "stage_pass_rate": 0.20,
            "worst_filter_pass_rate": 0.10,
            "pass_rates": {"pass_common_alerts": 0.10, "pass_lilly": 0.95},
        }
    }
    assert score_structural(evidence, BASE_CONFIG) == pytest.approx(18.0)


def test_overall_renormalizes_missing_components():
    evidence = {
        "yield": {
            "available": True,
            "metrics": {"final_count": 100, "initial_count": 100},
            "warnings": [],
        },
        "physchem": {
            "available": True,
            "metrics": {"pass_rates": {"qed_pass": 1.0}},
            "warnings": [],
        },
        "structural": {"available": False, "metrics": {}, "warnings": []},
        "synthesis": {"available": False, "metrics": {}, "warnings": []},
        "docking_pose": {"available": False, "metrics": {}, "warnings": []},
        "diversity": {"available": False, "metrics": {}, "warnings": []},
    }
    score = score_model_evidence("ModelA", evidence, BASE_CONFIG)
    assert score["overall"] is not None
    assert score["components"]["structural"]["score"] is None
    assert score["components"]["structural"]["available"] is False
    assert score["candidate_pool_quality"]["overall"] is not None
    assert score["overall"] > 0


def test_hard_cap_limits_generator_score():
    evidence = {
        "yield": {
            "available": True,
            "metrics": {
                "final_count": 300,
                "initial_count": 1000,
                "final_retention_rate": 0.30,
            },
            "warnings": [],
        },
        "physchem": {
            "available": True,
            "metrics": {"all_pass_rate": 0.90, "pass_rates": {"qed_pass": 0.90}},
            "warnings": [],
        },
        "structural": {
            "available": True,
            "metrics": {
                "stage_pass_rate": 0.10,
                "worst_filter_pass_rate": 0.10,
                "pass_rates": {"pass_common_alerts": 0.10},
            },
            "warnings": [],
        },
        "synthesis": {
            "available": True,
            "metrics": {"solve_rate": 1.0, "median_sa": 1.0, "median_ra": 1.0},
            "warnings": [],
        },
        "docking_pose": {
            "available": True,
            "metrics": {
                "median_affinity": -9.0,
                "median_cnnscore": 0.85,
                "median_cnnaffinity": 6.5,
                "pose_pass_rate": 1.0,
            },
            "warnings": [],
        },
        "diversity": {
            "available": True,
            "metrics": {"IntDiv1": 1.0, "IntDiv2": 1.0},
            "warnings": [],
        },
    }
    score = score_model_evidence("ModelA", evidence, BASE_CONFIG)
    assert score["uncapped_overall"] > 60.0
    assert score["overall"] == 60.0
    assert score["hard_caps_applied"][0]["metric"] == "structural_stage_pass_rate"


def test_grade_assignment():
    assert assign_grade(92) == "Excellent"
    assert assign_grade(85) == "Strong"
    assert assign_grade(70) == "Moderate"
    assert assign_grade(40) == "Weak"


def test_confidence_assignment():
    assert assign_confidence(100, 5, BASE_CONFIG) == "High"
    assert assign_confidence(30, 4, BASE_CONFIG) == "Medium"
    assert assign_confidence(29, 6, BASE_CONFIG) == "Low"


def test_compute_weighted_scores_respects_run_false(tmp_path):
    result = compute_weighted_scores(
        tmp_path,
        ["ModelA"],
        {**BASE_CONFIG, "run": False},
    )
    assert result == {}


def test_physchem_prefers_initial_descriptor_flags(tmp_path):
    initial_dir = tmp_path / "stages" / "01_descriptors_initial" / "filtered"
    final_dir = tmp_path / "stages" / "07_descriptors_final" / "filtered"
    initial_dir.mkdir(parents=True)
    final_dir.mkdir(parents=True)

    pd.DataFrame(
        {
            "model_name": ["ModelA", "ModelA"],
            "qed_pass": [True, False],
            "logp_pass": [True, True],
        }
    ).to_csv(initial_dir / "pass_flags.csv", index=False)
    pd.DataFrame(
        {
            "model_name": ["ModelA", "ModelA"],
            "qed_pass": [True, True],
            "logp_pass": [True, True],
        }
    ).to_csv(final_dir / "pass_flags.csv", index=False)

    evidence = collect_physchem_evidence(tmp_path, "ModelA")

    assert evidence["available"] is True
    assert score_physchem(evidence, BASE_CONFIG) == pytest.approx(50.0)


def test_physchem_falls_back_when_initial_flags_are_empty(tmp_path):
    initial_dir = tmp_path / "stages" / "01_descriptors_initial" / "filtered"
    final_dir = tmp_path / "stages" / "07_descriptors_final" / "filtered"
    initial_dir.mkdir(parents=True)
    final_dir.mkdir(parents=True)

    pd.DataFrame(columns=["model_name", "qed_pass", "logp_pass"]).to_csv(
        initial_dir / "pass_flags.csv", index=False
    )
    pd.DataFrame(
        {
            "model_name": ["ModelA", "ModelA"],
            "qed_pass": [True, True],
            "logp_pass": [True, True],
        }
    ).to_csv(final_dir / "pass_flags.csv", index=False)

    evidence = collect_physchem_evidence(tmp_path, "ModelA")

    assert evidence["available"] is True
    assert any("Empty CSV" in warning for warning in evidence["warnings"])
    assert score_physchem(evidence, BASE_CONFIG) == pytest.approx(100.0)


def test_physchem_filters_initial_flags_by_model(tmp_path):
    initial_dir = tmp_path / "stages" / "01_descriptors_initial" / "filtered"
    initial_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "model_name": ["ModelA", "ModelA", "ModelB", "ModelB"],
            "qed_pass": [True, False, True, True],
            "logp_pass": [False, False, True, True],
        }
    ).to_csv(initial_dir / "pass_flags.csv", index=False)

    model_a = collect_physchem_evidence(tmp_path, "ModelA")
    model_b = collect_physchem_evidence(tmp_path, "ModelB")

    assert score_physchem(model_a, BASE_CONFIG) == pytest.approx(0.0)
    assert score_physchem(model_b, BASE_CONFIG) == pytest.approx(100.0)
