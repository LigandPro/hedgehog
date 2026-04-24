"""Integration tests for Generator Reality Assessment report data."""

from __future__ import annotations

import pandas as pd
import pytest
import yaml

from hedgehog.reporting.report_generator import ReportGenerator

MODEL_NAME = "LiTr_moses_full_75_aug"


def _write_weighted_config(path):
    config = {
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
    path.write_text(yaml.safe_dump(config))


def _build_weighted_run(tmp_path, with_config=True):
    (tmp_path / "input").mkdir()
    (tmp_path / "output").mkdir()
    (tmp_path / "configs").mkdir()
    (tmp_path / "stages" / "03_structural_filters_post").mkdir(parents=True)
    (tmp_path / "stages" / "04_synthesis").mkdir(parents=True)
    (tmp_path / "stages" / "06_docking_filters").mkdir(parents=True)
    (tmp_path / "stages" / "01_descriptors_initial" / "filtered").mkdir(parents=True)
    (tmp_path / "stages" / "01_descriptors_initial" / "metrics").mkdir(parents=True)
    (tmp_path / "stages" / "07_descriptors_final" / "filtered").mkdir(parents=True)
    (tmp_path / "stages" / "07_descriptors_final" / "metrics").mkdir(parents=True)

    initial = 120
    final = 105
    pd.DataFrame(
        {
            "smiles": ["CCO"] * initial,
            "model_name": [MODEL_NAME] * initial,
            "mol_idx": list(range(initial)),
        }
    ).to_csv(tmp_path / "input" / "sampled_molecules.csv", index=False)

    pd.DataFrame(
        {
            "smiles": ["CCO"] * final,
            "model_name": [MODEL_NAME] * final,
            "mol_idx": list(range(final)),
            "gnina_affinity": [-8.5] * final,
            "gnina_cnnscore": [0.75] * final,
            "gnina_cnnaffinity": [6.0] * final,
            "pass_search_box": [True] * final,
            "pass_pose_quality": [True] * final,
            "pass_interactions": [True] * final,
            "pass_conformer_deviation": [True] * final,
            "pass": [True] * final,
        }
    ).to_csv(tmp_path / "output" / "final_molecules.csv", index=False)

    initial_pass_flags = pd.DataFrame(
        {
            "model_name": [MODEL_NAME] * initial,
            "qed_pass": [True] * initial,
            "logp_pass": [True] * 60 + [False] * 60,
            "tpsa_pass": [True] * initial,
        }
    )
    initial_pass_flags.to_csv(
        tmp_path / "stages" / "01_descriptors_initial" / "filtered" / "pass_flags.csv",
        index=False,
    )
    pd.DataFrame(
        {
            "model_name": [MODEL_NAME] * initial,
            "QED": [0.7] * initial,
            "MolWt": [350.0] * initial,
            "LogP": [2.5] * initial,
            "TPSA": [80.0] * initial,
            "Fsp3": [0.4] * initial,
        }
    ).to_csv(
        tmp_path
        / "stages"
        / "01_descriptors_initial"
        / "metrics"
        / "descriptors_all.csv",
        index=False,
    )

    final_pass_flags = pd.DataFrame(
        {
            "model_name": [MODEL_NAME] * final,
            "qed_pass": [True] * final,
            "logp_pass": [True] * final,
            "tpsa_pass": [True] * final,
        }
    )
    final_pass_flags.to_csv(
        tmp_path / "stages" / "07_descriptors_final" / "filtered" / "pass_flags.csv",
        index=False,
    )
    pd.DataFrame(
        {
            "model_name": [MODEL_NAME] * final,
            "QED": [0.7] * final,
            "MolWt": [350.0] * final,
            "LogP": [2.5] * final,
            "TPSA": [80.0] * final,
            "Fsp3": [0.4] * final,
        }
    ).to_csv(
        tmp_path
        / "stages"
        / "07_descriptors_final"
        / "metrics"
        / "descriptors_all.csv",
        index=False,
    )

    pd.DataFrame(
        {
            "model_name": [MODEL_NAME] * final,
        }
    ).to_csv(
        tmp_path / "stages" / "03_structural_filters_post" / "filtered_molecules.csv",
        index=False,
    )
    pd.DataFrame(
        {
            "model_name": [MODEL_NAME] * 15,
            "pass_common_alerts": [False] * 15,
            "pass_molgraph_stats": [True] * 15,
            "pass_molcomplexity": [True] * 15,
            "pass_NIBR": [True] * 15,
            "pass_bredt": [True] * 15,
            "pass_lilly": [True] * 15,
            "pass_protecting_groups": [True] * 15,
            "pass_ring_infraction": [True] * 15,
            "pass_stereo_center": [True] * 15,
            "pass_halogenicity": [True] * 15,
        }
    ).to_csv(
        tmp_path / "stages" / "03_structural_filters_post" / "failed_molecules.csv",
        index=False,
    )

    pd.DataFrame(
        {
            "model_name": [MODEL_NAME] * final,
            "sa_score": [2.5] * final,
            "syba_score": [20.0] * final,
            "ra_score": [0.8] * final,
            "solved": [True] * 90 + [False] * 15,
            "search_time": [25.0] * final,
        }
    ).to_csv(
        tmp_path / "stages" / "04_synthesis" / "filtered_molecules.csv",
        index=False,
    )

    config = {}
    if with_config:
        config_path = tmp_path / "configs" / "config_weighted_score.yml"
        _write_weighted_config(config_path)
        config["config_weighted_score"] = str(config_path)

    return ReportGenerator(
        base_path=tmp_path,
        stages=[],
        config=config,
        initial_count=initial,
        final_count=final,
    )


def test_collect_data_includes_weighted_scores(tmp_path):
    gen = _build_weighted_run(tmp_path)

    data = gen._collect_data()

    assert "weighted_scores" in data
    assert MODEL_NAME in data["weighted_scores"]["models"]
    score = data["weighted_scores"]["models"][MODEL_NAME]
    assert 0 <= score["overall"] <= 100
    assert score["confidence"] == "High"
    assert set(score.keys()) >= {
        "overall",
        "grade",
        "confidence",
        "candidate_pool_quality",
        "hard_caps_applied",
        "components",
        "bottlenecks",
        "warnings",
    }
    assert score["components"]["physchem"]["score"] == pytest.approx(50.0)
    assert score["components"]["physchem"]["evidence"]["all_pass_rate"] == 0.5
    assert score["components"]["physchem"]["evidence"][
        "mean_flag_pass_rate"
    ] == pytest.approx(100 / 120)
    assert score["components"]["physchem"]["evidence"]["worst_flag"] == "logp_pass"
    assert score["components"]["structural"]["score"] == pytest.approx(87.5)
    assert score["components"]["structural"]["evidence"]["stage_pass_rate"] == (
        105 / 120
    )
    assert 0 <= score["candidate_pool_quality"]["overall"] <= 100
    assert score["candidate_pool_quality"]["grade"] in {
        "Excellent",
        "Strong",
        "Moderate",
        "Weak",
    }


def test_model_stats_use_output_final_molecules_fallback(tmp_path):
    gen = _build_weighted_run(tmp_path)

    assert not (tmp_path / "final_molecules.csv").exists()
    model_stats = gen._get_model_stats()

    assert model_stats
    assert model_stats[0]["model_name"] == MODEL_NAME
    assert model_stats[0]["final"] == 105


def test_collect_data_has_model_stats_with_output_final_molecules(tmp_path):
    gen = _build_weighted_run(tmp_path)

    data = gen._collect_data()

    assert data["models"] != []


def test_generate_report_renders_generator_reality_section(tmp_path):
    gen = _build_weighted_run(tmp_path)

    report_path = gen.generate()

    html = report_path.read_text()
    assert "Generator Reality Assessment" in html
    assert "Final Candidate Pool Quality" in html
    assert MODEL_NAME in html
    assert "High" in html


def test_missing_weighted_score_config_disables_score_without_failure(tmp_path):
    gen = _build_weighted_run(tmp_path, with_config=False)

    data = gen._collect_data()

    assert data["weighted_scores"] == {}
