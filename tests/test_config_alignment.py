from pathlib import Path

import pandas as pd
import pytest
import yaml
from rdkit import Chem

from hedgehog._constants import KEY_ALIGNMENT_SKIP_FINAL_DESCRIPTORS
from hedgehog.config_alignment import (
    GeneratedConfigShapeError,
    _align_descriptor_config,
    _align_docking_config,
    _align_structural_filter_config,
    _dump_aligned_stage_yaml,
    _dump_generated_yaml,
    _read_structural_rule_masks,
    _yaml_number,
    create_aligned_stage_config,
    create_probe_config,
    descriptor_bounds_mode_from_master,
    set_probe_molprep_allowed_atoms,
    validate_descriptor_bounds_mode,
    validate_target_coverage_percent,
)
from hedgehog.main import _alignment_target_coverage_from_config


def _write_yaml(path: Path, value: dict) -> None:
    path.write_text(yaml.safe_dump(value, sort_keys=False), encoding="utf-8")


def _base_master(tmp_path: Path) -> dict:
    molprep = tmp_path / "molprep.yml"
    descriptors = tmp_path / "descriptors.yml"
    synthesis = tmp_path / "synthesis.yml"
    docking = tmp_path / "docking.yml"
    docking_filters = tmp_path / "docking_filters.yml"
    struct_filters = tmp_path / "struct_filters.yml"

    _write_yaml(
        molprep,
        {
            "run": True,
            "set_allowed_atoms_from_targets": True,
            "columns": {"smiles": "smiles"},
            "filters": {"allowed_atoms": ["C", "N", "O"]},
        },
    )
    _write_yaml(
        descriptors,
        {
            "run": True,
            "filter_data": True,
            "structural_constraints": {
                "enabled": True,
                "type_limits": {"Hal": 3},
                "element_limits": {"N": 6},
                "max_n_or_o_atoms": 10,
            },
            "borders": {"metric_min": -10, "metric_max": 200},
        },
    )
    _write_yaml(
        synthesis,
        {
            "run": True,
            "enabled_scores": None,
            "run_retrosynthesis": True,
            "filter_solved_only": True,
            "sa_score_min": 1,
            "sa_score_max": 10,
            "score_filters": {"sync_score": {"min": 0.5, "max": 1}},
        },
    )
    _write_yaml(
        docking,
        {
            "run": True,
            "calculate_score_thresholds_from_targets": True,
            "score_thresholds": {},
            "tools": "all",
        },
    )
    _write_yaml(
        docking_filters,
        {
            "run": True,
            "search_box": {
                "max_outside_fraction": 0,
                "short_circuit": True,
            },
            "pose_quality": {
                "clash_cutoff": 0.75,
                "volume_clash_cutoff": 0.075,
                "max_distance": 5.0,
            },
            "interactions": {"min_hbonds": 1, "required_residues": ["ASP12"]},
            "shepherd_score": {"min_shape_score": 0.5},
            "conformer_deviation": {
                "max_rmsd_to_conformer": 3,
                "early_stop_on_pass": True,
            },
            "aggregation": {"mode": "all", "save_metrics": True},
        },
    )
    alerts = tmp_path / "alerts.csv"
    pd.DataFrame({"rule_set_name": ["Rule A", "Rule B"], "smarts": ["C", "N"]}).to_csv(
        alerts, index=False
    )
    _write_yaml(
        struct_filters,
        {
            "run": True,
            "filter_data": True,
            "write_per_filter_outputs": False,
            "generate_plots": True,
            "generate_failure_analysis": True,
            "calculate_common_alerts": False,
            "calculate_bredt": False,
            "calculate_lilly": False,
            "calculate_ring_infraction": False,
            "ring_infraction_hetcycle_min_size": 4,
            "calculate_stereo_center": False,
            "stereo_max_centers": 4,
            "stereo_max_undefined": 2,
            "calculate_halogenicity": False,
            "halogenicity_thresh_F": 6,
            "halogenicity_thresh_Br": 3,
            "halogenicity_thresh_Cl": 3,
            "calculate_symmetry": False,
            "symmetry_threshold": 0.8,
            "alerts_data_path": str(alerts),
            "include_rulesets": [],
            "exclude_smarts": ["ignored"],
        },
    )

    return {
        "generated_mols_path": str(tmp_path / "candidates.csv"),
        "target_mols_path": str(tmp_path / "targets.csv"),
        "alignment": {
            "enabled": True,
            "target_coverage_percent": 95,
        },
        "folder_to_save": str(tmp_path / "run"),
        "sample_size": 10,
        "config_mol_prep": str(molprep),
        "config_descriptors": str(descriptors),
        "config_structFilters": str(struct_filters),
        "config_synthesis": str(synthesis),
        "config_docking": str(docking),
        "config_docking_filters": str(docking_filters),
    }


@pytest.mark.parametrize(
    "value,outward,integer_values,expected",
    [
        (1.234, "min", False, 1.23),
        (1.231, "max", False, 1.24),
        (-7.123, "max", False, -7.12),
        (4.01, "max", True, 5),
        (112.01000000000002, "max", False, 112.02),
        (1.15, "min", False, 1.15),
    ],
)
def test_generated_threshold_numbers_are_normalized(
    value, outward, integer_values, expected
):
    assert (
        _yaml_number(value, outward=outward, integer_values=integer_values) == expected
    )


def test_generated_yaml_preserves_source_shape_comments_and_disabled_fields(tmp_path):
    source = tmp_path / "source.yml"
    source_text = (
        "# source comment\n"
        "run: true\n"
        "n_jobs: 96  # keep this comment\n"
        "borders:\n"
        "  metric_min: 0\n"
        "  metric_max: 10\n"
        "  # optional_min: 1\n"
        "  # optional_max: 2\n"
    )
    source.write_text(source_text, encoding="utf-8")
    target = tmp_path / "generated.yml"

    generated = _dump_generated_yaml(
        {
            "run": True,
            "n_jobs": 96,
            "borders": {
                "metric_min": 2,
                "metric_max": 8,
                "unexpected_min": -100,
            },
            "unexpected": "discarded",
        },
        target,
        source=source,
    )

    assert source.read_text(encoding="utf-8") == source_text
    assert list(generated) == ["run", "n_jobs", "borders"]
    assert list(generated["borders"]) == ["metric_min", "metric_max"]
    assert target.read_text(encoding="utf-8") == (
        "# source comment\n"
        "run: true\n"
        "n_jobs: 96  # keep this comment\n"
        "borders:\n"
        "  metric_min: 2\n"
        "  metric_max: 8\n"
        "  # optional_min: 1\n"
        "  # optional_max: 2\n"
    )


def test_generated_yaml_keeps_separator_after_replaced_block_list(tmp_path):
    source = tmp_path / "source.yml"
    source.write_text(
        "include_rulesets:\n"
        "  - Dundee\n"
        "  - BMS\n"
        "exclude_smarts:\n"
        "  - ignored\n",
        encoding="utf-8",
    )
    target = tmp_path / "generated.yml"

    generated = _dump_generated_yaml(
        {
            "include_rulesets": ["AlphaScreen-Hitters"],
            "exclude_smarts": ["kept"],
        },
        target,
        source=source,
    )

    assert yaml.safe_load(target.read_text(encoding="utf-8")) == generated
    assert target.read_text(encoding="utf-8") == (
        "include_rulesets:\n"
        "  - AlphaScreen-Hitters\n"
        "exclude_smarts:\n"
        "  - kept\n"
    )


def test_generated_yaml_preserves_key_indent_after_indentless_block_list(tmp_path):
    source = tmp_path / "source.yml"
    source.write_text(
        "filters:\n"
        "  allowed_atoms:\n"
        "  - C\n"
        "  - N\n"
        "  reject_radicals: true\n",
        encoding="utf-8",
    )
    target = tmp_path / "generated.yml"

    generated = _dump_generated_yaml(
        {
            "filters": {
                "allowed_atoms": ["C", "N", "Cl"],
                "reject_radicals": True,
            }
        },
        target,
        source=source,
    )

    assert yaml.safe_load(target.read_text(encoding="utf-8")) == generated
    assert target.read_text(encoding="utf-8") == (
        "filters:\n"
        "  allowed_atoms:\n"
        "  - C\n"
        "  - N\n"
        "  - Cl\n"
        "  reject_radicals: true\n"
    )


def test_aligned_stage_yaml_rejects_non_threshold_changes(tmp_path):
    source = tmp_path / "descriptors.yml"
    _write_yaml(
        source,
        {
            "run": True,
            "filter_data": False,
            "borders": {"metric_min": 0, "metric_max": 10},
        },
    )

    with pytest.raises(
        GeneratedConfigShapeError,
        match=r"non-threshold source fields: filter_data",
    ):
        _dump_aligned_stage_yaml(
            {
                "run": True,
                "filter_data": True,
                "borders": {"metric_min": 2, "metric_max": 8},
            },
            tmp_path / "generated.yml",
            source=source,
            config_key="config_descriptors",
        )


def test_aligned_stage_yaml_copies_everything_except_threshold_values(tmp_path):
    source = tmp_path / "descriptors.yml"
    source_text = (
        "# original\n"
        "run: true\n"
        "filter_data: false  # unchanged\n"
        "borders:\n"
        "  metric_min: 0\n"
        "  metric_max: 10\n"
    )
    source.write_text(source_text, encoding="utf-8")
    target = tmp_path / "generated.yml"

    _dump_aligned_stage_yaml(
        {
            "run": True,
            "filter_data": False,
            "borders": {"metric_min": 2, "metric_max": 8},
        },
        target,
        source=source,
        config_key="config_descriptors",
    )

    assert target.read_text(encoding="utf-8") == source_text.replace(
        "metric_min: 0", "metric_min: 2"
    ).replace("metric_max: 10", "metric_max: 8")


def test_aligned_docking_yaml_allows_existing_threshold_values_to_change(tmp_path):
    source = tmp_path / "docking.yml"
    _write_yaml(
        source,
        {
            "run": True,
            "score_thresholds": {
                "gnina": {
                    "score_property": "minimizedAffinity",
                    "max": -5,
                }
            },
            "tools": "gnina",
        },
    )
    target = tmp_path / "generated.yml"

    generated = _dump_aligned_stage_yaml(
        {
            "run": True,
            "score_thresholds": {
                "gnina": {
                    "score_property": "minimizedAffinity",
                    "max": -7,
                }
            },
            "tools": "gnina",
        },
        target,
        source=source,
        config_key="config_docking",
    )

    assert generated["score_thresholds"]["gnina"]["max"] == -7
    assert generated["run"] is True
    assert generated["tools"] == "gnina"


@pytest.mark.parametrize(
    ("scores", "percentile", "calibrated_maximum", "effective_maximum"),
    [
        ([-10.0, -9.0, -8.0], 50, -9, -6.5),
        ([-6.0, -5.0, -4.0], 100, -4, -4),
    ],
)
def test_docking_alignment_can_relax_but_not_tighten_configured_cutoff(
    scores, percentile, calibrated_maximum, effective_maximum
):
    config = {
        "calculate_score_thresholds_from_targets": True,
        "tools": "smina",
        "score_thresholds": {
            "smina": {
                "score_property": "minimizedAffinity",
                "max": -6.5,
            }
        },
    }

    audit = _align_docking_config(
        config,
        pd.DataFrame({"smina": scores}),
        percentile,
    )

    assert audit["calibrated_score_thresholds"]["smina"]["max"] == calibrated_maximum
    assert audit["score_thresholds"]["smina"]["max"] == effective_maximum
    assert config["score_thresholds"]["smina"]["max"] == effective_maximum


def test_descriptor_alignment_reads_ring_size_lists_without_replacing_them():
    config = {"borders": {"ring_size_min": 3, "ring_size_max": 12}}
    metrics = pd.DataFrame(
        {"ring_size": ["[3, 6]", "[4, 7]", "[5, 8]", "[6, 9]", "[]"]}
    )
    original = metrics["ring_size"].copy()

    thresholds = _align_descriptor_config(config, metrics, 80)

    assert metrics["ring_size"].equals(original)
    assert thresholds == {"ring_size_min": 3, "ring_size_max": 12}
    assert config["borders"] == {"ring_size_min": 3, "ring_size_max": 12}


def test_descriptor_alignment_expand_never_narrows_source_borders():
    config = {"borders": {"metric_min": 0, "metric_max": 10}}
    metrics = pd.DataFrame({"metric": [2, 4, 8]})

    thresholds = _align_descriptor_config(
        config, metrics, 100, bounds_mode="expand"
    )

    assert thresholds == {"metric_min": 0, "metric_max": 10}
    assert config["borders"] == thresholds


def test_descriptor_alignment_expand_expands_only_required_sides():
    config = {"borders": {"metric_min": 0, "metric_max": 10}}
    metrics = pd.DataFrame({"metric": [-5, 2, 8]})

    thresholds = _align_descriptor_config(
        config, metrics, 100, bounds_mode="expand"
    )

    assert thresholds == {"metric_min": -5, "metric_max": 10}
    assert config["borders"] == thresholds


def test_descriptor_alignment_target_keeps_existing_replacement_behavior():
    config = {"borders": {"metric_min": 0, "metric_max": 10}}
    metrics = pd.DataFrame({"metric": [2, 4, 8]})

    thresholds = _align_descriptor_config(
        config, metrics, 100, bounds_mode="target"
    )

    assert thresholds == {"metric_min": 2, "metric_max": 8}


@pytest.mark.parametrize("value", ["full", "", None, True, "expand_only", "target_only"])
def test_validate_descriptor_bounds_mode_rejects_unknown_values(value):
    with pytest.raises(ValueError, match="descriptor_bounds_mode"):
        validate_descriptor_bounds_mode(value)


def test_descriptor_bounds_mode_defaults_to_expand():
    assert descriptor_bounds_mode_from_master({}) == "expand"
    assert descriptor_bounds_mode_from_master({"alignment": {}}) == "expand"


@pytest.mark.parametrize("value", [0, -1, 100.1, float("inf"), True, "95"])
def test_validate_target_coverage_percent_rejects_invalid_values(value):
    with pytest.raises(ValueError):
        validate_target_coverage_percent(value)


def test_alignment_disabled_in_master_config_returns_none():
    config = {
        "alignment": {
            "enabled": False,
            "target_coverage_percent": 95,
        }
    }

    assert _alignment_target_coverage_from_config(config, None) is None


def test_alignment_enabled_in_master_config_uses_configured_target_coverage():
    config = {
        "alignment": {
            "enabled": True,
            "target_coverage_percent": 87.5,
        }
    }

    assert _alignment_target_coverage_from_config(config, None) == 87.5


def test_cli_alignment_target_coverage_overrides_disabled_master_config():
    config = {
        "alignment": {
            "enabled": False,
            "target_coverage_percent": 80,
        }
    }

    assert _alignment_target_coverage_from_config(config, 92) == 92


@pytest.mark.parametrize(
    "settings, message",
    [
        ({"enabled": True}, "target_coverage_percent is required"),
        ({"enabled": "yes", "target_coverage_percent": 95}, "must be true or false"),
        ({"enabled": True, "target_coverage_percent": 0}, "greater than 0"),
    ],
)
def test_invalid_master_alignment_config_is_rejected(settings, message):
    with pytest.raises(ValueError, match=message):
        _alignment_target_coverage_from_config({"alignment": settings}, None)


def test_invalid_descriptor_bounds_mode_in_master_config_is_rejected():
    config = {
        "alignment": {
            "enabled": True,
            "target_coverage_percent": 95,
            "descriptor_bounds_mode": "replace_if_needed",
        }
    }

    with pytest.raises(ValueError, match="descriptor_bounds_mode"):
        _alignment_target_coverage_from_config(config, None)


def test_legacy_retention_percentile_is_still_accepted():
    config = {
        "alignment": {
            "enabled": True,
            "retention_percentile": 75,
        }
    }

    assert _alignment_target_coverage_from_config(config, None) == 75


def test_conflicting_coverage_and_legacy_values_are_rejected():
    config = {
        "alignment": {
            "enabled": True,
            "target_coverage_percent": 80,
            "retention_percentile": 90,
        }
    }

    with pytest.raises(ValueError, match="conflicts with deprecated"):
        _alignment_target_coverage_from_config(config, None)


def test_probe_config_bypasses_synthesis_without_changing_its_policy(tmp_path):
    master = _base_master(tmp_path)
    master["_run_stage_selection_override"] = ["mol_prep", "descriptors"]
    targets = tmp_path / "targets.csv"
    targets.write_text("smiles\nCCO\n", encoding="utf-8")

    probe = create_probe_config(master, str(targets), tmp_path / "alignment")
    allowed_atoms = set_probe_molprep_allowed_atoms(
        probe,
        pd.DataFrame({"smiles": ["CCO", "C[SiH3]", "[Na+].[Cl-]", "[H]C"]}),
    )

    assert probe[KEY_ALIGNMENT_SKIP_FINAL_DESCRIPTORS] is True
    assert probe["_run_stage_selection_override"] == ["mol_prep", "descriptors"]

    original_molprep = yaml.safe_load(Path(master["config_mol_prep"]).read_text())
    probe_molprep = yaml.safe_load(Path(probe["config_mol_prep"]).read_text())
    original_descriptors = yaml.safe_load(
        Path(master["config_descriptors"]).read_text()
    )
    original_synthesis = yaml.safe_load(Path(master["config_synthesis"]).read_text())
    probe_descriptors = yaml.safe_load(Path(probe["config_descriptors"]).read_text())
    probe_synthesis = yaml.safe_load(Path(probe["config_synthesis"]).read_text())
    probe_struct = yaml.safe_load(Path(probe["config_structFilters"]).read_text())
    probe_docking = yaml.safe_load(Path(probe["config_docking"]).read_text())
    probe_filters = yaml.safe_load(Path(probe["config_docking_filters"]).read_text())

    assert allowed_atoms == ["H", "C", "O", "Na", "Si", "Cl"]
    assert original_molprep["filters"]["allowed_atoms"] == ["C", "N", "O"]
    assert probe_molprep["filters"]["allowed_atoms"] == allowed_atoms
    assert original_descriptors["borders"] == {
        "metric_min": -10,
        "metric_max": 200,
    }
    assert probe_descriptors["borders"] == {}
    assert probe_descriptors["structural_constraints"]["enabled"] is False
    assert probe_synthesis == original_synthesis | {"run": False}
    assert Path(probe["config_synthesis"]).parent.name == (
        "calibration_configs_unfiltered"
    )
    assert Path(probe["folder_to_save"]).name == "calibration_target_run"
    source_master = tmp_path / "alignment" / "source_configs" / "source_config.yml"
    assert source_master.is_file()
    source = yaml.safe_load(source_master.read_text())
    assert "_run_stage_selection_override" not in source
    source_synthesis = yaml.safe_load(Path(source["config_synthesis"]).read_text())
    assert source_synthesis == original_synthesis
    assert probe_struct["run"] is True
    assert probe_struct["filter_data"] is False
    assert probe_struct["write_per_filter_outputs"] is True
    assert probe_struct["generate_plots"] is False
    assert probe_struct["calculate_common_alerts"] is True
    assert probe_struct["calculate_bredt"] is True
    assert probe_struct["calculate_lilly"] is True
    assert probe_struct["calculate_symmetry"] is True
    assert probe_struct["include_rulesets"] == "all"
    assert probe_struct["exclude_smarts"] == []
    assert probe_docking["calculate_score_thresholds_from_targets"] is True
    assert probe_docking["score_thresholds"] == {}
    assert probe_filters["run"] is False
    assert probe["sample_size"] is None


def test_molprep_config_is_created_after_stage_with_all_target_atoms(tmp_path):
    master = _base_master(tmp_path)
    targets = tmp_path / "targets.csv"
    targets.write_text("smiles\nCCO\nC[SiH3]\n[Na+].[Cl-]\n", encoding="utf-8")
    target_run = tmp_path / "target_run"
    sampled_path = target_run / "input" / "sampled_molecules.csv"
    sampled_path.parent.mkdir(parents=True)
    pd.DataFrame({"smiles": ["CCO", "C[SiH3]", "[Na+].[Cl-]"]}).to_csv(
        sampled_path, index=False
    )

    result = create_aligned_stage_config(
        master,
        target_run,
        tmp_path / "alignment",
        str(targets),
        90,
        "mol_prep",
    )

    assert result is not None
    aligned, _master_path, thresholds_path = result
    aligned_molprep = yaml.safe_load(Path(aligned["config_mol_prep"]).read_text())
    assert aligned_molprep["filters"]["allowed_atoms"] == [
        "C",
        "O",
        "Na",
        "Si",
        "Cl",
    ]
    summary = yaml.safe_load(thresholds_path.read_text())
    assert summary["stages"]["mol_prep"] == {
        "thresholds": {"filters.allowed_atoms": ["C", "O", "Na", "Si", "Cl"]},
        "status": "ready",
        "note": "Allowed atoms are derived from all target molecules.",
    }
    assert not (
        tmp_path / "alignment" / "aligned_configs" / "config_descriptors.yml"
    ).exists()


def test_incremental_descriptor_alignment_propagates_expand_mode(tmp_path):
    master = _base_master(tmp_path)
    master["alignment"]["descriptor_bounds_mode"] = "expand"
    targets = tmp_path / "targets.csv"
    targets.write_text("smiles\nCCO\n", encoding="utf-8")
    metrics_dir = (
        tmp_path / "target_run" / "stages" / "02_descriptors_initial" / "metrics"
    )
    metrics_dir.mkdir(parents=True)
    pd.DataFrame({"metric": [1, 2, 100]}).to_csv(
        metrics_dir / "descriptors_all.csv", index=False
    )

    result = create_aligned_stage_config(
        master,
        tmp_path / "target_run",
        tmp_path / "alignment",
        str(targets),
        100,
        "descriptors",
    )

    assert result is not None
    aligned, _master_path, thresholds_path = result
    descriptor_config = yaml.safe_load(Path(aligned["config_descriptors"]).read_text())
    summary = yaml.safe_load(thresholds_path.read_text())
    assert descriptor_config["borders"] == {"metric_min": -10, "metric_max": 200}
    assert summary["descriptor_bounds_mode"] == "expand"
    assert summary["stages"]["descriptors"]["bounds_mode"] == "expand"
    assert summary["stages"]["descriptors"]["thresholds"] == {
        "metric_min": -10,
        "metric_max": 200,
    }


def test_aligned_configs_are_created_one_completed_stage_at_a_time(tmp_path):
    master = _base_master(tmp_path)
    targets = tmp_path / "targets.csv"
    targets.write_text("smiles\nCCO\n", encoding="utf-8")
    target_run = tmp_path / "target_run"
    descriptors_dir = target_run / "stages" / "02_descriptors_initial" / "metrics"
    descriptors_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "metric": [1, 2, 3],
            "n_N_atoms": [1, 2, 3],
            "n_NO_atoms": [2, 3, 4],
        }
    ).to_csv(descriptors_dir / "descriptors_all.csv", index=False)

    descriptor_result = create_aligned_stage_config(
        master,
        target_run,
        tmp_path / "alignment",
        str(targets),
        100,
        "descriptors",
    )

    assert descriptor_result is not None
    aligned_dir = tmp_path / "alignment" / "aligned_configs"
    assert (aligned_dir / "config_descriptors.yml").exists()
    assert not (aligned_dir / "config_synthesis.yml").exists()
    assert not (aligned_dir / "config_docking.yml").exists()
    first_summary = yaml.safe_load(
        (aligned_dir / "alignment_thresholds.yml").read_text()
    )
    assert first_summary["stages"]["descriptors"]["status"] == "ready"
    assert "synthesis" not in first_summary["stages"]


def test_completed_stage_can_be_recalibrated_in_place(tmp_path):
    master = _base_master(tmp_path)
    targets = tmp_path / "targets.csv"
    targets.write_text("smiles\nCCO\n", encoding="utf-8")
    target_run = tmp_path / "target_run"
    metrics_dir = target_run / "stages" / "02_descriptors_initial" / "metrics"
    metrics_dir.mkdir(parents=True)
    pd.DataFrame({"metric": [1.0, 2.0, 100.0]}).to_csv(
        metrics_dir / "descriptors_all.csv", index=False
    )
    first = create_aligned_stage_config(
        master, target_run, tmp_path / "alignment", str(targets), 100, "descriptors"
    )
    assert first is not None

    second = create_aligned_stage_config(
        first[0], target_run, tmp_path / "alignment", str(targets), 50, "descriptors"
    )

    assert second is not None
    summary = yaml.safe_load(second[2].read_text())
    assert summary["target_coverage_percent"] == 50
    assert summary["stages"]["descriptors"]["status"] == "ready"


def test_structural_target_coverage_preserves_numeric_parameters():
    config = {
        "calculate_ring_infraction": True,
        "ring_infraction_hetcycle_min_size": 4,
        "calculate_stereo_center": True,
        "stereo_max_centers": 4,
        "stereo_max_undefined": 2,
        "calculate_halogenicity": True,
        "halogenicity_thresh_F": 6,
        "halogenicity_thresh_Br": 3,
        "halogenicity_thresh_Cl": 3,
        "calculate_symmetry": True,
        "symmetry_threshold": 0.8,
    }
    masks = pd.DataFrame(
        {
            "model_name": ["target"] * 10,
            "mol_idx": list(range(10)),
            "smiles": ["C"] * 10,
            "ring_infraction": [True] * 10,
            "stereo_center": [True] * 10,
            "halogenicity": [True] * 10,
            "symmetry": [True] * 10,
            "_metric_stereo_centers": list(range(10)),
            "_metric_stereo_undefined": [0] * 10,
            "_metric_halogen_F": list(range(10)),
            "_metric_halogen_Br": [0] * 10,
            "_metric_halogen_Cl": [0] * 10,
            "_metric_symmetry": [value / 10 for value in range(1, 11)],
            "_metric_ring_hard_failure": [False] * 10,
            "_metric_ring_problem_size": [3, 3, 4, 4, 5, 5, 6, 6, 7, 7],
        }
    )

    source_config = yaml.safe_load(yaml.safe_dump(config))
    audit = _align_structural_filter_config(config, masks, 80)

    assert audit["parameters"] == {}
    assert audit["retained_molecules"] == 10
    assert audit["policy"] == "source_config_preserved"
    assert config == source_config


def test_structural_audit_uses_source_hard_rules_without_changing_policy():
    config = {
        "calculate_common_alerts": True,
        "filter_common_alerts": True,
        "include_rulesets": ["PAINS", "Dundee"],
        "common_alerts_filter_include_rulesets": ["PAINS"],
        "common_alerts_filter_exclude_rulesets": [],
        "calculate_molgraph_stats": True,
        "filter_molgraph_stats": True,
        "calculate_NIBR": True,
        "filter_NIBR": True,
        "calculate_lilly": True,
        "filter_lilly": True,
        "calculate_bredt": True,
        "filter_bredt": False,
        "calculate_stereo_center": True,
        "filter_stereo_center": False,
        "filter_undefined_stereo_center": True,
    }
    masks = pd.DataFrame(
        {
            "model_name": ["target"] * 10,
            "mol_idx": list(range(10)),
            "smiles": ["CCO"] * 10,
            "common_alerts:PAINS": [False, *([True] * 9)],
            "common_alerts:Dundee": [False, False, *([True] * 8)],
            "molgraph_stats": [True] * 10,
            "NIBR": [True, False, *([True] * 8)],
            "lilly": [True, True, False, *([True] * 7)],
            "bredt": [False] * 10,
            "stereo_center": [False] * 10,
            "undefined_stereo_center": [True] * 10,
        }
    )

    source_config = yaml.safe_load(yaml.safe_dump(config))
    audit = _align_structural_filter_config(config, masks, 80)

    assert audit["retained_molecules"] == 7
    assert audit["coverage_met"] is False
    assert set(audit["enabled_rules"]) == {
        "common_alerts:PAINS",
        "molgraph_stats",
        "NIBR",
        "lilly",
        "undefined_stereo_center",
    }
    assert audit["disabled_rules"] == []
    assert "common_alerts:Dundee" not in audit["rules"]
    assert "bredt" not in audit["rules"]
    assert "stereo_center" not in audit["rules"]
    assert config == source_config


def test_generated_structural_config_preserves_source_hard_policy(tmp_path):
    master = _base_master(tmp_path)
    struct_path = Path(master["config_structFilters"])
    struct = yaml.safe_load(struct_path.read_text())
    struct.update(
        {
            "calculate_common_alerts": True,
            "filter_common_alerts": True,
            "include_rulesets": ["PAINS", "Dundee"],
            "common_alerts_filter_include_rulesets": ["PAINS"],
            "common_alerts_filter_exclude_rulesets": [],
            "calculate_NIBR": True,
            "filter_NIBR": True,
            "calculate_lilly": True,
            "filter_lilly": True,
            "calculate_stereo_center": True,
            "filter_stereo_center": False,
            "filter_undefined_stereo_center": True,
        }
    )
    _write_yaml(struct_path, struct)

    stage_dir = tmp_path / "target_run/stages/03_structural_filters_post"
    identities = {
        "smiles": ["CCO"] * 10,
        "model_name": ["target"] * 10,
        "mol_idx": [f"mol-{index}" for index in range(10)],
    }
    for name, values in {
        "NIBR": [True, False, *([True] * 8)],
        "lilly": [True, True, False, *([True] * 7)],
    }.items():
        output = stage_dir / name
        output.mkdir(parents=True)
        pd.DataFrame({**identities, "pass": values}).to_csv(
            output / "extended.csv", index=False
        )
    common = stage_dir / "common_alerts"
    common.mkdir(parents=True)
    pd.DataFrame(
        {
            **identities,
            "pass_PAINS": [False, *([True] * 9)],
            "pass_Dundee": [False, False, *([True] * 8)],
            "pass_any": [True] * 10,
        }
    ).to_csv(common / "extended.csv", index=False)
    stereo = stage_dir / "stereo_center"
    stereo.mkdir(parents=True)
    pd.DataFrame(
        {
            **identities,
            "pass": [True] * 10,
            "undefined_stereo_pass": [True] * 10,
        }
    ).to_csv(stereo / "extended.csv", index=False)

    result = create_aligned_stage_config(
        master,
        tmp_path / "target_run",
        tmp_path / "alignment",
        str(tmp_path / "targets.csv"),
        80,
        "struct_filters",
    )

    assert result is not None
    aligned, _master_path, thresholds_path = result
    generated = yaml.safe_load(Path(aligned["config_structFilters"]).read_text())
    audit = yaml.safe_load(thresholds_path.read_text())["stages"]["struct_filters"]
    assert generated == struct
    assert generated["filter_lilly"] is True
    assert audit["status"] == "source_config_preserved"
    assert audit["thresholds"]["retained_molecules"] == 7
    assert audit["thresholds"]["coverage_met"] is False


def test_structural_metrics_keep_undefined_stereo_as_separate_policy(tmp_path):
    stage_dir = tmp_path / "03_structural_filters_post"
    stereo_dir = stage_dir / "stereo_center"
    stereo_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "smiles": ["CCO", "CCN"],
            "model_name": ["target", "target"],
            "mol_idx": ["one", "two"],
            "pass": [True, True],
            "undefined_stereo_pass": [True, False],
        }
    ).to_csv(stereo_dir / "extended.csv", index=False)

    masks = _read_structural_rule_masks(stage_dir)

    assert masks is not None
    assert masks["stereo_center"].tolist() == [True, True]
    assert masks["undefined_stereo_center"].tolist() == [True, False]


def test_leftover_synthesis_bounds_mode_is_stripped_from_generated_master(tmp_path):
    master = _base_master(tmp_path)
    master["alignment"]["synthesis_bounds_mode"] = "expand"
    alignment_root = tmp_path / "alignment"
    source_dir = alignment_root / "source_configs"
    source_dir.mkdir(parents=True)
    source_master = source_dir / "source_config.yml"
    source_master.write_text(
        "alignment:\n"
        "  enabled: true\n"
        "  target_coverage_percent: 95\n"
        "  descriptor_bounds_mode: expand\n"
        "  synthesis_bounds_mode: leftover\n"
        "config_synthesis: unused.yml\n",
        encoding="utf-8",
    )
    targets = tmp_path / "targets.csv"
    targets.write_text("smiles\nCCO\n", encoding="utf-8")
    metrics_dir = (
        tmp_path
        / "target_run"
        / "stages"
        / "02_descriptors_initial"
        / "metrics"
    )
    metrics_dir.mkdir(parents=True)
    pd.DataFrame({"metric": [1.0]}).to_csv(
        metrics_dir / "descriptors_all.csv", index=False
    )

    result = create_aligned_stage_config(
        master,
        tmp_path / "target_run",
        alignment_root,
        str(targets),
        90,
        "descriptors",
    )

    assert result is not None
    aligned, master_path, thresholds_path = result
    written = yaml.safe_load(master_path.read_text())
    summary = yaml.safe_load(thresholds_path.read_text())
    assert "synthesis_bounds_mode" not in aligned.get("alignment", {})
    assert "synthesis_bounds_mode" not in written.get("alignment", {})
    assert "synthesis_bounds_mode" not in master_path.read_text()
    assert "synthesis_bounds_mode" not in summary
    assert (
        _alignment_target_coverage_from_config(
            {
                "alignment": {
                    "enabled": True,
                    "target_coverage_percent": 95,
                    "synthesis_bounds_mode": "leftover",
                }
            },
            None,
        )
        == 95
    )


def test_structural_alignment_preserves_disabled_source_rules(tmp_path):
    master = _base_master(tmp_path)
    targets = tmp_path / "targets.csv"
    targets.write_text("smiles\nCCO\n", encoding="utf-8")
    target_run = tmp_path / "target_run"
    stage_dir = target_run / "stages" / "03_structural_filters_post"

    identities = {
        "model_name": ["target"] * 10,
        "mol_idx": list(range(10)),
        "smiles": ["CCO"] * 10,
    }
    common_alerts = stage_dir / "common_alerts"
    bredt = stage_dir / "bredt"
    lilly = stage_dir / "lilly"
    common_alerts.mkdir(parents=True)
    bredt.mkdir(parents=True)
    lilly.mkdir(parents=True)

    pd.DataFrame(
        {
            **identities,
            "pass_Rule A": [True] * 10,
            "pass_Rule B": [False, False, *([True] * 8)],
            "pass": [False, False, *([True] * 8)],
            "pass_any": [True] * 10,
        }
    ).to_csv(common_alerts / "extended.csv", index=False)
    pd.DataFrame({**identities, "pass": [True, True, False, *([True] * 7)]}).to_csv(
        bredt / "extended.csv", index=False
    )
    pd.DataFrame(
        {**identities, "pass": [True, True, True, *([False] * 5), True, True]}
    ).to_csv(lilly / "extended.csv", index=False)

    result = create_aligned_stage_config(
        master,
        target_run,
        tmp_path / "alignment",
        str(targets),
        80,
        "struct_filters",
    )

    assert result is not None
    aligned, _master_path, thresholds_path = result
    aligned_dir = tmp_path / "alignment" / "aligned_configs"
    struct_config = yaml.safe_load(Path(aligned["config_structFilters"]).read_text())
    summary = yaml.safe_load(thresholds_path.read_text())
    thresholds = summary["stages"]["struct_filters"]["thresholds"]

    assert (aligned_dir / "config_structFilters.yml").exists()
    assert not (aligned_dir / "config_synthesis.yml").exists()
    assert struct_config["calculate_common_alerts"] is False
    assert struct_config["include_rulesets"] == []
    assert struct_config["calculate_bredt"] is False
    assert struct_config["calculate_lilly"] is False
    assert struct_config["calculate_symmetry"] is False
    assert thresholds["required_retained_molecules"] == 8
    assert thresholds["retained_molecules"] == 10
    assert thresholds["enabled_rules"] == []
    assert thresholds["rules"] == {}
    failure_audit = pd.read_csv(aligned_dir / "structural_filter_failures.csv")
    assert len(failure_audit) == 8
    assert set(failure_audit["rule"]) == {
        "common_alerts:Rule B",
        "bredt",
        "lilly",
    }


def test_docking_config_is_created_from_all_tool_scores(tmp_path):
    master = _base_master(tmp_path)
    targets = tmp_path / "targets.csv"
    targets.write_text("smiles\nCCO\n", encoding="utf-8")
    docking_dir = tmp_path / "target_run" / "stages" / "05_docking"
    docking_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "smiles": ["CCO"] * 10,
            "model_name": ["target"] * 10,
            "mol_idx": [f"mol-{index}" for index in range(10)],
        }
    ).to_csv(docking_dir / "input_molecules.csv", index=False)

    writer = Chem.SDWriter(str(docking_dir / "docking_out.sdf"))
    try:
        offsets = {"smina": -12, "gnina": -11, "matcha": -10}
        for index in range(10):
            for tool, offset in offsets.items():
                mol = Chem.MolFromSmiles("CCO")
                mol.SetProp("mol_idx", f"mol-{index}")
                mol.SetProp("docking_tool", tool)
                mol.SetDoubleProp("minimizedAffinity", offset + index)
                writer.write(mol)
                if tool == "gnina":
                    worse_pose = Chem.MolFromSmiles("CCO")
                    worse_pose.SetProp("mol_idx", f"mol-{index}")
                    worse_pose.SetProp("docking_tool", tool)
                    worse_pose.SetDoubleProp("minimizedAffinity", offset + index + 100)
                    writer.write(worse_pose)
    finally:
        writer.close()

    result = create_aligned_stage_config(
        master,
        tmp_path / "target_run",
        tmp_path / "alignment",
        str(targets),
        80,
        "docking",
    )

    assert result is not None
    aligned, _master_path, thresholds_path = result
    docking_config = yaml.safe_load(Path(aligned["config_docking"]).read_text())
    assert docking_config["calculate_score_thresholds_from_targets"] is True
    assert docking_config["score_thresholds"] == {
        "smina": {"score_property": "minimizedAffinity", "max": -5},
        "gnina": {"score_property": "minimizedAffinity", "max": -4},
        "matcha": {"score_property": "minimizedAffinity", "max": -3},
    }
    summary = yaml.safe_load(thresholds_path.read_text())
    audit = summary["stages"]["docking"]["thresholds"]
    assert summary["stages"]["docking"]["status"] == "ready"
    assert audit["target_molecules"] == 10
    assert audit["required_retained_molecules"] == 8
    assert audit["retained_molecules"] == 8
    assert audit["combination"] == "all_configured_tools_must_pass"
    assert not (
        tmp_path / "alignment" / "aligned_configs" / "config_docking_filters.yml"
    ).exists()


def test_docking_backend_aligns_matcha_with_minimized_affinity(tmp_path):
    master = _base_master(tmp_path)
    docking_config_path = Path(master["config_docking"])
    docking_config = yaml.safe_load(docking_config_path.read_text())
    docking_config["matcha_config"] = {"backend": "docking"}
    _write_yaml(docking_config_path, docking_config)

    targets = tmp_path / "targets.csv"
    targets.write_text("smiles\nCCO\n", encoding="utf-8")
    docking_dir = tmp_path / "target_run" / "stages" / "05_docking"
    docking_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "smiles": ["CCO"] * 10,
            "model_name": ["target"] * 10,
            "mol_idx": [f"mol-{index}" for index in range(10)],
        }
    ).to_csv(docking_dir / "input_molecules.csv", index=False)

    writer = Chem.SDWriter(str(docking_dir / "docking_out.sdf"))
    try:
        for index in range(10):
            for tool, score in {
                "smina": -12 + index,
                "gnina": -11 + index,
            }.items():
                mol = Chem.MolFromSmiles("CCO")
                mol.SetProp("mol_idx", f"mol-{index}")
                mol.SetProp("docking_tool", tool)
                mol.SetDoubleProp("minimizedAffinity", score)
                writer.write(mol)

            matcha = Chem.MolFromSmiles("CCO")
            matcha.SetProp("mol_idx", f"mol-{index}")
            matcha.SetProp("docking_tool", "matcha")
            matcha.SetDoubleProp("balmus_score", -300 + index * 10)
            matcha.SetDoubleProp("minimizedAffinity", -10 + index)
            writer.write(matcha)
    finally:
        writer.close()

    result = create_aligned_stage_config(
        master,
        tmp_path / "target_run",
        tmp_path / "alignment",
        str(targets),
        80,
        "docking",
    )

    assert result is not None
    aligned, _master_path, _thresholds_path = result
    aligned_docking = yaml.safe_load(Path(aligned["config_docking"]).read_text())
    assert aligned_docking["score_thresholds"] == {
        "smina": {"score_property": "minimizedAffinity", "max": -5},
        "gnina": {"score_property": "minimizedAffinity", "max": -4},
        "matcha": {"score_property": "minimizedAffinity", "max": -3},
    }


def test_docking_score_alignment_requires_explicit_config_switch(tmp_path):
    master = _base_master(tmp_path)
    docking_config_path = Path(master["config_docking"])
    docking_config = yaml.safe_load(docking_config_path.read_text())
    docking_config["calculate_score_thresholds_from_targets"] = False
    _write_yaml(docking_config_path, docking_config)

    targets = tmp_path / "targets.csv"
    targets.write_text("smiles\nCCO\n", encoding="utf-8")
    probe = create_probe_config(master, str(targets), tmp_path / "alignment")
    probe_docking = yaml.safe_load(Path(probe["config_docking"]).read_text())

    assert probe_docking["run"] is False
    assert (
        create_aligned_stage_config(
            master,
            tmp_path / "target_run",
            tmp_path / "alignment",
            str(targets),
            90,
            "docking",
        )
        is None
    )
    summary = yaml.safe_load(
        (
            tmp_path / "alignment" / "aligned_configs" / "alignment_thresholds.yml"
        ).read_text()
    )
    assert summary["stages"]["docking"]["status"] == "disabled_by_config"


def test_stage_configs_write_coverage_thresholds_without_docking_filters(tmp_path):
    master = _base_master(tmp_path)
    targets = tmp_path / "targets.csv"
    targets.write_text("smiles\nCCO\n", encoding="utf-8")
    target_run = tmp_path / "target_run"

    descriptors_dir = target_run / "stages" / "02_descriptors_initial" / "metrics"
    docking_filters_dir = target_run / "stages" / "06_docking_filters"
    descriptors_dir.mkdir(parents=True)
    docking_filters_dir.mkdir(parents=True)

    values = list(range(101))
    pd.DataFrame(
        {
            "metric": values,
            "n_N_atoms": values,
            "n_NO_atoms": values,
            "Hal": values,
        }
    ).to_csv(descriptors_dir / "descriptors_all.csv", index=False)
    pd.DataFrame(
        {
            "gnina_minimizedAffinity": [-value / 10 for value in values],
            "gnina_CNNscore": [value / 100 for value in values],
            "frac_atoms_outside_box": [value / 100 for value in values],
            "n_hbonds": values,
            "min_conformer_rmsd": [value / 10 for value in values],
        }
    ).to_csv(docking_filters_dir / "metrics.csv", index=False)

    descriptor_result = create_aligned_stage_config(
        master, target_run, tmp_path / "alignment", str(targets), 90, "descriptors"
    )
    docking_filter_result = create_aligned_stage_config(
        master,
        target_run,
        tmp_path / "alignment",
        str(targets),
        90,
        "docking_filters",
    )

    assert descriptor_result is not None
    assert docking_filter_result is None
    aligned, master_path, thresholds_path = descriptor_result

    descriptor_config = yaml.safe_load(Path(aligned["config_descriptors"]).read_text())
    summary = yaml.safe_load(thresholds_path.read_text())

    assert aligned["alignment"]["enabled"] is False
    assert KEY_ALIGNMENT_SKIP_FINAL_DESCRIPTORS not in aligned
    assert descriptor_config["borders"] == {
        "metric_min": -10,
        "metric_max": 200,
    }
    assert descriptor_config["structural_constraints"] == {
        "enabled": True,
        "type_limits": {"Hal": 3},
        "element_limits": {"N": 6},
        "max_n_or_o_atoms": 10,
    }
    assert not (
        tmp_path / "alignment" / "aligned_configs" / "config_docking_filters.yml"
    ).exists()
    assert master_path.exists()
    assert summary["target_coverage_percent"] == 90
    assert summary["selection_method"] == "stage_specific_target_coverage"
    assert summary["stages"]["descriptors"]["note"] == (
        "Min/max ranges use one shared stage-level reference subset."
    )
    assert set(summary["stages"]) == {
        "mol_prep",
        "descriptors",
        "struct_filters",
        "docking",
        "docking_filters",
        "final_descriptors",
    }
    assert summary["stages"]["docking_filters"] == {
        "thresholds": {},
        "status": "not_aligned",
        "note": "Docking alignment calibrates raw score thresholds only.",
    }
    assert summary["stages"]["final_descriptors"] == {
        "thresholds": {},
        "status": "skipped_redundant",
        "note": "Initial descriptor measurements are reused.",
    }
