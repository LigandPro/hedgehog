from pathlib import Path

import pandas as pd
import pytest
import yaml
from rdkit import Chem

from hedgehog._constants import KEY_ALIGNMENT_SKIP_FINAL_DESCRIPTORS
from hedgehog.config_alignment import (
    _align_descriptor_config,
    _align_structural_filter_config,
    _yaml_number,
    create_aligned_stage_config,
    create_probe_config,
    set_probe_molprep_allowed_atoms,
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
            "pose_quality": {"max_clashes": 2, "max_strain_energy": 50},
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
            "exclude_descriptions": {"Rule A": ["ignored"]},
        },
    )

    return {
        "generated_mols_path": str(tmp_path / "candidates.csv"),
        "target_mols_path": str(tmp_path / "targets.csv"),
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


def test_descriptor_alignment_reads_ring_size_lists_without_replacing_them():
    config = {"borders": {"ring_size_min": 3, "ring_size_max": 12}}
    metrics = pd.DataFrame(
        {"ring_size": ["[3, 6]", "[4, 7]", "[5, 8]", "[6, 9]", "[]"]}
    )
    original = metrics["ring_size"].copy()

    thresholds = _align_descriptor_config(config, metrics, 80)

    assert metrics["ring_size"].equals(original)
    assert thresholds == {"ring_size_min": 3, "ring_size_max": 8}
    assert config["borders"] == {"ring_size_min": 3, "ring_size_max": 8}


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


def test_probe_config_keeps_thresholds_and_disables_filtering_explicitly(tmp_path):
    master = _base_master(tmp_path)
    targets = tmp_path / "targets.csv"
    targets.write_text("smiles\nCCO\n", encoding="utf-8")

    probe = create_probe_config(master, str(targets), tmp_path / "alignment")
    allowed_atoms = set_probe_molprep_allowed_atoms(
        probe,
        pd.DataFrame({"smiles": ["CCO", "C[SiH3]", "[Na+].[Cl-]", "[H]C"]}),
    )

    assert probe[KEY_ALIGNMENT_SKIP_FINAL_DESCRIPTORS] is True

    original_molprep = yaml.safe_load(Path(master["config_mol_prep"]).read_text())
    probe_molprep = yaml.safe_load(Path(probe["config_mol_prep"]).read_text())
    original_descriptors = yaml.safe_load(
        Path(master["config_descriptors"]).read_text()
    )
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
    assert probe_synthesis["filter_solved_only"] is False
    assert probe_synthesis["run_retrosynthesis"] is False
    assert probe_synthesis["enabled_scores"] is None
    assert probe_synthesis["alignment_measurement_mode"] is True
    assert probe_synthesis["apply_score_filters"] is False
    assert probe_synthesis["sa_score_min"] == 1
    assert probe_synthesis["sa_score_max"] == 10
    assert probe_synthesis["score_filters"]["sync_score"] == {
        "min": 0.5,
        "max": 1,
    }
    assert Path(probe["config_synthesis"]).parent.name == (
        "calibration_configs_unfiltered"
    )
    assert Path(probe["folder_to_save"]).name == "calibration_target_run"
    source_master = tmp_path / "alignment" / "source_configs" / "source_config.yml"
    assert source_master.is_file()
    source = yaml.safe_load(source_master.read_text())
    source_synthesis = yaml.safe_load(Path(source["config_synthesis"]).read_text())
    assert source_synthesis["run_retrosynthesis"] is True
    assert source_synthesis["sa_score_min"] == 1
    assert probe_struct["run"] is True
    assert probe_struct["filter_data"] is False
    assert probe_struct["write_per_filter_outputs"] is True
    assert probe_struct["generate_plots"] is False
    assert probe_struct["calculate_common_alerts"] is True
    assert probe_struct["calculate_bredt"] is True
    assert probe_struct["calculate_lilly"] is True
    assert probe_struct["calculate_symmetry"] is True
    assert probe_struct["include_rulesets"] == ["Rule A", "Rule B"]
    assert probe_struct["exclude_descriptions"] == {}
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
    assert first_summary["stages"]["synthesis"]["status"] == "pending_metrics"

    synthesis_dir = target_run / "stages" / "04_synthesis"
    synthesis_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "sa_score": [1.0, 2.0, 3.0],
            "sync_score": [0.1, 0.2, 0.3],
        }
    ).to_csv(synthesis_dir / "synthesis_scores.csv", index=False)

    synthesis_result = create_aligned_stage_config(
        master,
        target_run,
        tmp_path / "alignment",
        str(targets),
        100,
        "synthesis",
    )

    assert synthesis_result is not None
    assert (aligned_dir / "config_descriptors.yml").exists()
    assert (aligned_dir / "config_synthesis.yml").exists()
    assert not (aligned_dir / "config_docking.yml").exists()
    second_summary = yaml.safe_load(
        (aligned_dir / "alignment_thresholds.yml").read_text()
    )
    assert second_summary["stages"]["descriptors"]["status"] == "ready"
    assert second_summary["stages"]["synthesis"]["status"] == "ready"


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


def test_synthesis_alignment_uses_only_source_enabled_scorers(tmp_path):
    master = _base_master(tmp_path)
    synthesis_config = yaml.safe_load(Path(master["config_synthesis"]).read_text())
    synthesis_config["enabled_scores"] = ["gasa"]
    synthesis_config["score_filters"]["gasa_score"] = {"min": 0.2, "max": 0.8}
    _write_yaml(Path(master["config_synthesis"]), synthesis_config)

    targets = tmp_path / "targets.csv"
    targets.write_text("smiles\nCCO\n", encoding="utf-8")
    synthesis_dir = tmp_path / "target_run" / "stages" / "04_synthesis"
    synthesis_dir.mkdir(parents=True)
    values = list(range(101))
    pd.DataFrame(
        {
            "sa_score": values,
            "sync_score": [value / 100 for value in reversed(values)],
            "gasa_score": values,
        }
    ).to_csv(synthesis_dir / "synthesis_scores.csv", index=False)

    result = create_aligned_stage_config(
        master,
        tmp_path / "target_run",
        tmp_path / "alignment",
        str(targets),
        90,
        "synthesis",
    )

    assert result is not None
    aligned, _master_path, thresholds_path = result
    generated = yaml.safe_load(Path(aligned["config_synthesis"]).read_text())
    summary = yaml.safe_load(thresholds_path.read_text())
    assert generated["enabled_scores"] == ["gasa"]
    assert generated["sa_score_min"] == 1
    assert generated["sa_score_max"] == 10
    assert generated["score_filters"]["sync_score"] == {"min": 0.5, "max": 1}
    assert generated["score_filters"]["gasa_score"] == {"min": 5, "max": 95}
    thresholds = summary["stages"]["synthesis"]["thresholds"]
    assert thresholds == {"gasa_score": {"min": 5, "max": 95}}


def test_synthesis_alignment_disables_unavailable_target_score_criteria(tmp_path):
    master = _base_master(tmp_path)
    synthesis_config = yaml.safe_load(Path(master["config_synthesis"]).read_text())
    synthesis_config["enabled_scores"] = ["sa", "sync"]
    _write_yaml(Path(master["config_synthesis"]), synthesis_config)
    targets = tmp_path / "targets.csv"
    targets.write_text("smiles\nCCO\n", encoding="utf-8")
    synthesis_dir = tmp_path / "target_run" / "stages" / "04_synthesis"
    synthesis_dir.mkdir(parents=True)
    pd.DataFrame({"sa_score": [1.0, 2.0, 3.0]}).to_csv(
        synthesis_dir / "synthesis_scores.csv", index=False
    )

    result = create_aligned_stage_config(
        master,
        tmp_path / "target_run",
        tmp_path / "alignment",
        str(targets),
        90,
        "synthesis",
    )

    assert result is not None
    aligned, _master_path, _thresholds_path = result
    synthesis_config = yaml.safe_load(Path(aligned["config_synthesis"]).read_text())
    assert synthesis_config["run_retrosynthesis"] is True
    assert synthesis_config["score_filters"]["sync_score"] == {
        "min": None,
        "max": None,
    }


def test_structural_numeric_parameters_are_aligned_before_rule_selection():
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

    audit = _align_structural_filter_config(config, masks, 80)

    assert audit["parameters"] == {
        "stereo_max_centers": 8,
        "stereo_max_undefined": 1,
        "halogenicity_thresh_F": 7,
        "halogenicity_thresh_Br": 0,
        "halogenicity_thresh_Cl": 0,
        "symmetry_threshold": 0.8,
        "ring_infraction_hetcycle_min_size": 3,
    }
    assert audit["retained_molecules"] == 8
    assert config["calculate_ring_infraction"] is False
    assert config["calculate_stereo_center"] is True
    assert config["calculate_halogenicity"] is True
    assert config["calculate_symmetry"] is True


def test_structural_config_turns_rules_on_and_off_for_combined_percentile(tmp_path):
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
    assert struct_config["calculate_common_alerts"] is True
    assert struct_config["include_rulesets"] == ["Rule A"]
    assert struct_config["calculate_bredt"] is True
    assert struct_config["calculate_lilly"] is False
    assert struct_config["calculate_symmetry"] is False
    assert thresholds["required_retained_molecules"] == 8
    assert thresholds["retained_molecules"] == 9
    assert thresholds["rules"]["common_alerts:Rule B"] == {
        "enabled": False,
        "failed_molecules": 2,
        "failed_percent": 20.0,
        "combined_retained_if_enabled": 7,
    }
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
    synthesis_dir = target_run / "stages" / "04_synthesis"
    docking_filters_dir = target_run / "stages" / "06_docking_filters"
    descriptors_dir.mkdir(parents=True)
    synthesis_dir.mkdir(parents=True)
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
            "sa_score": values,
            "sync_score": [value / 100 for value in values],
        }
    ).to_csv(synthesis_dir / "synthesis_scores.csv", index=False)
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
    synthesis_result = create_aligned_stage_config(
        master, target_run, tmp_path / "alignment", str(targets), 90, "synthesis"
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
    assert synthesis_result is not None
    assert docking_filter_result is None
    aligned, master_path, thresholds_path = synthesis_result

    descriptor_config = yaml.safe_load(Path(aligned["config_descriptors"]).read_text())
    synthesis_config = yaml.safe_load(Path(aligned["config_synthesis"]).read_text())
    summary = yaml.safe_load(thresholds_path.read_text())

    assert aligned["alignment"]["enabled"] is False
    assert KEY_ALIGNMENT_SKIP_FINAL_DESCRIPTORS not in aligned
    assert descriptor_config["borders"] == {
        "metric_min": 5,
        "metric_max": 95,
        "Hal_min": 5,
        "Hal_max": 95,
        "n_NO_atoms_min": 5,
        "n_NO_atoms_max": 95,
        "n_N_atoms_min": 5,
        "n_N_atoms_max": 95,
    }
    assert descriptor_config["structural_constraints"] == {
        "enabled": False,
        "type_limits": {"Hal": 3},
        "element_limits": {"N": 6},
        "max_n_or_o_atoms": 10,
    }
    assert synthesis_config["enabled_scores"] is None
    assert synthesis_config["sa_score_min"] == 5
    assert synthesis_config["sa_score_max"] == 95
    assert synthesis_config["score_filters"]["sync_score"] == {
        "min": 0.5,
        "max": 1,
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
        "synthesis",
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
