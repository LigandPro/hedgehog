from pathlib import Path
from unittest.mock import MagicMock

import pandas as pd
import yaml

from hedgehog.struct_filters import main as structfilters_main


def _write_yaml(path: Path, payload: dict) -> None:
    path.write_text(yaml.safe_dump(payload), encoding="utf-8")


def _build_base_config(tmp_path: Path, struct_cfg_payload: dict) -> dict:
    input_df = pd.DataFrame({"smiles": ["CCO"], "model_name": ["m1"], "mol_idx": [0]})
    input_csv = tmp_path / "input.csv"
    input_df.to_csv(input_csv, index=False)

    struct_cfg = tmp_path / "config_structFilters.yml"
    _write_yaml(struct_cfg, struct_cfg_payload)

    return {
        "sample_size": None,
        "folder_to_save": str(tmp_path),
        "generated_mols_path": str(input_csv),
        "config_structFilters": str(struct_cfg),
    }


def _mock_filter_processing(monkeypatch):
    prepared_payload = {
        "mols": [object()],
        "smiles_model_mols": [("CCO", None, object(), 0)],
        "base_df": pd.DataFrame(
            {"smiles": ["CCO"], "model_name": ["m1"], "mol_idx": [0]}
        ),
    }
    filter_df = pd.DataFrame(
        {
            "smiles": ["CCO"],
            "model_name": ["m1"],
            "mol_idx": [0],
            "mol": [object()],
            "pass": [True],
        }
    )
    metrics_df = pd.DataFrame(
        {"model_name": ["m1"], "num_mol": [1], "banned_ratio": [0.0]}
    )

    monkeypatch.setattr(
        structfilters_main,
        "prepare_structfilters_input",
        lambda *args, **kwargs: prepared_payload,
    )
    monkeypatch.setattr(
        structfilters_main,
        "process_prepared_payload",
        lambda *args, **kwargs: filter_df.copy(),
    )
    monkeypatch.setattr(
        structfilters_main,
        "get_basic_stats",
        lambda *args, **kwargs: (metrics_df.copy(), filter_df.copy()),
    )
    monkeypatch.setattr(
        structfilters_main,
        "inject_identity_columns_to_all_csvs",
        lambda *args, **kwargs: None,
    )


def test_flags_disable_outputs_and_plots(tmp_path, monkeypatch):
    _mock_filter_processing(monkeypatch)
    config = _build_base_config(
        tmp_path,
        {
            "filter_data": True,
            "calculate_bredt": True,
            "enforced_filters": [],
            "write_per_filter_outputs": False,
            "generate_plots": True,
            "generate_failure_analysis": True,
        },
    )

    save_mock = MagicMock()
    plot_stats_mock = MagicMock()
    plot_ratio_mock = MagicMock()
    fail_analysis_mock = MagicMock()
    combine_mock = MagicMock()

    monkeypatch.setattr(structfilters_main, "_save_filter_results", save_mock)
    monkeypatch.setattr(structfilters_main, "plot_calculated_stats", plot_stats_mock)
    monkeypatch.setattr(structfilters_main, "plot_restriction_ratios", plot_ratio_mock)
    monkeypatch.setattr(
        structfilters_main, "plot_filter_failures_analysis", fail_analysis_mock
    )
    monkeypatch.setattr(
        structfilters_main, "combine_filter_results_in_memory", combine_mock
    )

    structfilters_main.main(config, "StructFilters")

    save_mock.assert_not_called()
    plot_stats_mock.assert_not_called()
    plot_ratio_mock.assert_not_called()
    fail_analysis_mock.assert_not_called()
    combine_mock.assert_called_once()
    assert combine_mock.call_args.args[2] == {}


def test_flags_enable_outputs_and_plots(tmp_path, monkeypatch):
    _mock_filter_processing(monkeypatch)
    config = _build_base_config(
        tmp_path,
        {
            "filter_data": False,
            "calculate_bredt": True,
            "write_per_filter_outputs": True,
            "generate_plots": True,
            "generate_failure_analysis": True,
        },
    )

    save_mock = MagicMock()
    plot_stats_mock = MagicMock()
    plot_ratio_mock = MagicMock()
    fail_analysis_mock = MagicMock()

    monkeypatch.setattr(structfilters_main, "_save_filter_results", save_mock)
    monkeypatch.setattr(structfilters_main, "plot_calculated_stats", plot_stats_mock)
    monkeypatch.setattr(structfilters_main, "plot_restriction_ratios", plot_ratio_mock)
    monkeypatch.setattr(
        structfilters_main, "plot_filter_failures_analysis", fail_analysis_mock
    )
    monkeypatch.setattr(
        structfilters_main,
        "combine_filter_results_in_memory",
        lambda *args, **kwargs: None,
    )

    structfilters_main.main(config, "StructFilters")

    save_mock.assert_called_once()
    plot_stats_mock.assert_called_once()
    plot_ratio_mock.assert_called_once()
    fail_analysis_mock.assert_called_once()


def test_progress_uses_real_molecule_totals(tmp_path, monkeypatch):
    _mock_filter_processing(monkeypatch)

    input_df = pd.DataFrame(
        {
            "smiles": ["CCO", "CCC"],
            "model_name": ["m1", "m1"],
            "mol_idx": [0, 1],
        }
    )
    input_csv = tmp_path / "input.csv"
    input_df.to_csv(input_csv, index=False)

    struct_cfg = tmp_path / "config_structFilters.yml"
    _write_yaml(
        struct_cfg,
        {
            "filter_data": False,
            "calculate_bredt": True,
            "write_per_filter_outputs": False,
            "generate_plots": False,
            "generate_failure_analysis": False,
        },
    )

    config = {
        "sample_size": None,
        "folder_to_save": str(tmp_path),
        "generated_mols_path": str(input_csv),
        "config_structFilters": str(struct_cfg),
    }

    monkeypatch.setattr(
        structfilters_main,
        "combine_filter_results_in_memory",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        structfilters_main, "_save_filter_results", lambda *args, **kwargs: None
    )

    reporter = MagicMock()
    structfilters_main.main(config, "StructFilters", reporter=reporter)

    assert reporter.progress.call_count >= 2
    totals = [call.args[1] for call in reporter.progress.call_args_list]
    currents = [call.args[0] for call in reporter.progress.call_args_list]

    assert all(total == 2 for total in totals)
    assert 0 in currents
    assert 2 in currents


def test_enforced_filters_keep_other_calculations_diagnostic(tmp_path, monkeypatch):
    _mock_filter_processing(monkeypatch)
    config = _build_base_config(
        tmp_path,
        {
            "filter_data": True,
            "calculate_bredt": True,
            "calculate_protecting_groups": True,
            "enforced_filters": ["bredt"],
            "write_per_filter_outputs": False,
            "generate_plots": False,
            "generate_failure_analysis": False,
        },
    )
    captured = {}

    def capture_combine(_output_dir, _input_df, pass_masks):
        captured.update(pass_masks)

    monkeypatch.setattr(
        structfilters_main, "combine_filter_results_in_memory", capture_combine
    )

    structfilters_main.main(config, "StructFilters")

    assert set(captured) == {"bredt"}


def test_per_filter_flags_control_survival_independently_of_calculation(
    tmp_path, monkeypatch
):
    _mock_filter_processing(monkeypatch)
    config = _build_base_config(
        tmp_path,
        {
            "filter_data": True,
            "calculate_bredt": True,
            "calculate_protecting_groups": True,
            "filter_bredt": False,
            "filter_protecting_groups": True,
            "write_per_filter_outputs": False,
            "generate_plots": False,
            "generate_failure_analysis": False,
        },
    )
    captured = {}

    def capture_combine(_output_dir, _input_df, pass_masks):
        captured.update(pass_masks)

    monkeypatch.setattr(
        structfilters_main, "combine_filter_results_in_memory", capture_combine
    )

    structfilters_main.main(config, "StructFilters")

    assert set(captured) == {"protecting_groups"}
    profile = pd.read_csv(
        tmp_path / "StructFilters" / "structural_liability_profile.csv"
    )
    assert bool(profile["bredt__filter_enabled"].iloc[0]) is False
    assert bool(profile["protecting_groups__filter_enabled"].iloc[0]) is True


def test_undefined_stereo_is_hard_while_total_is_diagnostic(tmp_path, monkeypatch):
    _mock_filter_processing(monkeypatch)
    stereo_df = pd.DataFrame(
        {
            "smiles": ["CCO"],
            "model_name": ["m1"],
            "mol_idx": [0],
            "mol": [object()],
            "pass": [True],
            "n_stereo_centers": [3],
            "n_undefined_stereo_centers": [3],
            "stereo_max_centers": [4],
            "stereo_max_undefined": [2],
            "undefined_stereo_pass": [False],
            "undefined_stereo_reason": ["undefined_stereocenters=3 > maximum=2"],
        }
    )
    metrics_df = pd.DataFrame(
        {"model_name": ["m1"], "num_mol": [1], "banned_ratio": [0.0]}
    )
    monkeypatch.setattr(
        structfilters_main,
        "get_basic_stats",
        lambda *args, **kwargs: (metrics_df.copy(), stereo_df.copy()),
    )
    config = _build_base_config(
        tmp_path,
        {
            "filter_data": True,
            "calculate_stereo_center": True,
            "filter_stereo_center": False,
            "filter_undefined_stereo_center": True,
            "write_per_filter_outputs": False,
            "write_structural_liability_profile": True,
            "generate_plots": False,
            "generate_failure_analysis": False,
        },
    )
    captured = {}

    def capture_combine(_output_dir, _input_df, pass_masks):
        captured.update(pass_masks)

    monkeypatch.setattr(
        structfilters_main, "combine_filter_results_in_memory", capture_combine
    )

    structfilters_main.main(config, "StructFilters")

    assert set(captured) == {"undefined_stereo_center"}
    assert captured["undefined_stereo_center"]["pass"].tolist() == [False]
    profile = pd.read_csv(
        tmp_path / "StructFilters" / "structural_liability_profile.csv"
    )
    assert bool(profile["stereo_center__filter_enabled"].iloc[0]) is False
    assert bool(profile["undefined_stereo_center__filter_enabled"].iloc[0]) is True
    assert bool(profile["stage3_hard_pass"].iloc[0]) is False
    assert profile["hard_failed_filters"].iloc[0] == "undefined_stereo_center"
