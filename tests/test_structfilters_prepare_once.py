from pathlib import Path

import pandas as pd
import yaml

from hedgehog.stages.structFilters import main as structfilters_main


def _write_yaml(path: Path, payload: dict) -> None:
    path.write_text(yaml.safe_dump(payload), encoding="utf-8")


def test_prepare_payload_called_once_for_multiple_filters(tmp_path, monkeypatch):
    input_df = pd.DataFrame(
        {
            "smiles": ["CCO", "CCN", "CCC"],
            "model_name": ["m1", "m1", "m1"],
            "mol_idx": [0, 1, 2],
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
            "calculate_halogenicity": True,
            "parse_input_n_jobs": 1,
            "write_per_filter_outputs": False,
            "generate_plots": False,
            "generate_failure_analysis": False,
            "combine_in_memory": True,
        },
    )

    config = {
        "sample_size": None,
        "folder_to_save": str(tmp_path),
        "generated_mols_path": str(input_csv),
        "config_structFilters": str(struct_cfg),
    }

    call_counts = {"prepare": 0, "process": 0}
    prepared_payload = {
        "mols": [object(), object(), object()],
        "smiles_model_mols": [
            ("CCO", None, object(), 0),
            ("CCN", None, object(), 1),
            ("CCC", None, object(), 2),
        ],
        "base_df": input_df.copy(),
    }

    def fake_prepare(df, subsample, parse_n_jobs, progress_cb=None):
        del subsample, parse_n_jobs, progress_cb
        call_counts["prepare"] += 1
        assert len(df) == 3
        return prepared_payload

    def fake_process(config, payload, apply_filter, progress_cb=None):
        del config, payload, apply_filter, progress_cb
        call_counts["process"] += 1
        return pd.DataFrame(
            {
                "smiles": ["CCO", "CCN", "CCC"],
                "model_name": ["m1", "m1", "m1"],
                "mol_idx": [0, 1, 2],
                "mol": [object(), object(), object()],
                "pass": [True, False, True],
            }
        )

    def fake_get_basic_stats(
        config_struct_filters, filter_results, model_name, filter_name
    ):
        del config_struct_filters, model_name, filter_name
        metrics = pd.DataFrame({"model_name": ["m1"], "num_mol": [len(filter_results)]})
        return metrics, filter_results.copy()

    monkeypatch.setattr(structfilters_main, "prepare_structfilters_input", fake_prepare)
    monkeypatch.setattr(structfilters_main, "process_prepared_payload", fake_process)
    monkeypatch.setattr(structfilters_main, "get_basic_stats", fake_get_basic_stats)
    monkeypatch.setattr(
        structfilters_main, "_save_filter_results", lambda *args, **kwargs: None
    )
    monkeypatch.setattr(
        structfilters_main,
        "inject_identity_columns_to_all_csvs",
        lambda *args, **kwargs: None,
    )

    structfilters_main.main(config, "stages/03_structural_filters_post")

    assert call_counts["prepare"] == 1
    assert call_counts["process"] == 2
