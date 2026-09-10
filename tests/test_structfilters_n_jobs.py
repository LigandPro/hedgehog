"""Structural filters stage n_jobs is declared and honored for all shipped profiles."""

from pathlib import Path
from unittest.mock import MagicMock

import pandas as pd
import pytest
import yaml

from hedgehog._constants import CFG_STRUCT_FILTERS
from hedgehog.struct_filters import main as structfilters_main
from hedgehog.struct_filters import utils as structfilters_utils
from hedgehog.utils.parallel import resolve_n_jobs

CONFIGS_DIR = Path(__file__).resolve().parents[1] / "src" / "hedgehog" / "configs"

PROFILE_FILES = (
    "config_structFilters.yml",
    "config_structFilters_exploration.yml",
    "config_structFilters_strict.yml",
)


@pytest.mark.parametrize("profile_name", PROFILE_FILES)
def test_shipped_struct_filter_profiles_declare_n_jobs(profile_name):
    cfg = yaml.safe_load((CONFIGS_DIR / profile_name).read_text(encoding="utf-8"))
    assert cfg["n_jobs"] == 16
    # Stage value wins over a conflicting master setting.
    assert resolve_n_jobs(cfg, {"n_jobs": 4}) == 16
    opts = structfilters_main._resolve_stage_options(cfg, {"n_jobs": 4})
    assert opts["n_jobs"] == 16


def test_main_passes_stage_n_jobs_to_prepare(tmp_path, monkeypatch):
    input_df = pd.DataFrame(
        {
            "smiles": ["CCO", "CCN"],
            "model_name": ["m1", "m1"],
            "mol_idx": [0, 1],
        }
    )
    input_csv = tmp_path / "input.csv"
    input_df.to_csv(input_csv, index=False)

    struct_cfg = tmp_path / "config_structFilters.yml"
    struct_cfg.write_text(
        yaml.safe_dump(
            {
                "n_jobs": 7,
                "filter_data": False,
                "calculate_bredt": True,
                "write_per_filter_outputs": False,
                "generate_plots": False,
                "generate_failure_analysis": False,
            }
        ),
        encoding="utf-8",
    )

    config = {
        "n_jobs": 99,
        "sample_size": None,
        "folder_to_save": str(tmp_path),
        "generated_mols_path": str(input_csv),
        CFG_STRUCT_FILTERS: str(struct_cfg),
    }

    prepared_payload = {
        "mols": [object(), object()],
        "smiles_model_mols": [
            ("CCO", None, object(), 0),
            ("CCN", None, object(), 1),
        ],
        "base_df": input_df.copy(),
    }
    seen = {}

    def fake_prepare(df, subsample, n_jobs, progress_cb=None):
        del df, subsample, progress_cb
        seen["n_jobs"] = n_jobs
        return prepared_payload

    def fake_process(config_arg, payload, apply_filter, progress_cb=None):
        del config_arg, payload, apply_filter, progress_cb
        return pd.DataFrame(
            {
                "smiles": ["CCO", "CCN"],
                "model_name": ["m1", "m1"],
                "mol_idx": [0, 1],
                "mol": [object(), object()],
                "pass": [True, True],
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
    assert seen["n_jobs"] == 7


def test_common_alerts_resolve_uses_stage_n_jobs(monkeypatch):
    captured = {}

    def _fake_parallel_map(func, items, n_jobs, **kwargs):
        del func, items, kwargs
        captured["n_jobs"] = n_jobs
        return []

    monkeypatch.setattr(structfilters_utils, "parallel_map", _fake_parallel_map)
    monkeypatch.setattr(
        structfilters_utils,
        "filter_alerts",
        lambda _cfg: pd.DataFrame(
            {
                "rule_set_name": ["PAINS"],
                "smarts": ["[#6]"],
                "description": ["x"],
                "mincount": [1],
            }
        ),
    )
    monkeypatch.setattr(
        structfilters_utils,
        "load_config",
        lambda _path: {"n_jobs": 5, "calculate_common_alerts": True},
    )
    monkeypatch.setattr(structfilters_utils, "_compile_alert_smarts", lambda _data: [])
    monkeypatch.setattr(
        structfilters_utils,
        "_setup_alerts_progress",
        lambda *_a, **_k: (None, MagicMock(), None),
    )
    monkeypatch.setattr(
        structfilters_utils,
        "_build_alerts_results",
        lambda *_a, **_k: pd.DataFrame(),
    )

    structfilters_utils.apply_structural_alerts(
        {CFG_STRUCT_FILTERS: "unused.yml", "n_jobs": 99},
        mols=["m1", "m2"],
    )
    assert captured["n_jobs"] == 5
