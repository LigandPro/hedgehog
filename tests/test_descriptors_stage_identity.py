"""Regression tests for descriptors stage identity preparation."""

from pathlib import Path

import pandas as pd
import yaml

from hedgehog.descriptors import stage as descriptors_stage


def _write_yaml(path: Path, payload: dict) -> None:
    path.write_text(yaml.safe_dump(payload), encoding="utf-8")


def _build_config(tmp_path: Path) -> dict:
    cfg_path = tmp_path / "config_descriptors.yml"
    _write_yaml(cfg_path, {"run": True, "filter_data": False, "borders": {}})
    return {
        "folder_to_save": str(tmp_path),
        "config_descriptors": str(cfg_path),
    }


def test_stage_adds_model_name_and_mol_idx_for_minimal_input(tmp_path, monkeypatch):
    config = _build_config(tmp_path)
    data = pd.DataFrame({"smiles": ["CCO", "CCC"]})
    captured = {}

    def _fake_compute_metrics(df, *args, **kwargs):
        captured["input_df"] = df.copy()
        return df.copy()

    monkeypatch.setattr(descriptors_stage, "compute_metrics", _fake_compute_metrics)
    monkeypatch.setattr(descriptors_stage, "run_waves", lambda *_args, **_kwargs: None)

    descriptors_stage.run(data, config)

    prepared = captured["input_df"]
    assert prepared["model_name"].tolist() == ["single", "single"]
    assert prepared["mol_idx"].astype(str).str.startswith("LP-").all()


def test_stage_preserves_existing_mol_idx_values(tmp_path, monkeypatch):
    config = _build_config(tmp_path)
    data = pd.DataFrame(
        {
            "smiles": ["CCO", "CCC"],
            "model_name": ["m1", "m1"],
            "mol_idx": ["keep-1", "keep-2"],
        }
    )
    captured = {}

    def _fake_compute_metrics(df, *args, **kwargs):
        captured["input_df"] = df.copy()
        return df.copy()

    monkeypatch.setattr(descriptors_stage, "compute_metrics", _fake_compute_metrics)
    monkeypatch.setattr(descriptors_stage, "run_waves", lambda *_args, **_kwargs: None)

    descriptors_stage.run(data, config)

    prepared = captured["input_df"]
    assert prepared["mol_idx"].tolist() == ["keep-1", "keep-2"]
