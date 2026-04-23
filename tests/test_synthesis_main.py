"""Tests for synthesis/main.py path resolution and setup fallback."""

from __future__ import annotations

import importlib
from pathlib import Path

import pytest

import hedgehog.synthesis.main as synthesis_main


def test_resolve_retrosynthesis_config_expands_tilde(monkeypatch):
    monkeypatch.setenv("HOME", "/tmp/fake-home")
    result = synthesis_main._resolve_retrosynthesis_config(
        {"config_retrosynthesis": "~/retro/config.yml"}
    )
    expected = Path("/tmp/fake-home/retro/config.yml").resolve(strict=False)
    assert result == expected


def test_resolve_retrosynthesis_config_uses_project_root_for_relative():
    result = synthesis_main._resolve_retrosynthesis_config(
        {"config_retrosynthesis": "configs/retro.yml"}
    )
    expected_root = synthesis_main._project_root()
    assert result == (expected_root / "configs/retro.yml").resolve(strict=False)


def test_main_raises_for_unsupported_python_autoinstall(monkeypatch, tmp_path: Path):
    output_dir = tmp_path / "results"
    output_dir.mkdir(parents=True, exist_ok=True)
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("smiles,model_name,mol_idx\nCCO,m,0\n", encoding="utf-8")
    synth_cfg = tmp_path / "config_synthesis.yml"
    synth_cfg.write_text(
        "run: true\nrun_retrosynthesis: true\nfilter_solved_only: true\n",
        encoding="utf-8",
    )
    missing_retro_cfg = tmp_path / "missing-retro.yml"

    def _mock_load_config(path):
        if Path(path) == synth_cfg:
            return {"run_retrosynthesis": True, "filter_solved_only": True}
        return {}

    def _copy_df(input_df, *_args, **_kwargs):
        return input_df.copy()

    monkeypatch.setattr(synthesis_main, "load_config", _mock_load_config)
    monkeypatch.setattr(
        synthesis_main,
        "get_input_path",
        lambda config, folder_to_save: str(input_csv),
    )
    monkeypatch.setattr(synthesis_main, "calculate_synthesis_scores", _copy_df)
    monkeypatch.setattr(synthesis_main, "apply_synthesis_score_filters", _copy_df)
    setup_module = importlib.import_module("hedgehog.setup")
    unsupported_error = "AiZynthFinder upstream supports Python 3.10-3.12"
    monkeypatch.setattr(
        setup_module,
        "ensure_aizynthfinder",
        lambda project_root: (_ for _ in ()).throw(RuntimeError(unsupported_error)),
    )

    config = {
        "folder_to_save": str(output_dir),
        "config_synthesis": str(synth_cfg),
        "generated_mols_path": str(input_csv),
        "config_retrosynthesis": str(missing_retro_cfg),
    }

    with pytest.raises(RuntimeError, match="supports Python 3.10-3.12"):
        synthesis_main.main(config)
