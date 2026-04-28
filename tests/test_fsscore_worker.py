"""Tests for FSScore isolated worker wrapper."""

from __future__ import annotations

import argparse
from types import SimpleNamespace

from hedgehog.workers import fsscore_worker


def _args() -> argparse.Namespace:
    return argparse.Namespace(
        worker_python="/tmp/fsscore/bin/python",
        model_path="/tmp/model.ckpt",
        input_csv="/tmp/input.csv",
        output_csv="/tmp/output.csv",
        smiles_column="smiles",
        batch_size=128,
        num_workers=None,
        graph_datapath=None,
    )


def test_main_sets_cpu_default_for_subprocess(monkeypatch) -> None:
    """Worker should disable CUDA by default for broader host compatibility."""
    monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
    monkeypatch.setattr(fsscore_worker, "_parse_args", _args)

    captured: dict[str, object] = {}

    def _fake_run(command, **kwargs):  # noqa: ANN001
        captured["command"] = command
        captured["env"] = kwargs.get("env")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(fsscore_worker.subprocess, "run", _fake_run)

    assert fsscore_worker.main() == 0
    assert captured["command"][0:3] == [
        "/tmp/fsscore/bin/python",
        "-m",
        "fsscore.score",
    ]
    env = captured["env"]
    assert isinstance(env, dict)
    assert env.get("CUDA_VISIBLE_DEVICES") == ""


def test_main_preserves_explicit_cuda_env(monkeypatch) -> None:
    """Explicit CUDA_VISIBLE_DEVICES should be preserved, not overridden."""
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "0")
    monkeypatch.setattr(fsscore_worker, "_parse_args", _args)

    captured: dict[str, object] = {}

    def _fake_run(_command, **kwargs):  # noqa: ANN001
        captured["env"] = kwargs.get("env")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(fsscore_worker.subprocess, "run", _fake_run)

    assert fsscore_worker.main() == 0
    env = captured["env"]
    assert isinstance(env, dict)
    assert env.get("CUDA_VISIBLE_DEVICES") == "0"
