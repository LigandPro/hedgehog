"""Tests for isolated GASA worker launcher."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

from hedgehog.workers import gasa_worker


def test_gasa_worker_writes_output_csv(tmp_path: Path, monkeypatch) -> None:
    repo_path = tmp_path / "modules" / "gasa"
    repo_path.mkdir(parents=True, exist_ok=True)

    input_csv = tmp_path / "input.csv"
    output_csv = tmp_path / "output.csv"
    pd.DataFrame({"smiles": ["CCO", "c1ccccc1"]}).to_csv(input_csv, index=False)

    def _fake_run(command, **kwargs):
        assert command[0] == "/tmp/gasa-worker-python"
        assert command[1] == "-c"
        assert kwargs["cwd"] == str(repo_path)
        assert kwargs["env"]["PYTHONPATH"].split(os.pathsep)[0] == str(repo_path)

        source = pd.read_csv(command[3])
        source["gasa_score"] = [0.2, 0.8]
        source.to_csv(command[4], index=False)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(gasa_worker.subprocess, "run", _fake_run)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "gasa_worker.py",
            "--worker-python",
            "/tmp/gasa-worker-python",
            "--repo-path",
            str(repo_path),
            "--input-csv",
            str(input_csv),
            "--output-csv",
            str(output_csv),
            "--smiles-column",
            "smiles",
        ],
    )

    exit_code = gasa_worker.main()
    assert exit_code == 0

    out_df = pd.read_csv(output_csv)
    assert out_df["gasa_score"].tolist() == [0.2, 0.8]


def test_gasa_worker_fails_when_repo_is_missing(tmp_path: Path, monkeypatch) -> None:
    input_csv = tmp_path / "input.csv"
    pd.DataFrame({"smiles": ["CCO"]}).to_csv(input_csv, index=False)

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "gasa_worker.py",
            "--worker-python",
            "/tmp/gasa-worker-python",
            "--repo-path",
            str(tmp_path / "missing-repo"),
            "--input-csv",
            str(input_csv),
            "--output-csv",
            str(tmp_path / "output.csv"),
        ],
    )

    assert gasa_worker.main() == 1
