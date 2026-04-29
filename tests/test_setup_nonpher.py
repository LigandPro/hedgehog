"""Tests for optional Nonpher setup/check helpers."""

from __future__ import annotations

import csv
import subprocess

import pytest

from hedgehog.setup import _nonpher as nonpher_setup
from hedgehog.workers import nonpher_worker


def test_check_nonpher_runtime_uses_external_python(monkeypatch):
    """When --python is provided, helper should probe that interpreter."""

    def _fake_external(python_bin: str, probe_smiles: str):
        return nonpher_setup.NonpherCheckResult(
            available=True,
            detail=f"{python_bin}:{probe_smiles}",
        )

    monkeypatch.setattr(nonpher_setup, "_check_nonpher_external_python", _fake_external)
    monkeypatch.setattr(
        nonpher_setup,
        "_check_nonpher_current_env",
        lambda probe_smiles: nonpher_setup.NonpherCheckResult(False, probe_smiles),
    )

    result = nonpher_setup.check_nonpher_runtime(
        python_bin="/tmp/nonpher/bin/python",
        probe_smiles="CCO",
    )

    assert result.available is True
    assert result.detail == "/tmp/nonpher/bin/python:CCO"


def test_check_nonpher_runtime_uses_env_python(monkeypatch):
    """When env var is set, helper should probe that interpreter."""

    monkeypatch.setenv(nonpher_setup.NONPHER_PYTHON_ENV_VAR, "/tmp/env/nonpher-python")
    monkeypatch.setattr(
        nonpher_setup,
        "_check_nonpher_external_python",
        lambda python_bin, probe_smiles: nonpher_setup.NonpherCheckResult(
            available=True,
            detail=f"{python_bin}:{probe_smiles}",
        ),
    )

    result = nonpher_setup.check_nonpher_runtime(probe_smiles="CCC")
    assert result.available is True
    assert result.detail == "/tmp/env/nonpher-python:CCC"


def test_check_nonpher_runtime_reports_missing_dependency(monkeypatch):
    """Missing dependency in current env should return unavailable result."""

    def _raise_missing():
        raise ModuleNotFoundError("No module named 'nonpher'")

    monkeypatch.setattr(
        nonpher_setup, "create_nonpher_complexity_filter", _raise_missing
    )

    result = nonpher_setup.check_nonpher_runtime()
    assert result.available is False
    assert "ModuleNotFoundError" in result.detail


def test_check_nonpher_runtime_external_stderr(monkeypatch):
    """External probe should surface first stderr line in failure detail."""

    def _fake_run(*args, **kwargs):
        del kwargs
        return subprocess.CompletedProcess(
            args=args[0],
            returncode=1,
            stdout="",
            stderr="No module named nonpher\nextra line",
        )

    monkeypatch.setattr(nonpher_setup.subprocess, "run", _fake_run)

    result = nonpher_setup.check_nonpher_runtime(
        python_bin="/tmp/nonpher/bin/python",
        probe_smiles="CCO",
    )
    assert result.available is False
    assert "No module named nonpher" in result.detail


def test_check_nonpher_runtime_external_probe_uses_newlines(monkeypatch):
    """External probe script should pass real newlines, not escaped sequences."""

    captured: dict[str, str] = {}

    def _fake_run(command, **kwargs):
        captured["script"] = command[2]
        del kwargs
        return subprocess.CompletedProcess(
            args=command,
            returncode=0,
            stdout="ok\n",
            stderr="",
        )

    monkeypatch.setattr(nonpher_setup.subprocess, "run", _fake_run)

    result = nonpher_setup.check_nonpher_runtime(
        python_bin="/tmp/nonpher/bin/python",
        probe_smiles="CCO",
    )
    assert result.available is True
    assert "\\n" not in captured["script"]
    assert "\n" in captured["script"]


def test_ensure_nonpher_external_runtime_uses_existing_python(monkeypatch):
    """ensure helper should reuse configured external interpreter when probe succeeds."""

    monkeypatch.setenv(
        nonpher_setup.NONPHER_PYTHON_ENV_VAR, "/tmp/nonpher-existing/bin/python"
    )
    monkeypatch.setattr(
        nonpher_setup,
        "_check_nonpher_external_python",
        lambda python_bin, probe_smiles: nonpher_setup.NonpherCheckResult(
            available=True,
            detail=f"ok:{python_bin}:{probe_smiles}",
        ),
    )

    result = nonpher_setup.ensure_nonpher_external_runtime(probe_smiles="CCN")
    assert result.available is True
    assert result.python_bin == "/tmp/nonpher-existing/bin/python"
    assert result.install_attempted is False


def test_ensure_nonpher_external_runtime_reports_conda_blocker(monkeypatch):
    """ensure helper should expose solver blocker when conda create fails."""

    def _fake_run(command: list[str], timeout: int = 3600):
        del timeout
        return subprocess.CompletedProcess(
            args=command,
            returncode=1,
            stdout="",
            stderr="LibMambaUnsatisfiableError: nothing provides libboost 1.65.*",
        )

    monkeypatch.setattr(nonpher_setup, "_run_command", _fake_run)

    result = nonpher_setup.ensure_nonpher_external_runtime(
        env_prefix="/tmp/nonpher-test-env"
    )
    assert result.available is False
    assert result.install_attempted is True
    assert "LibMambaUnsatisfiableError" in result.detail


def test_ensure_nonpher_uv_runtime_uses_existing_python(monkeypatch, tmp_path):
    """uv helper should return immediately when existing prefix python passes probe."""

    env_prefix = tmp_path / "nonpher-uv"
    python_path = env_prefix / "bin" / "python"
    python_path.parent.mkdir(parents=True, exist_ok=True)
    python_path.write_text("#!/usr/bin/env python\n", encoding="utf-8")
    (env_prefix / "pyvenv.cfg").write_text("home = /tmp\n", encoding="utf-8")

    monkeypatch.setattr(
        nonpher_setup,
        "_check_nonpher_external_python",
        lambda *_: nonpher_setup.NonpherCheckResult(True, "ok"),
    )
    monkeypatch.setattr(
        nonpher_setup,
        "_run_command",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("uv bootstrap commands must not run")
        ),
    )

    result = nonpher_setup.ensure_nonpher_uv_runtime(env_prefix=str(env_prefix))
    assert result.available is True
    assert result.install_attempted is False
    assert result.python_bin == str(python_path)


def test_ensure_nonpher_uv_runtime_reports_linker_blocker(monkeypatch, tmp_path):
    """uv helper should surface native linker blocker for molpher build."""

    env_prefix = tmp_path / "nonpher-uv"
    python_path = env_prefix / "bin" / "python"

    monkeypatch.setattr(nonpher_setup, "resolve_uv_binary", lambda: "uv")
    monkeypatch.setattr(
        nonpher_setup,
        "_check_nonpher_external_python",
        lambda *_: nonpher_setup.NonpherCheckResult(False, "missing"),
    )

    steps = iter(
        [
            subprocess.CompletedProcess(
                args=["uv", "venv"], returncode=0, stdout="", stderr=""
            ),
            subprocess.CompletedProcess(
                args=["uv", "pip", "install", "numpy<2"],
                returncode=0,
                stdout="",
                stderr="",
            ),
            subprocess.CompletedProcess(
                args=["uv", "pip", "install", "rdkit-pypi==2022.9.5"],
                returncode=0,
                stdout="",
                stderr="",
            ),
            subprocess.CompletedProcess(
                args=["uv", "pip", "install", "setuptools<81", "wheel"],
                returncode=0,
                stdout="",
                stderr="",
            ),
            subprocess.CompletedProcess(
                args=["uv", "pip", "install", "--no-deps", "nonpher"],
                returncode=0,
                stdout="",
                stderr="",
            ),
            subprocess.CompletedProcess(
                args=["uv", "pip", "install", "--no-build-isolation", "molpher-lib"],
                returncode=1,
                stdout="",
                stderr="... /usr/bin/ld: cannot find -lmolpher ...",
            ),
        ]
    )
    monkeypatch.setattr(
        nonpher_setup, "_run_command", lambda *_args, **_kwargs: next(steps)
    )

    result = nonpher_setup.ensure_nonpher_uv_runtime(env_prefix=str(env_prefix))
    assert result.available is False
    assert result.install_attempted is True
    assert result.python_bin == str(python_path)
    assert "cannot find -lmolpher" in result.detail


def test_ensure_nonpher_uv_runtime_uses_fallback_after_linker_blocker(
    monkeypatch, tmp_path
):
    """uv linker failure should not force NaN when a validated fallback exists."""

    env_prefix = tmp_path / "nonpher-uv"
    fallback_python = tmp_path / "nonpher-hybrid" / "bin" / "python"
    fallback_python.parent.mkdir(parents=True)
    fallback_python.write_text("#!/usr/bin/env python\n", encoding="utf-8")

    monkeypatch.setattr(nonpher_setup, "resolve_uv_binary", lambda: "uv")

    def _fake_check(python_bin: str, _probe_smiles: str):
        return nonpher_setup.NonpherCheckResult(
            available=python_bin == str(fallback_python),
            detail="fallback ok" if python_bin == str(fallback_python) else "missing",
        )

    monkeypatch.setattr(nonpher_setup, "_check_nonpher_external_python", _fake_check)

    steps = iter(
        [
            subprocess.CompletedProcess(
                args=["uv", "venv"], returncode=0, stdout="", stderr=""
            ),
            subprocess.CompletedProcess(
                args=["uv", "pip", "install", "numpy<2"],
                returncode=0,
                stdout="",
                stderr="",
            ),
            subprocess.CompletedProcess(
                args=["uv", "pip", "install", "rdkit-pypi==2022.9.5"],
                returncode=0,
                stdout="",
                stderr="",
            ),
            subprocess.CompletedProcess(
                args=["uv", "pip", "install", "setuptools<81", "wheel"],
                returncode=0,
                stdout="",
                stderr="",
            ),
            subprocess.CompletedProcess(
                args=["uv", "pip", "install", "--no-deps", "nonpher"],
                returncode=0,
                stdout="",
                stderr="",
            ),
            subprocess.CompletedProcess(
                args=["uv", "pip", "install", "--no-build-isolation", "molpher-lib"],
                returncode=1,
                stdout="",
                stderr="... /usr/bin/ld: cannot find -lmolpher ...",
            ),
        ]
    )
    monkeypatch.setattr(
        nonpher_setup, "_run_command", lambda *_args, **_kwargs: next(steps)
    )

    result = nonpher_setup.ensure_nonpher_uv_runtime(
        env_prefix=str(env_prefix),
        fallback_python=str(fallback_python),
    )

    assert result.available is True
    assert result.install_attempted is True
    assert result.python_bin == str(fallback_python)
    assert "fallback" in result.detail


def test_ensure_nonpher_uv_runtime_preserves_non_venv_prefix(monkeypatch, tmp_path):
    """uv helper must not delete arbitrary directories passed as env_prefix."""

    env_prefix = tmp_path / "not-a-venv"
    env_prefix.mkdir()
    keep_file = env_prefix / "keep.txt"
    keep_file.write_text("keep", encoding="utf-8")

    monkeypatch.setattr(
        nonpher_setup,
        "resolve_uv_binary",
        lambda: (_ for _ in ()).throw(AssertionError("uv should not be resolved")),
    )

    result = nonpher_setup.ensure_nonpher_uv_runtime(env_prefix=str(env_prefix))

    assert result.available is False
    assert result.install_attempted is False
    assert "Refusing to overwrite non-virtualenv path" in result.detail
    assert keep_file.read_text(encoding="utf-8") == "keep"


def test_nonpher_lobachevsky_setup_commands_include_isolated_env():
    """Guidance should prefer uv-only setup and include external fallback."""
    commands = nonpher_setup.nonpher_lobachevsky_setup_commands(
        optional_env_root="~/work/nonpher-opt-envs",
        fallback_python="/tmp/nonpher-hybrid/bin/python",
    )

    assert any(
        "export HEDGEHOG_OPTIONAL_ENV_ROOT=~/work/nonpher-opt-envs" in command
        for command in commands
    )
    assert any("ensure_nonpher_uv_runtime" in command for command in commands)
    assert any(
        f"export {nonpher_setup.NONPHER_PYTHON_ENV_VAR}=/tmp/nonpher-hybrid/bin/python"
        in command
        for command in commands
    )
    assert commands[-1].endswith("--python /tmp/nonpher-hybrid/bin/python")


def test_nonpher_batch_worker_command_has_expected_shape():
    """Command builder should target external python + worker script."""

    command = nonpher_setup.nonpher_batch_worker_command(
        worker_python="/tmp/nonpher/bin/python",
        input_csv="/tmp/in.csv",
        output_csv="/tmp/out.csv",
    )
    assert command[0] == "/tmp/nonpher/bin/python"
    assert command[1].endswith("/workers/nonpher_worker.py")
    assert "--input-csv" in command
    assert "--output-csv" in command


def test_run_nonpher_batch_external_uses_env_python(monkeypatch):
    """Batch launcher should use HEDGEHOG_NONPHER_PYTHON when python_bin omitted."""

    captured: dict[str, list[str]] = {}

    def _fake_run(command: list[str], timeout: int = 3600):
        del timeout
        captured["command"] = command
        return subprocess.CompletedProcess(
            args=command, returncode=0, stdout="", stderr=""
        )

    monkeypatch.setenv(nonpher_setup.NONPHER_PYTHON_ENV_VAR, "/tmp/nonpher/env/python")
    monkeypatch.setattr(nonpher_setup, "_run_command", _fake_run)

    nonpher_setup.run_nonpher_batch_external(
        input_csv="/tmp/in.csv",
        output_csv="/tmp/out.csv",
    )
    assert captured["command"][0] == "/tmp/nonpher/env/python"


def test_run_nonpher_batch_external_surfaces_worker_failure(monkeypatch):
    """Batch launcher should include worker stderr first line in raised error."""

    monkeypatch.setattr(
        nonpher_setup,
        "_run_command",
        lambda command, timeout=3600: subprocess.CompletedProcess(
            args=command,
            returncode=1,
            stdout="",
            stderr="worker failed\nextra diagnostics",
        ),
    )

    with pytest.raises(RuntimeError, match="worker failed"):
        nonpher_setup.run_nonpher_batch_external(
            input_csv="/tmp/in.csv",
            output_csv="/tmp/out.csv",
            python_bin="/tmp/nonpher/bin/python",
        )


def test_nonpher_worker_main_scores_rows_without_pandas(monkeypatch, tmp_path):
    """Worker should score rows via stdlib CSV and mocked Nonpher/RDKit hooks."""

    input_csv = tmp_path / "input.csv"
    output_csv = tmp_path / "output.csv"
    input_csv.write_text("smiles,id\nCCO,1\nCCCCC,2\nbad,3\n,4\n", encoding="utf-8")

    class _FakeFilter:
        def isTooComplex(self, mol):
            return bool(mol == "complex")

    monkeypatch.setattr(
        nonpher_worker, "_create_complexity_filter", lambda: _FakeFilter()
    )
    monkeypatch.setattr(
        nonpher_worker,
        "_mol_from_smiles",
        lambda smiles: (
            None if smiles == "bad" else ("complex" if smiles == "CCCCC" else "simple")
        ),
    )

    code = nonpher_worker.main(
        [
            "--input-csv",
            str(input_csv),
            "--output-csv",
            str(output_csv),
            "--smiles-column",
            "smiles",
        ]
    )
    assert code == 0

    with output_csv.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))

    assert rows[0]["nonpher_complexity_score"] == "0.0"
    assert rows[1]["nonpher_complexity_score"] == "1.0"
    assert rows[2]["nonpher_complexity_score"] == ""
    assert rows[3]["nonpher_complexity_score"] == ""


def test_nonpher_worker_main_missing_smiles_column_returns_error(tmp_path):
    """Worker should fail fast when expected smiles column is missing."""

    input_csv = tmp_path / "input.csv"
    output_csv = tmp_path / "output.csv"
    input_csv.write_text("id,value\n1,a\n", encoding="utf-8")

    code = nonpher_worker.main(
        [
            "--input-csv",
            str(input_csv),
            "--output-csv",
            str(output_csv),
            "--smiles-column",
            "smiles",
        ]
    )
    assert code == 1
