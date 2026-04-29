"""Tests for optional GASA setup helper."""

from __future__ import annotations

from pathlib import Path

import pytest

from hedgehog.setup._gasa import _resolve_python_binary, ensure_gasa_worker


def _write_min_checkout(repo_path: Path) -> None:
    (repo_path / "model").mkdir(parents=True, exist_ok=True)
    (repo_path / "gasa.py").write_text("# stub", encoding="utf-8")
    (repo_path / "model" / "gasa.pth").write_bytes(b"x")
    (repo_path / "model" / "gasa.json").write_text("{}", encoding="utf-8")


def _write_venv_python(venv_dir: Path) -> Path:
    python_path = venv_dir / "bin" / "python"
    python_path.parent.mkdir(parents=True, exist_ok=True)
    python_path.write_text("", encoding="utf-8")
    return python_path


def test_resolve_python_binary_prefers_py310(monkeypatch: pytest.MonkeyPatch) -> None:
    """Auto-resolve should prefer Python 3.10 for GASA dependency compatibility."""
    candidates = {
        "python3.10": "/usr/bin/python3.10",
        "python3.11": "/usr/bin/python3.11",
        "python3.12": "/usr/bin/python3.12",
    }
    monkeypatch.setattr(
        "hedgehog.setup._gasa.shutil.which",
        lambda name: candidates.get(name),
    )

    assert _resolve_python_binary(None) == "/usr/bin/python3.10"


def test_resolve_python_binary_falls_back_to_uv_managed_py310(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """When system python3.10/3.11/3.12 are unavailable, return uv-managed target."""
    monkeypatch.setattr("hedgehog.setup._gasa.shutil.which", lambda _name: None)
    assert _resolve_python_binary(None) == "3.10"
    assert _resolve_python_binary("3.10") == "3.10"


def test_ensure_gasa_worker_clones_and_installs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr("hedgehog.setup._gasa.confirm_download", lambda *_: True)
    monkeypatch.setattr("hedgehog.setup._gasa.resolve_uv_binary", lambda: "uv")
    monkeypatch.setattr(
        "hedgehog.setup._gasa._resolve_python_binary", lambda *_: "python3.10"
    )

    calls: list[tuple[list[str], int]] = []

    def _fake_run(cmd: list[str], cwd: Path, timeout: int = 1800) -> None:
        calls.append((cmd, timeout))
        if cmd[:2] == ["git", "clone"]:
            _write_min_checkout(Path(cmd[-1]))
        elif cmd[:2] == ["uv", "venv"]:
            _write_venv_python(Path(cmd[-1]))

    monkeypatch.setattr("hedgehog.setup._gasa._run", _fake_run)

    result = ensure_gasa_worker(tmp_path)
    assert result.repo_path == tmp_path / "modules" / "gasa"
    assert result.worker_python == tmp_path / ".venv-gasa-worker" / "bin" / "python"
    assert any(cmd[:2] == ["git", "clone"] for cmd, _timeout in calls)
    assert any(
        cmd[:4] == ["uv", "pip", "install", "--python"] and "dgl==1.1.3" in cmd
        for cmd, _timeout in calls
    )
    assert any(
        cmd[:4] == ["uv", "pip", "install", "--python"] and "setuptools<81" in cmd
        for cmd, _timeout in calls
    )
    assert any(cmd[:3] == ["uv", "venv", "--python"] for cmd, _timeout in calls)
    assert any(
        cmd[1] == "-c"
        and "import numpy, pandas, sklearn, hyperopt, rdkit, torch, dgl, dgllife"
        in cmd[2]
        for cmd, _timeout in calls
    )
    assert any(cmd[1] == "-c" and timeout >= 600 for cmd, timeout in calls)


def test_ensure_gasa_worker_honors_optional_env_root(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    repo_path = tmp_path / "modules" / "gasa"
    _write_min_checkout(repo_path)

    custom_root = tmp_path / "shared" / "optional-envs"
    worker_python = _write_venv_python(custom_root / "gasa")
    monkeypatch.setenv("HEDGEHOG_OPTIONAL_ENV_ROOT", str(custom_root))

    monkeypatch.setattr("hedgehog.setup._gasa._run", lambda *args, **kwargs: None)

    result = ensure_gasa_worker(tmp_path)
    assert result.repo_path == repo_path
    assert result.worker_python == worker_python


def test_ensure_gasa_worker_recreates_existing_env_when_stale(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    repo_path = tmp_path / "modules" / "gasa"
    _write_min_checkout(repo_path)
    _write_venv_python(tmp_path / ".venv-gasa-worker")

    monkeypatch.setattr("hedgehog.setup._gasa.confirm_download", lambda *_: True)
    monkeypatch.setattr("hedgehog.setup._gasa.resolve_uv_binary", lambda: "uv")
    monkeypatch.setattr(
        "hedgehog.setup._gasa._resolve_python_binary", lambda *_: "python3.10"
    )

    verify_calls = {"count": 0}

    def _fake_verify(_python_path: Path) -> None:
        verify_calls["count"] += 1
        if verify_calls["count"] == 1:
            raise RuntimeError("stale env")

    monkeypatch.setattr("hedgehog.setup._gasa._verify_gasa_dependencies", _fake_verify)

    commands: list[list[str]] = []
    removed_paths: list[Path] = []

    def _fake_rmtree(path: Path) -> None:
        removed_paths.append(path)
        python_path = path / "bin" / "python"
        if python_path.exists():
            python_path.unlink()

    def _fake_run(cmd: list[str], cwd: Path, timeout: int = 1800) -> None:
        del cwd, timeout
        commands.append(cmd)
        if cmd[:2] == ["uv", "venv"]:
            _write_venv_python(Path(cmd[-1]))

    monkeypatch.setattr("hedgehog.setup._gasa.shutil.rmtree", _fake_rmtree)
    monkeypatch.setattr("hedgehog.setup._gasa._run", _fake_run)

    result = ensure_gasa_worker(tmp_path)
    assert result.worker_python == tmp_path / ".venv-gasa-worker" / "bin" / "python"
    assert removed_paths == [tmp_path / ".venv-gasa-worker"]
    assert any(
        cmd[:3] == ["uv", "venv", "--python"] and "--clear" not in cmd
        for cmd in commands
    )


def test_ensure_gasa_worker_raises_for_invalid_existing_checkout(
    tmp_path: Path,
) -> None:
    bad_repo = tmp_path / "modules" / "gasa"
    bad_repo.mkdir(parents=True, exist_ok=True)
    (bad_repo / "README.md").write_text("invalid", encoding="utf-8")

    with pytest.raises(RuntimeError, match="missing required files"):
        ensure_gasa_worker(tmp_path)
