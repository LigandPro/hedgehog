"""Tests for FSScore setup helpers."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from hedgehog.setup._fsscore import (
    ensure_fsscore_checkout,
    ensure_fsscore_runtime,
)


def _write_min_checkout(repo_path: Path) -> None:
    (repo_path / "src" / "fsscore").mkdir(parents=True, exist_ok=True)
    (repo_path / "src" / "fsscore" / "score.py").write_text("", encoding="utf-8")
    (repo_path / "models").mkdir(parents=True, exist_ok=True)
    (
        repo_path / "models" / "pretrain_graph_GGLGGL_ep242_best_valloss.ckpt"
    ).write_bytes(b"x")


class TestEnsureFsscoreCheckout:
    """Tests for ensure_fsscore_checkout()."""

    def test_uses_existing_checkout(self, tmp_path: Path) -> None:
        repo_path = tmp_path / "modules" / "fsscore"
        _write_min_checkout(repo_path)

        result = ensure_fsscore_checkout(tmp_path)
        assert result == repo_path

    def test_clones_when_missing(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr("hedgehog.setup._fsscore.confirm_download", lambda *_: True)

        calls: list[list[str]] = []

        def _fake_run(command, **kwargs):
            calls.append(command)
            repo_path = Path(command[-1])
            _write_min_checkout(repo_path)
            return SimpleNamespace(returncode=0)

        monkeypatch.setattr("hedgehog.setup._fsscore.subprocess.run", _fake_run)

        result = ensure_fsscore_checkout(tmp_path)
        assert result == tmp_path / "modules" / "fsscore"
        assert calls, "git clone should be called"
        assert calls[0][0:2] == ["git", "clone"]

    def test_raises_when_declined(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr(
            "hedgehog.setup._fsscore.confirm_download", lambda *_: False
        )

        with pytest.raises(RuntimeError, match="declined"):
            ensure_fsscore_checkout(tmp_path)

    def test_raises_when_existing_checkout_is_invalid(self, tmp_path: Path) -> None:
        invalid_repo = tmp_path / "modules" / "fsscore"
        invalid_repo.mkdir(parents=True, exist_ok=True)
        (invalid_repo / "README.md").write_text("broken checkout", encoding="utf-8")

        with pytest.raises(RuntimeError, match="missing required files"):
            ensure_fsscore_checkout(tmp_path)


class TestEnsureFsscoreRuntime:
    """Tests for ensure_fsscore_runtime()."""

    def test_reuses_existing_runtime(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        checkout_path = tmp_path / "modules" / "fsscore"
        _write_min_checkout(checkout_path)
        env_path = tmp_path / ".venv-fsscore-worker"
        worker_python = env_path / "bin" / "python"
        worker_python.parent.mkdir(parents=True, exist_ok=True)
        worker_python.write_text("#!/usr/bin/env python\n", encoding="utf-8")

        monkeypatch.setattr(
            "hedgehog.setup._fsscore._check_fsscore_runtime",
            lambda *_: True,
        )

        install_called = {"value": False}

        def _fail_install(*_args, **_kwargs):
            install_called["value"] = True
            raise AssertionError("runtime install should not run")

        monkeypatch.setattr(
            "hedgehog.setup._fsscore._install_fsscore_runtime", _fail_install
        )

        runtime = ensure_fsscore_runtime(tmp_path)
        assert runtime.checkout_path == checkout_path
        assert (
            runtime.model_path
            == checkout_path
            / "models"
            / "pretrain_graph_GGLGGL_ep242_best_valloss.ckpt"
        )
        assert runtime.worker_python == worker_python
        assert runtime.env_path == env_path
        assert install_called["value"] is False

    def test_prefers_optional_env_root(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        checkout_path = tmp_path / "modules" / "fsscore"
        _write_min_checkout(checkout_path)
        opt_root = tmp_path / "optional-envs"
        worker_python = opt_root / "fsscore" / "bin" / "python"
        worker_python.parent.mkdir(parents=True, exist_ok=True)
        worker_python.write_text("#!/usr/bin/env python\n", encoding="utf-8")

        monkeypatch.setenv("HEDGEHOG_OPTIONAL_ENV_ROOT", str(opt_root))
        monkeypatch.setattr(
            "hedgehog.setup._fsscore._check_fsscore_runtime",
            lambda *_: True,
        )

        runtime = ensure_fsscore_runtime(tmp_path)
        assert runtime.env_path == opt_root / "fsscore"
        assert runtime.worker_python == worker_python

    def test_installs_when_runtime_probe_fails(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        checkout_path = tmp_path / "modules" / "fsscore"
        _write_min_checkout(checkout_path)
        env_path = tmp_path / ".venv-fsscore-worker"
        worker_python = env_path / "bin" / "python"

        probes = iter([False, True])
        monkeypatch.setattr(
            "hedgehog.setup._fsscore._check_fsscore_runtime",
            lambda *_: next(probes),
        )
        monkeypatch.setattr("hedgehog.setup._fsscore.confirm_download", lambda *_: True)

        calls: list[tuple[Path, Path]] = []

        def _fake_install(repo_path: Path, target_env_path: Path) -> Path:
            calls.append((repo_path, target_env_path))
            worker_python.parent.mkdir(parents=True, exist_ok=True)
            worker_python.write_text("#!/usr/bin/env python\n", encoding="utf-8")
            return worker_python

        monkeypatch.setattr(
            "hedgehog.setup._fsscore._install_fsscore_runtime", _fake_install
        )

        runtime = ensure_fsscore_runtime(tmp_path)
        assert calls == [(checkout_path, env_path)]
        assert runtime.worker_python == worker_python

    def test_raises_when_install_declined(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        checkout_path = tmp_path / "modules" / "fsscore"
        _write_min_checkout(checkout_path)

        monkeypatch.setattr(
            "hedgehog.setup._fsscore._check_fsscore_runtime",
            lambda *_: False,
        )
        monkeypatch.setattr(
            "hedgehog.setup._fsscore.confirm_download", lambda *_: False
        )

        with pytest.raises(RuntimeError, match="declined"):
            ensure_fsscore_runtime(tmp_path)


def test_install_runtime_commands_include_setuptools_pin(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Runtime install should create venv and pin setuptools<81."""
    checkout_path = tmp_path / "modules" / "fsscore"
    checkout_path.mkdir(parents=True, exist_ok=True)
    env_path = tmp_path / ".venv-fsscore-worker"

    calls: list[list[str]] = []

    def _fake_run(command, **kwargs):
        calls.append(command)
        if command[0:2] == ["uv", "venv"]:
            worker_python = Path(command[-1]) / "bin" / "python"
            worker_python.parent.mkdir(parents=True, exist_ok=True)
            worker_python.write_text("#!/usr/bin/env python\n", encoding="utf-8")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr("hedgehog.setup._fsscore.subprocess.run", _fake_run)
    monkeypatch.setattr("hedgehog.setup._fsscore.resolve_uv_binary", lambda: "uv")

    from hedgehog.setup._fsscore import _install_fsscore_runtime

    worker_python = _install_fsscore_runtime(checkout_path, env_path)
    assert worker_python == env_path / "bin" / "python"
    assert calls[0][0:3] == ["uv", "venv", "--python"]
    assert calls[1][0:4] == ["uv", "pip", "install", "--python"]
    assert "-e" in calls[1]
    assert calls[2][-1] == "setuptools<81"
