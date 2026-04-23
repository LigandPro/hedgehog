"""Tests for hedgehog.setup._aizynthfinder auto-installation."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from hedgehog.setup._aizynthfinder import ensure_aizynthfinder

_AIZYNTHFINDER_PROBE = (
    "import importlib.util, sys; "
    "sys.exit(0 if importlib.util.find_spec('aizynthfinder') else 1)"
)


def _install_root(root: Path) -> Path:
    """Return the expected AiZynthFinder workspace root for a project root."""
    return root / "modules" / "aizynthfinder"


def _config_path(root: Path) -> Path:
    """Return the expected config.yml path for a given project root."""
    return _install_root(root) / "public" / "config.yml"


def _make_config(root: Path) -> Path:
    """Create the config.yml file so the function thinks it's installed."""
    cfg = _config_path(root)
    cfg.parent.mkdir(parents=True, exist_ok=True)
    cfg.write_text("version: 1\n", encoding="utf-8")
    return cfg


def _cmd_as_list(cmd: list[str] | tuple[str, ...]) -> list[str]:
    """Normalize subprocess command into a list of strings."""
    return [str(part) for part in cmd]


def _is_probe_command(cmd: list[str] | tuple[str, ...], uv_bin: str) -> bool:
    """Return True when command is the import probe for aizynthfinder."""
    cmd_list = _cmd_as_list(cmd)
    return cmd_list[:4] == [uv_bin, "run", "python", "-c"]


def _is_public_data_download(cmd: list[str] | tuple[str, ...]) -> bool:
    """Return True when command triggers public data download."""
    return any("download_public_data" in str(part) for part in cmd)


def _write_logging_yaml(root: Path) -> None:
    """Create synthesis/logging.yml source file for copy step."""
    logging_src = root / "src" / "hedgehog" / "synthesis" / "logging.yml"
    logging_src.parent.mkdir(parents=True, exist_ok=True)
    logging_src.write_text("version: 1\n", encoding="utf-8")


class TestEnsureAizynthfinder:
    """Tests for the ensure_aizynthfinder function."""

    @pytest.fixture(autouse=True)
    def _clear_uv_env(self, monkeypatch):
        monkeypatch.delenv("UV", raising=False)

    @pytest.fixture(autouse=True)
    def _supported_python(self, monkeypatch):
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.sys.version_info", (3, 12, 0)
        )

    def test_already_installed(self, tmp_path: Path, monkeypatch):
        """If config.yml and package exist, return immediately."""
        cfg = _make_config(tmp_path)
        _write_logging_yaml(tmp_path)
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.resolve_uv_binary",
            lambda: "/usr/bin/uv",
        )
        calls: list[list[str]] = []

        def _mock_run(cmd, *, cwd=None, check=False, timeout=None):
            del cwd, check, timeout
            calls.append([str(part) for part in cmd])
            return subprocess.CompletedProcess(cmd, 0)

        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.subprocess.run",
            _mock_run,
        )

        result = ensure_aizynthfinder(tmp_path)

        assert result == cfg
        assert len(calls) == 1
        assert calls[0][:4] == [
            "/usr/bin/uv",
            "run",
            "python",
            "-c",
        ]
        assert (
            _install_root(tmp_path) / "aizynthfinder" / "data" / "logging.yml"
        ).exists()

    def test_existing_config_reinstalls_missing_package(self, tmp_path: Path, monkeypatch):
        """Existing public data should still install the package if it is missing."""
        cfg = _make_config(tmp_path)
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.resolve_uv_binary",
            lambda: "/usr/bin/uv",
        )
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.confirm_download",
            lambda *_a, **_kw: True,
        )

        calls: list[list[str]] = []

        def _mock_run(cmd, *, cwd=None, check=False, timeout=None):
            del cwd, check, timeout
            cmd_list = _cmd_as_list(cmd)
            calls.append(cmd_list)
            if _is_probe_command(cmd, "/usr/bin/uv"):
                return subprocess.CompletedProcess(cmd, 1)
            return subprocess.CompletedProcess(cmd, 0)

        monkeypatch.setattr("hedgehog.setup._aizynthfinder.subprocess.run", _mock_run)

        result = ensure_aizynthfinder(tmp_path)

        assert result == cfg
        assert calls == [
            ["/usr/bin/uv", "run", "python", "-c", _AIZYNTHFINDER_PROBE],
            ["/usr/bin/uv", "sync", "--extra", "retrosynthesis"],
        ]

    def test_uv_not_found(self, tmp_path: Path, monkeypatch):
        """RuntimeError when uv is not on PATH."""
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.resolve_uv_binary",
            lambda: (_ for _ in ()).throw(RuntimeError("uv is not installed")),
        )

        with pytest.raises(RuntimeError, match="uv is not installed"):
            ensure_aizynthfinder(tmp_path)

    def test_user_declines(self, tmp_path: Path, monkeypatch):
        """RuntimeError when user declines the download prompt."""
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.resolve_uv_binary",
            lambda: "/usr/bin/uv",
        )
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.confirm_download",
            lambda *_a, **_kw: False,
        )

        with pytest.raises(RuntimeError, match="declined"):
            ensure_aizynthfinder(tmp_path)

    def test_unsupported_python(self, tmp_path: Path, monkeypatch):
        """Unsupported interpreter versions should fail before any setup runs."""
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.sys.version_info", (3, 13, 0)
        )

        with pytest.raises(RuntimeError, match="supports Python 3.10-3.12"):
            ensure_aizynthfinder(tmp_path)

    def test_uses_uv_from_environment(self, tmp_path: Path, monkeypatch):
        """Resolved uv binary should be used for probe, sync, and data download."""
        uv_bin = tmp_path / "custom-bin" / "uv"
        uv_bin.parent.mkdir(parents=True, exist_ok=True)
        uv_bin.write_text("#!/bin/sh\n", encoding="utf-8")
        uv_bin.chmod(0o755)
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.resolve_uv_binary",
            lambda: str(uv_bin),
        )
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.confirm_download",
            lambda *_a, **_kw: True,
        )

        calls: list[list[str]] = []

        def _mock_run(cmd, *, cwd=None, check=False, timeout=None):
            del cwd, check, timeout
            calls.append(_cmd_as_list(cmd))
            if _is_probe_command(cmd, str(uv_bin)):
                return subprocess.CompletedProcess(cmd, 1)
            if _is_public_data_download(cmd):
                _make_config(tmp_path)
            return subprocess.CompletedProcess(cmd, 0)

        monkeypatch.setattr("hedgehog.setup._aizynthfinder.subprocess.run", _mock_run)

        ensure_aizynthfinder(tmp_path)

        assert calls[0][:2] == [str(uv_bin), "run"]
        assert calls[1] == [str(uv_bin), "sync", "--extra", "retrosynthesis"]
        assert calls[2][0] == str(uv_bin)

    def test_full_install_cycle(self, tmp_path: Path, monkeypatch):
        """Full flow: probe, sync, download, copy logging, return config."""
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.resolve_uv_binary",
            lambda: "/usr/bin/uv",
        )
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.confirm_download",
            lambda *_a, **_kw: True,
        )

        _write_logging_yaml(tmp_path)

        subprocess_calls: list[tuple[list[str], Path | None]] = []

        def _mock_run(cmd, *, cwd=None, check=False, timeout=None):
            del check, timeout
            cmd_list = _cmd_as_list(cmd)
            subprocess_calls.append((cmd_list, cwd))
            if _is_probe_command(cmd, "/usr/bin/uv"):
                return subprocess.CompletedProcess(cmd, 1)
            if _is_public_data_download(cmd):
                _make_config(tmp_path)
            return subprocess.CompletedProcess(cmd, 0)

        monkeypatch.setattr("hedgehog.setup._aizynthfinder.subprocess.run", _mock_run)

        result = ensure_aizynthfinder(tmp_path)

        assert result == _config_path(tmp_path)
        cmds = [cmd for cmd, _cwd in subprocess_calls]
        assert len(cmds) == 3
        assert cmds[0][:4] == ["/usr/bin/uv", "run", "python", "-c"]
        assert cmds[1] == ["/usr/bin/uv", "sync", "--extra", "retrosynthesis"]
        assert "aizynthfinder.tools.download_public_data" in cmds[2]
        assert (
            _install_root(tmp_path) / "aizynthfinder" / "data" / "logging.yml"
        ).exists()

    def test_skips_sync_if_package_is_available(self, tmp_path: Path, monkeypatch):
        """uv sync is skipped when aizynthfinder is already importable."""
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.resolve_uv_binary",
            lambda: "/usr/bin/uv",
        )
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.confirm_download",
            lambda *_a, **_kw: True,
        )

        subprocess_calls: list[list[str]] = []

        def _mock_run(cmd, *, cwd=None, check=False, timeout=None):
            del cwd, check, timeout
            cmd_list = _cmd_as_list(cmd)
            subprocess_calls.append(cmd_list)
            if _is_probe_command(cmd, "/usr/bin/uv"):
                return subprocess.CompletedProcess(cmd, 0)
            if _is_public_data_download(cmd):
                _make_config(tmp_path)
            return subprocess.CompletedProcess(cmd, 0)

        monkeypatch.setattr("hedgehog.setup._aizynthfinder.subprocess.run", _mock_run)

        ensure_aizynthfinder(tmp_path)

        assert not any(cmd[:2] == ["/usr/bin/uv", "sync"] for cmd in subprocess_calls)

    def test_skips_download_if_public_has_files(self, tmp_path: Path, monkeypatch):
        """Existing config should skip download even though the package is probed."""
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.resolve_uv_binary",
            lambda: "/usr/bin/uv",
        )
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.confirm_download",
            lambda *_a, **_kw: True,
        )

        public = _install_root(tmp_path) / "public"
        public.mkdir(parents=True, exist_ok=True)
        (public / "config.yml").write_text("version: 1\n", encoding="utf-8")
        (public / "model.onnx").write_text("fake model", encoding="utf-8")

        subprocess_calls: list[list[str]] = []

        def _mock_run(cmd, *, cwd=None, check=False, timeout=None):
            del cwd, check, timeout
            subprocess_calls.append(_cmd_as_list(cmd))
            return subprocess.CompletedProcess(cmd, 0)

        monkeypatch.setattr("hedgehog.setup._aizynthfinder.subprocess.run", _mock_run)

        result = ensure_aizynthfinder(tmp_path)

        assert result == _config_path(tmp_path)
        assert subprocess_calls == [
            [
                "/usr/bin/uv",
                "run",
                "python",
                "-c",
                _AIZYNTHFINDER_PROBE,
            ]
        ]
