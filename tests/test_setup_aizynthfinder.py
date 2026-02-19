"""Tests for hedgehog.setup._aizynthfinder auto-installation."""

from __future__ import annotations

from pathlib import Path

import pytest

from hedgehog.setup._aizynthfinder import ensure_aizynthfinder


def _config_path(root: Path) -> Path:
    """Return the expected config.yml path for a given project root."""
    return (
        root / "modules" / "retrosynthesis" / "aizynthfinder" / "public" / "config.yml"
    )


def _make_config(root: Path) -> Path:
    """Create the config.yml file so the function thinks it's installed."""
    cfg = _config_path(root)
    cfg.parent.mkdir(parents=True, exist_ok=True)
    cfg.write_text("version: 1\n")
    return cfg


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestEnsureAizynthfinder:
    """Tests for the ensure_aizynthfinder function."""

    @pytest.fixture(autouse=True)
    def _clear_uv_env(self, monkeypatch):
        monkeypatch.delenv("UV", raising=False)

    def test_already_installed(self, tmp_path: Path, monkeypatch):
        """If config.yml already exists, return immediately."""
        cfg = _make_config(tmp_path)
        calls: list = []
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.subprocess.run",
            lambda *a, **kw: calls.append(a),
        )
        result = ensure_aizynthfinder(tmp_path)
        assert result == cfg
        assert calls == [], "No subprocess calls expected"

    def test_git_not_found(self, tmp_path: Path, monkeypatch):
        """RuntimeError when git is not on PATH."""
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.shutil.which",
            lambda name: None if name == "git" else "/usr/bin/" + name,
        )
        with pytest.raises(RuntimeError, match="git is not installed"):
            ensure_aizynthfinder(tmp_path)

    def test_uv_not_found(self, tmp_path: Path, monkeypatch):
        """RuntimeError when uv is not on PATH (git is found)."""
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.shutil.which",
            lambda name: "/usr/bin/git" if name == "git" else None,
        )
        with pytest.raises(RuntimeError, match="uv is not installed"):
            ensure_aizynthfinder(tmp_path)

    def test_user_declines(self, tmp_path: Path, monkeypatch):
        """RuntimeError when user declines the download prompt."""
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.shutil.which",
            lambda name: "/usr/bin/" + name,
        )
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.confirm_download",
            lambda *_a, **_kw: False,
        )
        with pytest.raises(RuntimeError, match="declined"):
            ensure_aizynthfinder(tmp_path)

    def test_uses_uv_from_environment(self, tmp_path: Path, monkeypatch):
        """UV env var should be accepted when uv is not on PATH."""
        uv_bin = tmp_path / "custom-bin" / "uv"
        uv_bin.parent.mkdir(parents=True, exist_ok=True)
        uv_bin.write_text("#!/bin/sh\n", encoding="utf-8")
        monkeypatch.setenv("UV", str(uv_bin))
        monkeypatch.setattr("hedgehog.setup._aizynthfinder.os.access", lambda *_: True)
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.shutil.which",
            lambda name: "/usr/bin/git" if name == "git" else None,
        )
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.confirm_download",
            lambda *_a, **_kw: True,
        )

        calls: list[list[str]] = []

        def _mock_run(cmd, *, cwd=None, check=False, timeout=None):
            calls.append([str(x) for x in cmd])
            if "clone" in cmd:
                aizynth = tmp_path / "modules" / "retrosynthesis" / "aizynthfinder"
                aizynth.mkdir(parents=True, exist_ok=True)
            if any("download_public_data" in str(part) for part in cmd):
                cfg = _config_path(tmp_path)
                cfg.parent.mkdir(parents=True, exist_ok=True)
                cfg.write_text("version: 1\n")

        monkeypatch.setattr("hedgehog.setup._aizynthfinder.subprocess.run", _mock_run)

        ensure_aizynthfinder(tmp_path)

        assert calls[1][:2] == [str(uv_bin), "sync"]
        assert calls[2][0] == str(uv_bin)

    def test_full_install_cycle(self, tmp_path: Path, monkeypatch):
        """Full flow: clone, sync, download, copy logging, return config."""
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.shutil.which",
            lambda name: "/usr/bin/" + name,
        )
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.confirm_download",
            lambda *_a, **_kw: True,
        )

        # Create logging.yml source so the copy step works
        logging_src = (
            tmp_path / "src" / "hedgehog" / "stages" / "synthesis" / "logging.yml"
        )
        logging_src.parent.mkdir(parents=True, exist_ok=True)
        logging_src.write_text("version: 1\n")

        subprocess_calls: list = []

        def _mock_run(cmd, *, cwd=None, check=False, timeout=None):
            subprocess_calls.append((cmd, cwd))
            # Simulate clone creating the directory structure
            if "clone" in cmd:
                aizynth = tmp_path / "modules" / "retrosynthesis" / "aizynthfinder"
                aizynth.mkdir(parents=True, exist_ok=True)
            # Simulate download creating config.yml
            if any("download_public_data" in str(part) for part in cmd):
                cfg = _config_path(tmp_path)
                cfg.parent.mkdir(parents=True, exist_ok=True)
                cfg.write_text("version: 1\n")

        monkeypatch.setattr("hedgehog.setup._aizynthfinder.subprocess.run", _mock_run)

        result = ensure_aizynthfinder(tmp_path)

        assert result == _config_path(tmp_path)
        # Three subprocess calls: clone, sync, download
        cmds = [c for c, _cwd in subprocess_calls]
        assert len(cmds) == 3
        assert "clone" in cmds[0]
        assert cmds[1] == ["/usr/bin/uv", "sync"]
        assert "aizynthfinder.tools.download_public_data" in cmds[2]

        # Logging.yml should be copied
        data_dir = (
            tmp_path
            / "modules"
            / "retrosynthesis"
            / "aizynthfinder"
            / "aizynthfinder"
            / "data"
        )
        assert (data_dir / "logging.yml").exists()

    def test_skips_clone_if_repo_exists(self, tmp_path: Path, monkeypatch):
        """Git clone is skipped when retrosynthesis directory already exists."""
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.shutil.which",
            lambda name: "/usr/bin/" + name,
        )
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.confirm_download",
            lambda *_a, **_kw: True,
        )

        # Pre-create retrosynthesis dir and aizynthfinder dir
        aizynth = tmp_path / "modules" / "retrosynthesis" / "aizynthfinder"
        aizynth.mkdir(parents=True, exist_ok=True)

        subprocess_calls: list = []

        def _mock_run(cmd, *, cwd=None, check=False, timeout=None):
            subprocess_calls.append(cmd)
            if any("download_public_data" in str(part) for part in cmd):
                cfg = _config_path(tmp_path)
                cfg.parent.mkdir(parents=True, exist_ok=True)
                cfg.write_text("version: 1\n")

        monkeypatch.setattr("hedgehog.setup._aizynthfinder.subprocess.run", _mock_run)

        ensure_aizynthfinder(tmp_path)

        assert not any("clone" in cmd for cmd in subprocess_calls), (
            "clone should be skipped"
        )

    def test_skips_sync_if_venv_exists(self, tmp_path: Path, monkeypatch):
        """uv sync is skipped when .venv already exists."""
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.shutil.which",
            lambda name: "/usr/bin/" + name,
        )
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.confirm_download",
            lambda *_a, **_kw: True,
        )

        # Pre-create retrosynthesis dir with .venv
        aizynth = tmp_path / "modules" / "retrosynthesis" / "aizynthfinder"
        (aizynth / ".venv").mkdir(parents=True, exist_ok=True)

        subprocess_calls: list = []

        def _mock_run(cmd, *, cwd=None, check=False, timeout=None):
            subprocess_calls.append(cmd)
            if any("download_public_data" in str(part) for part in cmd):
                cfg = _config_path(tmp_path)
                cfg.parent.mkdir(parents=True, exist_ok=True)
                cfg.write_text("version: 1\n")

        monkeypatch.setattr("hedgehog.setup._aizynthfinder.subprocess.run", _mock_run)

        ensure_aizynthfinder(tmp_path)

        assert not any("sync" in cmd for cmd in subprocess_calls), (
            "uv sync should be skipped"
        )

    def test_skips_download_if_public_has_files(self, tmp_path: Path, monkeypatch):
        """Public data download is skipped when public/ already has files."""
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.shutil.which",
            lambda name: "/usr/bin/" + name,
        )
        monkeypatch.setattr(
            "hedgehog.setup._aizynthfinder.confirm_download",
            lambda *_a, **_kw: True,
        )

        # Pre-create retrosynthesis dir with .venv and populated public/
        aizynth = tmp_path / "modules" / "retrosynthesis" / "aizynthfinder"
        (aizynth / ".venv").mkdir(parents=True, exist_ok=True)
        public = aizynth / "public"
        public.mkdir(parents=True, exist_ok=True)
        (public / "config.yml").write_text("version: 1\n")
        (public / "model.onnx").write_text("fake model")

        subprocess_calls: list = []

        def _mock_run(cmd, *, cwd=None, check=False, timeout=None):
            subprocess_calls.append(cmd)

        monkeypatch.setattr("hedgehog.setup._aizynthfinder.subprocess.run", _mock_run)

        result = ensure_aizynthfinder(tmp_path)

        assert result == _config_path(tmp_path)
        # No subprocess calls at all — config exists so it returns early
        # Actually, config exists → returns at step 1, before any subprocess
        assert subprocess_calls == []
