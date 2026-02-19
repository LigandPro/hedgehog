"""Tests for hedgehog.setup._shepherd_worker helper."""

from __future__ import annotations

from pathlib import Path

import pytest

from hedgehog.setup._shepherd_worker import ensure_shepherd_worker


class TestEnsureShepherdWorker:
    """Tests for the ensure_shepherd_worker function."""

    def test_raises_when_uv_is_missing(self, tmp_path: Path, monkeypatch) -> None:
        monkeypatch.setattr(
            "hedgehog.setup._shepherd_worker.shutil.which", lambda _: None
        )

        with pytest.raises(RuntimeError, match="uv is not installed"):
            ensure_shepherd_worker(tmp_path)

    def test_raises_when_no_supported_python(self, tmp_path: Path, monkeypatch) -> None:
        def _which(name: str) -> str | None:
            if name == "uv":
                return "/usr/bin/uv"
            return None

        monkeypatch.setattr("hedgehog.setup._shepherd_worker.shutil.which", _which)

        with pytest.raises(RuntimeError, match="No supported Python interpreter found"):
            ensure_shepherd_worker(tmp_path)

    def test_installs_with_explicit_python(self, tmp_path: Path, monkeypatch) -> None:
        explicit_python = tmp_path / "python3.12"
        explicit_python.write_text("#!/bin/sh\n", encoding="utf-8")

        monkeypatch.setattr(
            "hedgehog.setup._shepherd_worker.shutil.which",
            lambda name: "/usr/bin/uv" if name == "uv" else None,
        )
        monkeypatch.setattr(
            "hedgehog.setup._shepherd_worker.confirm_download",
            lambda *_a, **_kw: True,
        )

        calls: list[list[str]] = []

        def _mock_run(cmd, cwd=None, check=False, timeout=None):
            calls.append([str(x) for x in cmd])
            if cmd[:3] == [str(explicit_python), "-m", "venv"]:
                venv_dir = Path(cmd[3])
                bin_dir = venv_dir / "bin"
                bin_dir.mkdir(parents=True, exist_ok=True)
                (bin_dir / "python").write_text("#!/bin/sh\n", encoding="utf-8")
            if cmd[:3] == ["uv", "pip", "install"]:
                venv_python = Path(cmd[cmd.index("--python") + 1])
                worker_entry = venv_python.parent / "hedgehog-shepherd-worker"
                worker_entry.write_text("#!/bin/sh\n", encoding="utf-8")

        monkeypatch.setattr("hedgehog.setup._shepherd_worker.subprocess.run", _mock_run)

        result = ensure_shepherd_worker(tmp_path, python_bin=str(explicit_python))

        assert (
            result
            == tmp_path / ".venv-shepherd-worker" / "bin" / "hedgehog-shepherd-worker"
        )
        assert calls[0][:3] == [str(explicit_python), "-m", "venv"]
        assert calls[1][:3] == ["uv", "pip", "install"]
        assert ".[shepherd]" in calls[1]

    def test_selects_first_available_python(self, tmp_path: Path, monkeypatch) -> None:
        def _which(name: str) -> str | None:
            if name == "uv":
                return "/usr/bin/uv"
            if name == "python3.12":
                return "/usr/local/bin/python3.12"
            if name == "python3.11":
                return "/usr/local/bin/python3.11"
            if name == "python3.10":
                return "/usr/local/bin/python3.10"
            return None

        monkeypatch.setattr("hedgehog.setup._shepherd_worker.shutil.which", _which)
        monkeypatch.setattr(
            "hedgehog.setup._shepherd_worker.confirm_download",
            lambda *_a, **_kw: True,
        )

        calls: list[list[str]] = []

        def _mock_run(cmd, cwd=None, check=False, timeout=None):
            calls.append([str(x) for x in cmd])
            if cmd[:3] == ["/usr/local/bin/python3.12", "-m", "venv"]:
                venv_dir = Path(cmd[3])
                bin_dir = venv_dir / "bin"
                bin_dir.mkdir(parents=True, exist_ok=True)
                (bin_dir / "python").write_text("#!/bin/sh\n", encoding="utf-8")
            if cmd[:3] == ["uv", "pip", "install"]:
                venv_python = Path(cmd[cmd.index("--python") + 1])
                worker_entry = venv_python.parent / "hedgehog-shepherd-worker"
                worker_entry.write_text("#!/bin/sh\n", encoding="utf-8")

        monkeypatch.setattr("hedgehog.setup._shepherd_worker.subprocess.run", _mock_run)

        ensure_shepherd_worker(tmp_path)
        assert calls[0][:3] == ["/usr/local/bin/python3.12", "-m", "venv"]
