"""Tests for hedgehog.setup._nvmolkit_worker helper."""

from __future__ import annotations

from pathlib import Path

import pytest

from hedgehog.setup._nvmolkit_worker import ensure_nvmolkit_worker


class TestEnsureNvMolKitWorker:
    """Tests for the ensure_nvmolkit_worker function."""

    @pytest.fixture(autouse=True)
    def _clear_uv_env(self, monkeypatch) -> None:
        monkeypatch.delenv("UV", raising=False)

    def test_raises_when_uv_is_missing(self, tmp_path: Path, monkeypatch) -> None:
        monkeypatch.setattr(
            "hedgehog.setup._nvmolkit_worker.resolve_uv_binary",
            lambda: (_ for _ in ()).throw(RuntimeError("uv is not installed")),
        )

        with pytest.raises(RuntimeError, match="uv is not installed"):
            ensure_nvmolkit_worker(tmp_path)

    def test_raises_when_no_supported_python(self, tmp_path: Path, monkeypatch) -> None:
        monkeypatch.setattr(
            "hedgehog.setup._nvmolkit_worker.resolve_uv_binary",
            lambda: "/usr/bin/uv",
        )
        monkeypatch.setattr(
            "hedgehog.setup._nvmolkit_worker.shutil.which", lambda _: None
        )

        with pytest.raises(RuntimeError, match="No supported Python interpreter found"):
            ensure_nvmolkit_worker(tmp_path)

    def test_installs_with_explicit_python(self, tmp_path: Path, monkeypatch) -> None:
        explicit_python = tmp_path / "python3.12"
        explicit_python.write_text("#!/bin/sh\n", encoding="utf-8")

        monkeypatch.setattr(
            "hedgehog.setup._nvmolkit_worker.resolve_uv_binary",
            lambda: "/usr/bin/uv",
        )
        monkeypatch.setattr(
            "hedgehog.setup._nvmolkit_worker.shutil.which",
            lambda name: None,
        )
        monkeypatch.setattr(
            "hedgehog.setup._nvmolkit_worker.confirm_download",
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

        monkeypatch.setattr("hedgehog.setup._nvmolkit_worker.subprocess.run", _mock_run)

        result = ensure_nvmolkit_worker(tmp_path, python_bin=str(explicit_python))

        assert result == tmp_path / ".venv-nvmolkit-worker" / "bin" / "python"
        assert calls[0][:3] == [str(explicit_python), "-m", "venv"]
        assert calls[1][:3] == ["/usr/bin/uv", "pip", "install"]
        assert "nvmolkit" in calls[1]
        assert calls[2][:2] == [
            str(tmp_path / ".venv-nvmolkit-worker" / "bin" / "python"),
            "-c",
        ]
