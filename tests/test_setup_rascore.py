"""Tests for RAScore model auto-download helper."""

from __future__ import annotations

from pathlib import Path

import pytest

from hedgehog.setup._rascore import ensure_rascore_model


class TestEnsureRascoreModel:
    """Tests for ensure_rascore_model()."""

    def test_uses_existing_model(self, tmp_path: Path) -> None:
        """Return existing model when file is present and large enough."""
        model_path = (
            tmp_path
            / "modules"
            / "MolScore"
            / "molscore"
            / "data"
            / "models"
            / "RAScore"
            / "XGB_chembl_ecfp_counts"
            / "model.pkl"
        )
        model_path.parent.mkdir(parents=True, exist_ok=True)
        model_path.write_bytes(b"x" * 1_200_000)

        result = ensure_rascore_model(tmp_path)
        assert result == model_path

    def test_downloads_when_missing(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        """Download model into the expected MolScore path when missing."""
        calls: list[tuple[str, Path, str]] = []

        def _fake_download(url: str, dest: Path, description: str) -> None:
            calls.append((url, dest, description))
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(b"x" * 1_500_000)

        monkeypatch.setattr(
            "hedgehog.setup._rascore.download_with_progress", _fake_download
        )

        result = ensure_rascore_model(tmp_path)
        assert result.exists()
        assert result.name == "model.pkl"
        assert calls, "download_with_progress should be called"

    def test_raises_on_invalid_download(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        """Raise when downloaded file is too small/corrupt."""
        def _fake_download(url: str, dest: Path, description: str) -> None:
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(b"tiny")

        monkeypatch.setattr(
            "hedgehog.setup._rascore.download_with_progress", _fake_download
        )

        with pytest.raises(RuntimeError, match="invalid"):
            ensure_rascore_model(tmp_path)
