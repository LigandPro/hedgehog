"""Tests for hedgehog.setup._download utilities."""

from __future__ import annotations

import io
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from hedgehog.setup import _download
from hedgehog.setup._download import (
    confirm_download,
    download_with_progress,
    resolve_uv_binary,
)

# ---------------------------------------------------------------------------
# confirm_download
# ---------------------------------------------------------------------------


class TestConfirmDownload:
    """Tests for the confirm_download interactive prompt."""

    def test_auto_install_env_var(self, monkeypatch):
        """HEDGEHOG_AUTO_INSTALL=1 should return True without prompting."""
        monkeypatch.setenv("HEDGEHOG_AUTO_INSTALL", "1")
        # input() should never be called
        monkeypatch.setattr(
            "builtins.input", lambda _: (_ for _ in ()).throw(AssertionError)
        )
        assert confirm_download("gnina", "250 MB") is True

    def test_non_interactive_env_var(self, monkeypatch):
        """HEDGEHOG_NON_INTERACTIVE=1 should return False without prompting."""
        monkeypatch.setenv("HEDGEHOG_NON_INTERACTIVE", "1")
        monkeypatch.delenv("HEDGEHOG_AUTO_INSTALL", raising=False)
        monkeypatch.setattr(
            "builtins.input", lambda _: (_ for _ in ()).throw(AssertionError)
        )
        assert confirm_download("gnina", "250 MB") is False

    def test_user_accepts(self, monkeypatch):
        """User typing 'y' should return True."""
        monkeypatch.delenv("HEDGEHOG_AUTO_INSTALL", raising=False)
        monkeypatch.delenv("HEDGEHOG_NON_INTERACTIVE", raising=False)
        monkeypatch.setattr("builtins.input", lambda _: "y")
        assert confirm_download("gnina", "250 MB") is True

    def test_user_declines(self, monkeypatch):
        """User typing 'n' should return False."""
        monkeypatch.delenv("HEDGEHOG_AUTO_INSTALL", raising=False)
        monkeypatch.delenv("HEDGEHOG_NON_INTERACTIVE", raising=False)
        monkeypatch.setattr("builtins.input", lambda _: "n")
        assert confirm_download("gnina", "250 MB") is False

    def test_empty_input_accepts(self, monkeypatch):
        """Pressing Enter (empty input) should default to True."""
        monkeypatch.delenv("HEDGEHOG_AUTO_INSTALL", raising=False)
        monkeypatch.delenv("HEDGEHOG_NON_INTERACTIVE", raising=False)
        monkeypatch.setattr("builtins.input", lambda _: "")
        assert confirm_download("gnina", "250 MB") is True

    def test_eof_error_declines(self, monkeypatch):
        """EOFError (piped stdin) should return False."""
        monkeypatch.delenv("HEDGEHOG_AUTO_INSTALL", raising=False)
        monkeypatch.delenv("HEDGEHOG_NON_INTERACTIVE", raising=False)

        def _raise_eof(_prompt):
            raise EOFError

        monkeypatch.setattr("builtins.input", _raise_eof)
        assert confirm_download("gnina", "250 MB") is False


# ---------------------------------------------------------------------------
# download_with_progress
# ---------------------------------------------------------------------------


def _fake_urlopen(data: bytes, content_length: int | None = None):
    """Return a mock urlopen that yields *data* and reports Content-Length."""
    body = io.BytesIO(data)
    resp = MagicMock()
    resp.read = body.read
    resp.headers = {"Content-Length": str(content_length or len(data))}
    return resp


class TestDownloadWithProgress:
    """Tests for the download_with_progress function."""

    def test_successful_download(self, monkeypatch, tmp_path: Path):
        """File should contain the expected bytes after a successful download."""
        payload = b"hello world" * 100
        monkeypatch.setattr(
            "hedgehog.setup._download.urllib.request.urlopen",
            lambda _url: _fake_urlopen(payload),
        )

        dest = tmp_path / "tool.bin"
        download_with_progress("https://example.com/tool.bin", dest, "Downloading tool")

        assert dest.exists()
        assert dest.read_bytes() == payload
        # Temp file should be cleaned up
        assert not dest.with_suffix(dest.suffix + ".downloading").exists()

    def test_atomic_cleanup_on_error(self, monkeypatch, tmp_path: Path):
        """Temp .downloading file should be removed if the download fails."""

        def _exploding_urlopen(_url):
            raise ConnectionError("simulated network failure")

        monkeypatch.setattr(
            "hedgehog.setup._download.urllib.request.urlopen",
            _exploding_urlopen,
        )

        dest = tmp_path / "tool.bin"
        with pytest.raises(ConnectionError, match="simulated network failure"):
            download_with_progress(
                "https://example.com/tool.bin", dest, "Downloading tool"
            )

        assert not dest.exists()
        assert not dest.with_suffix(dest.suffix + ".downloading").exists()

    def test_creates_parent_dirs(self, monkeypatch, tmp_path: Path):
        """Parent directories should be created automatically."""
        payload = b"data"
        monkeypatch.setattr(
            "hedgehog.setup._download.urllib.request.urlopen",
            lambda _url: _fake_urlopen(payload),
        )

        dest = tmp_path / "a" / "b" / "c" / "tool.bin"
        download_with_progress("https://example.com/tool.bin", dest, "Downloading tool")

        assert dest.exists()
        assert dest.read_bytes() == payload


# ---------------------------------------------------------------------------
# resolve_uv_binary
# ---------------------------------------------------------------------------


def _write_executable(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/sh\n", encoding="utf-8")
    path.chmod(0o755)


def test_prefers_uv_env(monkeypatch, tmp_path: Path):
    """Resolve uv from the UV env var when it points to an executable."""
    uv_bin = tmp_path / "bin" / "uv"
    _write_executable(uv_bin)
    monkeypatch.setenv("UV", str(uv_bin))
    monkeypatch.setattr(_download.shutil, "which", lambda name: "/usr/bin/uv")

    assert resolve_uv_binary() == str(uv_bin)


def test_falls_back_to_path_when_uv_env_invalid(monkeypatch):
    """When UV points to a missing file, fall back to shutil.which."""
    monkeypatch.setenv("UV", "/nonexistent/uv")
    monkeypatch.setattr(_download.shutil, "which", lambda name: "/mock/path/uv")

    assert resolve_uv_binary() == "/mock/path/uv"


def test_raises_if_uv_absent(monkeypatch):
    """Missing UV env and PATH should raise a RuntimeError."""
    monkeypatch.delenv("UV", raising=False)
    monkeypatch.setattr(_download.shutil, "which", lambda name: None)

    with pytest.raises(RuntimeError, match="uv is not installed"):
        resolve_uv_binary()
