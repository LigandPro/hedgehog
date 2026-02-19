"""Tests for GNINA auto-installation module."""

from __future__ import annotations

import subprocess
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from hedgehog.setup._gnina import (
    _is_working_gnina,
    _resolve_gnina_download,
    ensure_gnina,
)


class TestEnsureGnina:
    """Tests for the ensure_gnina() resolution chain."""

    def test_found_in_path(self, monkeypatch):
        """If gnina is on PATH and working, return it immediately."""
        monkeypatch.setattr(
            "hedgehog.setup._gnina.shutil.which", lambda _: "/usr/bin/gnina"
        )
        monkeypatch.setattr("hedgehog.setup._gnina._is_working_gnina", lambda _: True)

        assert ensure_gnina() == "/usr/bin/gnina"

    def test_found_in_path_but_broken_falls_through(self, monkeypatch, tmp_path):
        """If PATH binary is broken, fall through to cache check."""
        monkeypatch.setattr(
            "hedgehog.setup._gnina.shutil.which", lambda _: "/usr/bin/gnina"
        )

        call_count = {"n": 0}

        def mock_is_working(path):
            call_count["n"] += 1
            if path == "/usr/bin/gnina":
                return False
            return True

        monkeypatch.setattr("hedgehog.setup._gnina._is_working_gnina", mock_is_working)

        # Set up a valid cached binary
        cache_dir = tmp_path / "bin"
        cache_dir.mkdir()
        cached = cache_dir / "gnina"
        cached.write_bytes(b"\x00" * 2_000_000)
        monkeypatch.setattr("hedgehog.setup._gnina._GNINA_CACHE_DIR", cache_dir)

        result = ensure_gnina()
        assert result == str(cached)
        assert call_count["n"] == 2  # Called once for PATH, once for cache

    def test_found_in_cache(self, monkeypatch, tmp_path):
        """If gnina is cached at ~/.hedgehog/bin/gnina, return it."""
        monkeypatch.setattr("hedgehog.setup._gnina.shutil.which", lambda _: None)

        cache_dir = tmp_path / "bin"
        cache_dir.mkdir()
        cached = cache_dir / "gnina"
        cached.write_bytes(b"\x00" * 2_000_000)

        monkeypatch.setattr("hedgehog.setup._gnina._GNINA_CACHE_DIR", cache_dir)
        monkeypatch.setattr("hedgehog.setup._gnina._is_working_gnina", lambda _: True)

        assert ensure_gnina() == str(cached)

    def test_cache_too_small_ignored(self, monkeypatch, tmp_path):
        """Cached binary smaller than 1MB should be ignored."""
        monkeypatch.setattr("hedgehog.setup._gnina.shutil.which", lambda _: None)

        cache_dir = tmp_path / "bin"
        cache_dir.mkdir()
        cached = cache_dir / "gnina"
        cached.write_bytes(b"\x00" * 100)  # Too small

        monkeypatch.setattr("hedgehog.setup._gnina._GNINA_CACHE_DIR", cache_dir)
        monkeypatch.setattr("hedgehog.setup._gnina.sys", MagicMock(platform="darwin"))

        with pytest.raises(RuntimeError, match="only available on Linux"):
            ensure_gnina()

    def test_not_linux_raises(self, monkeypatch):
        """On non-Linux platforms, raise RuntimeError."""
        monkeypatch.setattr("hedgehog.setup._gnina.shutil.which", lambda _: None)
        monkeypatch.setattr(
            "hedgehog.setup._gnina._GNINA_CACHE_DIR", Path("/nonexistent")
        )
        monkeypatch.setattr("hedgehog.setup._gnina.sys", MagicMock(platform="darwin"))

        with pytest.raises(RuntimeError, match="only available on Linux"):
            ensure_gnina()

    def test_user_declines_raises(self, monkeypatch, tmp_path):
        """If user declines download, raise RuntimeError."""
        monkeypatch.setattr("hedgehog.setup._gnina.shutil.which", lambda _: None)
        monkeypatch.setattr(
            "hedgehog.setup._gnina._GNINA_CACHE_DIR", Path("/nonexistent")
        )
        monkeypatch.setattr("hedgehog.setup._gnina.sys", MagicMock(platform="linux"))
        monkeypatch.setattr("hedgehog.setup._gnina.confirm_download", lambda *a: False)

        with pytest.raises(RuntimeError, match="declined"):
            ensure_gnina()

    def test_successful_download(self, monkeypatch, tmp_path):
        """Full download flow: which=None, no cache, linux, user accepts."""
        monkeypatch.setattr("hedgehog.setup._gnina.shutil.which", lambda _: None)

        cache_dir = tmp_path / "bin"
        monkeypatch.setattr("hedgehog.setup._gnina._GNINA_CACHE_DIR", cache_dir)
        monkeypatch.setattr("hedgehog.setup._gnina.sys", MagicMock(platform="linux"))
        monkeypatch.setattr("hedgehog.setup._gnina.confirm_download", lambda *a: True)
        monkeypatch.setattr(
            "hedgehog.setup._gnina._resolve_gnina_download",
            lambda: "https://example.com/gnina",
        )

        def fake_download(url, dest, desc):
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(b"\x00" * 2_000_000)

        monkeypatch.setattr(
            "hedgehog.setup._gnina.download_with_progress", fake_download
        )
        monkeypatch.setattr("hedgehog.setup._gnina._is_working_gnina", lambda _: True)

        result = ensure_gnina()
        assert result == str(cache_dir / "gnina")

    def test_download_fails_verification(self, monkeypatch, tmp_path):
        """If downloaded binary fails --version check, delete and raise."""
        monkeypatch.setattr("hedgehog.setup._gnina.shutil.which", lambda _: None)

        cache_dir = tmp_path / "bin"
        monkeypatch.setattr("hedgehog.setup._gnina._GNINA_CACHE_DIR", cache_dir)
        monkeypatch.setattr("hedgehog.setup._gnina.sys", MagicMock(platform="linux"))
        monkeypatch.setattr("hedgehog.setup._gnina.confirm_download", lambda *a: True)
        monkeypatch.setattr(
            "hedgehog.setup._gnina._resolve_gnina_download",
            lambda: "https://example.com/gnina",
        )

        def fake_download(url, dest, desc):
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(b"\x00" * 2_000_000)

        monkeypatch.setattr(
            "hedgehog.setup._gnina.download_with_progress", fake_download
        )
        monkeypatch.setattr("hedgehog.setup._gnina._is_working_gnina", lambda _: False)

        with pytest.raises(RuntimeError, match="failed verification"):
            ensure_gnina()

        # Binary should be cleaned up
        assert not (cache_dir / "gnina").exists()


class TestIsWorkingGnina:
    """Tests for the _is_working_gnina() helper."""

    def test_success(self, monkeypatch):
        """Returns True when gnina --version exits with 0."""
        mock_result = MagicMock(returncode=0)
        monkeypatch.setattr(
            "hedgehog.setup._gnina.subprocess.run", lambda *a, **kw: mock_result
        )

        assert _is_working_gnina("/usr/bin/gnina") is True

    def test_nonzero_exit(self, monkeypatch):
        """Returns False when gnina --version exits with non-zero."""
        mock_result = MagicMock(returncode=1)
        monkeypatch.setattr(
            "hedgehog.setup._gnina.subprocess.run", lambda *a, **kw: mock_result
        )

        assert _is_working_gnina("/usr/bin/gnina") is False

    def test_timeout(self, monkeypatch):
        """Returns False on TimeoutExpired."""

        def raise_timeout(*a, **kw):
            raise subprocess.TimeoutExpired(cmd="gnina", timeout=10)

        monkeypatch.setattr("hedgehog.setup._gnina.subprocess.run", raise_timeout)

        assert _is_working_gnina("/usr/bin/gnina") is False

    def test_os_error(self, monkeypatch):
        """Returns False on OSError (e.g. binary not found)."""

        def raise_os_error(*a, **kw):
            raise OSError("No such file")

        monkeypatch.setattr("hedgehog.setup._gnina.subprocess.run", raise_os_error)

        assert _is_working_gnina("/nonexistent/gnina") is False


class TestResolveGninaDownload:
    """Tests for _resolve_gnina_download() GitHub API interaction."""

    def test_selects_non_cuda_asset(self, monkeypatch):
        """Picks the first asset without 'cuda' in the name."""
        fake_response = {
            "assets": [
                {
                    "name": "gnina-v1.0-cuda11.tar.gz",
                    "browser_download_url": "https://example.com/cuda",
                },
                {
                    "name": "gnina",
                    "browser_download_url": "https://example.com/gnina-cpu",
                },
            ]
        }

        import json

        class FakeResponse:
            def read(self):
                return json.dumps(fake_response).encode()

            def __enter__(self):
                return self

            def __exit__(self, *a):
                pass

        monkeypatch.setattr(
            "hedgehog.setup._gnina.urllib.request.urlopen",
            lambda *a, **kw: FakeResponse(),
        )

        assert _resolve_gnina_download() == "https://example.com/gnina-cpu"

    def test_skips_archive_like_assets(self, monkeypatch):
        """Prefer a plain binary when archives are available."""
        fake_response = {
            "assets": [
                {
                    "name": "gnina-macosx.tar.gz",
                    "browser_download_url": "https://example.com/macos",
                },
                {
                    "name": "gnina-linux.tar.gz",
                    "browser_download_url": "https://example.com/linux-archive",
                },
                {
                    "name": "gnina.1.4.0",
                    "browser_download_url": "https://example.com/linux-binary",
                },
            ]
        }

        import json

        class FakeResponse:
            def read(self):
                return json.dumps(fake_response).encode()

            def __enter__(self):
                return self

            def __exit__(self, *a):
                pass

        monkeypatch.setattr(
            "hedgehog.setup._gnina.urllib.request.urlopen",
            lambda *a, **kw: FakeResponse(),
        )

        assert _resolve_gnina_download() == "https://example.com/linux-binary"

    def test_skips_other_platform_builds(self, monkeypatch):
        """Skip Mac/Windows builds even if they show up before Linux."""
        fake_response = {
            "assets": [
                {
                    "name": "gnina-mac",
                    "browser_download_url": "https://example.com/macos",
                },
                {
                    "name": "gnina-win64",
                    "browser_download_url": "https://example.com/windows",
                },
                {
                    "name": "gnina-linux",
                    "browser_download_url": "https://example.com/linux",
                },
            ]
        }

        import json

        class FakeResponse:
            def read(self):
                return json.dumps(fake_response).encode()

            def __enter__(self):
                return self

            def __exit__(self, *a):
                pass

        monkeypatch.setattr(
            "hedgehog.setup._gnina.urllib.request.urlopen",
            lambda *a, **kw: FakeResponse(),
        )

        assert _resolve_gnina_download() == "https://example.com/linux"

    def test_all_cuda_raises(self, monkeypatch):
        """RuntimeError if all assets contain 'cuda'."""
        fake_response = {
            "assets": [
                {
                    "name": "gnina-CUDA-only",
                    "browser_download_url": "https://example.com/cuda",
                },
            ]
        }

        import json

        class FakeResponse:
            def read(self):
                return json.dumps(fake_response).encode()

            def __enter__(self):
                return self

            def __exit__(self, *a):
                pass

        monkeypatch.setattr(
            "hedgehog.setup._gnina.urllib.request.urlopen",
            lambda *a, **kw: FakeResponse(),
        )

        with pytest.raises(RuntimeError, match="CUDA-specific"):
            _resolve_gnina_download()

    def test_api_failure_raises(self, monkeypatch):
        """RuntimeError if GitHub API call fails."""

        def raise_error(*a, **kw):
            raise OSError("Connection refused")

        monkeypatch.setattr(
            "hedgehog.setup._gnina.urllib.request.urlopen",
            raise_error,
        )

        with pytest.raises(RuntimeError, match="Failed to query"):
            _resolve_gnina_download()

    def test_falls_back_to_stable_tag_when_latest_too_large(self, monkeypatch):
        """When latest non-CUDA asset is oversized, fallback tag should be used."""
        import json

        latest_response = {
            "assets": [
                {
                    "name": "gnina.1.3.2",
                    "size": 1_426_790_536,
                    "browser_download_url": "https://example.com/latest-gnina",
                }
            ]
        }
        fallback_response = {
            "assets": [
                {
                    "name": "gnina",
                    "size": 306_470_832,
                    "browser_download_url": "https://example.com/v1.1-gnina",
                }
            ]
        }

        class FakeResponse:
            def __init__(self, payload):
                self.payload = payload

            def read(self):
                return json.dumps(self.payload).encode()

            def __enter__(self):
                return self

            def __exit__(self, *a):
                pass

        def fake_urlopen(req, *a, **kw):
            url = getattr(req, "full_url", str(req))
            if url.endswith("/latest"):
                return FakeResponse(latest_response)
            if url.endswith("/tags/v1.1"):
                return FakeResponse(fallback_response)
            raise AssertionError(f"Unexpected URL: {url}")

        monkeypatch.setattr(
            "hedgehog.setup._gnina.urllib.request.urlopen",
            fake_urlopen,
        )

        assert _resolve_gnina_download() == "https://example.com/v1.1-gnina"

    def test_size_limit_override_allows_large_latest_asset(self, monkeypatch):
        """Explicit size-limit override should allow selecting large latest asset."""
        import json

        latest_response = {
            "assets": [
                {
                    "name": "gnina.1.3.2",
                    "size": 1_426_790_536,
                    "browser_download_url": "https://example.com/latest-gnina",
                }
            ]
        }

        class FakeResponse:
            def read(self):
                return json.dumps(latest_response).encode()

            def __enter__(self):
                return self

            def __exit__(self, *a):
                pass

        monkeypatch.setenv("HEDGEHOG_GNINA_MAX_DOWNLOAD_BYTES", "3000000000")
        monkeypatch.setattr(
            "hedgehog.setup._gnina.urllib.request.urlopen",
            lambda *a, **kw: FakeResponse(),
        )

        assert _resolve_gnina_download() == "https://example.com/latest-gnina"

    def test_oversized_assets_raise_size_limit_error(self, monkeypatch):
        """Helpful error when all non-CUDA assets exceed max auto-install size."""
        import json

        latest_response = {
            "assets": [
                {
                    "name": "gnina.1.3.2",
                    "size": 1_426_790_536,
                    "browser_download_url": "https://example.com/latest-gnina",
                }
            ]
        }
        fallback_response = {
            "assets": [
                {
                    "name": "gnina",
                    "size": 1_100_000_000,
                    "browser_download_url": "https://example.com/v1.1-gnina",
                }
            ]
        }

        class FakeResponse:
            def __init__(self, payload):
                self.payload = payload

            def read(self):
                return json.dumps(self.payload).encode()

            def __enter__(self):
                return self

            def __exit__(self, *a):
                pass

        def fake_urlopen(req, *a, **kw):
            url = getattr(req, "full_url", str(req))
            if url.endswith("/latest"):
                return FakeResponse(latest_response)
            if url.endswith("/tags/v1.1"):
                return FakeResponse(fallback_response)
            raise AssertionError(f"Unexpected URL: {url}")

        monkeypatch.setattr(
            "hedgehog.setup._gnina.urllib.request.urlopen",
            fake_urlopen,
        )

        with pytest.raises(RuntimeError, match="size limit"):
            _resolve_gnina_download()

    def test_gpu_variant_prefers_cuda_asset(self, monkeypatch):
        """GPU variant should select a CUDA-labeled Linux asset when present."""
        import json

        fake_response = {
            "assets": [
                {
                    "name": "gnina.1.3.2",
                    "size": 1_426_790_536,
                    "browser_download_url": "https://example.com/cpu",
                },
                {
                    "name": "gnina.1.3.2.cuda12.8",
                    "size": 2_052_029_472,
                    "browser_download_url": "https://example.com/cuda12.8",
                },
            ]
        }

        class FakeResponse:
            def read(self):
                return json.dumps(fake_response).encode()

            def __enter__(self):
                return self

            def __exit__(self, *a):
                pass

        monkeypatch.setenv("HEDGEHOG_GNINA_VARIANT", "gpu")
        monkeypatch.setattr(
            "hedgehog.setup._gnina.urllib.request.urlopen",
            lambda *a, **kw: FakeResponse(),
        )

        assert _resolve_gnina_download() == "https://example.com/cuda12.8"

    def test_auto_variant_uses_gpu_when_detected(self, monkeypatch):
        """Auto variant should prefer CUDA asset when NVIDIA GPU is detected."""
        import json

        fake_response = {
            "assets": [
                {
                    "name": "gnina.1.3.2",
                    "size": 1_426_790_536,
                    "browser_download_url": "https://example.com/cpu",
                },
                {
                    "name": "gnina.1.3.2.cuda12.8",
                    "size": 2_052_029_472,
                    "browser_download_url": "https://example.com/cuda12.8",
                },
            ]
        }

        class FakeResponse:
            def read(self):
                return json.dumps(fake_response).encode()

            def __enter__(self):
                return self

            def __exit__(self, *a):
                pass

        monkeypatch.setenv("HEDGEHOG_GNINA_VARIANT", "auto")
        monkeypatch.setattr("hedgehog.setup._gnina._has_nvidia_gpu", lambda: True)
        monkeypatch.setattr(
            "hedgehog.setup._gnina.urllib.request.urlopen",
            lambda *a, **kw: FakeResponse(),
        )

        assert _resolve_gnina_download() == "https://example.com/cuda12.8"


class TestResolveDockingBinaryGninaFallback:
    """Tests for gnina fallback in _resolve_docking_binary()."""

    def test_gnina_fallback_called(self, monkeypatch):
        """When gnina is not on PATH, ensure_gnina() is called as fallback."""
        monkeypatch.setattr("shutil.which", lambda _: None)
        monkeypatch.setattr(
            "hedgehog.setup._gnina.ensure_gnina",
            lambda: "/home/user/.hedgehog/bin/gnina",
        )

        from hedgehog.stages.docking.utils import _resolve_docking_binary

        result = _resolve_docking_binary("gnina", "gnina")
        assert result == "/home/user/.hedgehog/bin/gnina"

    def test_gnina_fallback_error_becomes_file_not_found(self, monkeypatch):
        """RuntimeError from ensure_gnina becomes FileNotFoundError."""
        monkeypatch.setattr("shutil.which", lambda _: None)

        def raise_runtime(*a, **kw):
            raise RuntimeError("GNINA download declined by user.")

        monkeypatch.setattr("hedgehog.setup._gnina.ensure_gnina", raise_runtime)

        from hedgehog.stages.docking.utils import _resolve_docking_binary

        with pytest.raises(FileNotFoundError, match="declined"):
            _resolve_docking_binary("gnina", "gnina")

    def test_smina_still_raises(self, monkeypatch):
        """smina has no auto-install fallback, so FileNotFoundError directly."""
        monkeypatch.setattr("shutil.which", lambda _: None)

        from hedgehog.stages.docking.utils import _resolve_docking_binary

        with pytest.raises(FileNotFoundError, match="smina"):
            _resolve_docking_binary("smina", "smina")
