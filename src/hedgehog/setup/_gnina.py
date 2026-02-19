"""Auto-installation of GNINA docking binary."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path
from typing import Any

from hedgehog.setup._download import confirm_download, download_with_progress

_GNINA_CACHE_DIR = Path.home() / ".hedgehog" / "bin"
_GNINA_RELEASES_API_LATEST = "https://api.github.com/repos/gnina/gnina/releases/latest"
_GNINA_RELEASES_API_TAG = "https://api.github.com/repos/gnina/gnina/releases/tags"
_GNINA_FALLBACK_TAG = "v1.1"
_GNINA_DEFAULT_MAX_DOWNLOAD_BYTES = 800 * 1024 * 1024
_GNINA_ARCHIVE_SUFFIXES = (
    ".tar.gz",
    ".tgz",
    ".tar",
    ".zip",
    ".gz",
    ".bz2",
    ".xz",
    ".7z",
)
_GNINA_OS_BLACKLIST = ("darwin", "mac", "osx", "win", "windows")


def ensure_gnina() -> str:
    """Resolve a working GNINA binary, downloading if necessary.

    Resolution order:
      1. ``gnina`` on PATH (verified working)
      2. Cached binary at ``~/.hedgehog/bin/gnina``
      3. Download from GitHub (Linux only, with user confirmation)

    Returns:
        Absolute path to a working GNINA binary.

    Raises:
        RuntimeError: If GNINA cannot be found or installed.
    """
    # 1. Check PATH
    path_binary = shutil.which("gnina")
    if path_binary and _is_working_gnina(path_binary):
        return path_binary

    # 2. Check cache
    cached = _GNINA_CACHE_DIR / "gnina"
    if (
        cached.is_file()
        and cached.stat().st_size > 1_000_000
        and _is_working_gnina(str(cached))
    ):
        return str(cached)

    # 3. Platform gate
    if sys.platform != "linux":
        raise RuntimeError(
            "GNINA is only available on Linux. Please install it manually."
        )

    # 4. Download with user confirmation
    size_hint = f"up to {_format_size(_gnina_max_download_bytes())}"
    if not confirm_download("GNINA", size_hint):
        raise RuntimeError("GNINA download declined by user.")

    url = _resolve_gnina_download()
    cached.parent.mkdir(parents=True, exist_ok=True)
    download_with_progress(url, cached, "GNINA")
    os.chmod(cached, 0o755)

    if not _is_working_gnina(str(cached)):
        cached.unlink(missing_ok=True)
        raise RuntimeError(
            "Downloaded GNINA binary failed verification (--version check). "
            "Please install GNINA manually."
        )

    return str(cached)


def _is_working_gnina(path: str) -> bool:
    """Return True if *path* runs ``gnina --version`` successfully."""
    try:
        result = subprocess.run(
            [path, "--version"],
            capture_output=True,
            timeout=10,
            env=_gnina_env(),
        )
        return result.returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        return False


def _resolve_gnina_download() -> str:
    """Resolve a GNINA download URL from GitHub releases.

    Selection strategy:
      1. Check ``releases/latest`` first.
      2. If all non-CUDA assets are too large, fallback to a stable tag.
      3. Enforce a max auto-download size limit.

    Raises:
        RuntimeError: If the API call fails or no suitable asset is found.
    """
    max_bytes = _gnina_max_download_bytes()
    release_urls = [
        _GNINA_RELEASES_API_LATEST,
        f"{_GNINA_RELEASES_API_TAG}/{_GNINA_FALLBACK_TAG}",
    ]
    query_errors: list[str] = []
    saw_non_cuda_asset = False
    oversize_assets: list[int] = []

    for api_url in release_urls:
        try:
            data = _query_release_json(api_url)
        except Exception as exc:
            query_errors.append(str(exc))
            continue

        for asset in _iter_ranked_non_cuda_assets(data):
            saw_non_cuda_asset = True
            size = _asset_size_bytes(asset)
            if size > 0 and size > max_bytes:
                oversize_assets.append(size)
                continue
            return str(asset["browser_download_url"])

    if not saw_non_cuda_asset:
        if query_errors and len(query_errors) == len(release_urls):
            raise RuntimeError(
                f"Failed to query GNINA releases from GitHub: {'; '.join(query_errors)}"
            )
        raise RuntimeError(
            "No suitable GNINA release asset found on GitHub (all assets appear CUDA-specific)."
        )

    if oversize_assets:
        smallest = min(oversize_assets)
        raise RuntimeError(
            "No suitable GNINA release asset found on GitHub under auto-install size "
            f"limit ({_format_size(max_bytes)}). Smallest non-CUDA asset was "
            f"{_format_size(smallest)}. Set HEDGEHOG_GNINA_MAX_DOWNLOAD_BYTES to "
            "override or install GNINA manually."
        )

    if query_errors:
        raise RuntimeError(
            f"Failed to query GNINA releases from GitHub: {'; '.join(query_errors)}"
        )

    raise RuntimeError("No suitable GNINA release asset found on GitHub.")


def _query_release_json(api_url: str) -> dict[str, Any]:
    """Fetch a GNINA release object from GitHub API."""
    req = urllib.request.Request(api_url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as resp:  # noqa: S310
        return json.loads(resp.read())


def _iter_ranked_non_cuda_assets(release_data: dict[str, Any]) -> list[dict[str, Any]]:
    """Return non-CUDA assets ranked by CPU-like naming and smaller size."""
    assets: list[dict[str, Any]] = []
    for raw in release_data.get("assets", []):
        if not isinstance(raw, dict):
            continue
        name = str(raw.get("name", "")).strip()
        url = raw.get("browser_download_url")
        if not name or not url:
            continue
        lowered = name.lower()
        if "cuda" in lowered:
            continue
        if not _asset_is_linux_binary(lowered):
            continue
        assets.append(raw)

    def _rank_key(asset: dict[str, Any]) -> tuple[int, int, str]:
        name = str(asset.get("name", "")).lower()
        cpu_hint = 0 if any(tok in name for tok in ("cpu", "nocuda", "no_cuda")) else 1
        size = _asset_size_bytes(asset)
        size_key = size if size > 0 else sys.maxsize
        return (cpu_hint, size_key, name)

    return sorted(assets, key=_rank_key)


def _asset_size_bytes(asset: dict[str, Any]) -> int:
    """Parse release asset size in bytes (0 if missing/invalid)."""
    try:
        size = int(asset.get("size", 0))
        return size if size > 0 else 0
    except Exception:
        return 0


def _asset_is_linux_binary(name: str) -> bool:
    """Heuristic that rejects archives and non-Linux builds."""
    if not name:
        return False
    lower = name.lower()
    if any(lower.endswith(suffix) for suffix in _GNINA_ARCHIVE_SUFFIXES):
        return False
    if any(os_tag in lower for os_tag in _GNINA_OS_BLACKLIST):
        return False
    return True


def _gnina_max_download_bytes() -> int:
    """Resolve GNINA auto-install size limit from environment."""
    raw = os.environ.get("HEDGEHOG_GNINA_MAX_DOWNLOAD_BYTES")
    if not raw:
        return _GNINA_DEFAULT_MAX_DOWNLOAD_BYTES
    try:
        value = int(raw)
        if value > 0:
            return value
    except ValueError:
        pass
    return _GNINA_DEFAULT_MAX_DOWNLOAD_BYTES


def _format_size(size_bytes: int) -> str:
    """Format byte size for user-facing messages."""
    if size_bytes <= 0:
        return "unknown size"
    units = ("B", "KB", "MB", "GB", "TB")
    value = float(size_bytes)
    for unit in units:
        if value < 1024 or unit == units[-1]:
            if unit in ("B", "KB"):
                return f"{int(round(value))} {unit}"
            return f"{value:.1f} {unit}"
        value /= 1024.0
    return f"{size_bytes} B"


def _gnina_env() -> dict[str, str]:
    """Build an environment dict with extended ``LD_LIBRARY_PATH`` for GNINA.

    GNINA links against several shared libraries (libcudnn, libtorch, etc.)
    that may live inside a PyTorch installation or a Conda environment.
    """
    env = os.environ.copy()
    extra_paths: list[str] = []

    # PyTorch bundled libraries
    try:
        import torch

        torch_lib = Path(torch.__file__).parent / "lib"
        if torch_lib.is_dir():
            extra_paths.append(str(torch_lib))
    except ImportError:
        pass

    # Conda environment library directories
    _CONDA_ROOTS = ("miniforge3", "miniconda3", "mambaforge", "anaconda3")
    home = Path.home()
    for root_name in _CONDA_ROOTS:
        lib_dir = home / root_name / "lib"
        if lib_dir.is_dir():
            extra_paths.append(str(lib_dir))

    if extra_paths:
        existing = env.get("LD_LIBRARY_PATH", "")
        env["LD_LIBRARY_PATH"] = os.pathsep.join(
            extra_paths + ([existing] if existing else [])
        )

    return env
