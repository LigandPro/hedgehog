"""Auto-installation of GNINA docking binary."""

from __future__ import annotations

import json
import os
import re
import shutil
import site
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
_GNINA_DEFAULT_MAX_DOWNLOAD_BYTES_GPU = 3 * 1024 * 1024 * 1024
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
    variant = _gnina_variant()
    select_mode = _select_mode_for_variant(variant)
    size_hint = f"up to {_format_size(_gnina_max_download_bytes(select_mode))}"
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
    variant = _gnina_variant()
    select_mode = _select_mode_for_variant(variant)
    max_bytes = _gnina_max_download_bytes(select_mode)
    release_urls = [
        _GNINA_RELEASES_API_LATEST,
        f"{_GNINA_RELEASES_API_TAG}/{_GNINA_FALLBACK_TAG}",
    ]
    query_errors: list[str] = []
    saw_any_asset = False
    saw_cuda_asset = False
    saw_non_cuda_asset = False
    oversize_assets: list[int] = []

    for api_url in release_urls:
        try:
            data = _query_release_json(api_url)
        except Exception as exc:
            query_errors.append(str(exc))
            continue

        for asset in _iter_ranked_assets(data, select_mode):
            saw_any_asset = True
            is_cuda = _is_cuda_asset_name(str(asset.get("name", "")))
            if is_cuda:
                saw_cuda_asset = True
            else:
                saw_non_cuda_asset = True
            size = _asset_size_bytes(asset)
            if size > 0 and size > max_bytes:
                oversize_assets.append(size)
                continue
            return str(asset["browser_download_url"])

    if not saw_any_asset:
        if query_errors and len(query_errors) == len(release_urls):
            raise RuntimeError(
                f"Failed to query GNINA releases from GitHub: {'; '.join(query_errors)}"
            )
        if select_mode == "gpu":
            raise RuntimeError(
                "No suitable GNINA CUDA release asset found on GitHub. "
                "Set HEDGEHOG_GNINA_VARIANT=auto or cpu to allow CPU fallback."
            )
        if select_mode == "cpu":
            raise RuntimeError(
                "No suitable GNINA release asset found on GitHub (all assets appear CUDA-specific)."
            )
        raise RuntimeError("No suitable GNINA release asset found on GitHub.")

    if oversize_assets:
        smallest = min(oversize_assets)
        asset_label = (
            "CUDA-compatible asset"
            if select_mode in ("gpu", "auto_gpu")
            else "non-CUDA asset"
        )
        raise RuntimeError(
            "No suitable GNINA release asset found on GitHub under auto-install size "
            f"limit ({_format_size(max_bytes)}). Smallest {asset_label} was "
            f"{_format_size(smallest)}. Set HEDGEHOG_GNINA_MAX_DOWNLOAD_BYTES to "
            "override or install GNINA manually."
        )

    if select_mode == "gpu" and saw_non_cuda_asset and not saw_cuda_asset:
        raise RuntimeError(
            "GNINA GPU mode requested, but no CUDA asset was selected. "
            "Set HEDGEHOG_GNINA_VARIANT=auto or cpu to allow CPU fallback."
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


def _iter_ranked_assets(
    release_data: dict[str, Any], mode: str
) -> list[dict[str, Any]]:
    """Return Linux assets ranked according to selection mode.

    Modes:
      - ``cpu``: non-CUDA assets only
      - ``gpu``: CUDA assets only (prefer newest CUDA variant)
      - ``auto_gpu``: prefer CUDA assets first, then CPU fallback
    """
    assets: list[tuple[dict[str, Any], bool, tuple[int, int]]] = []
    for raw in release_data.get("assets", []):
        if not isinstance(raw, dict):
            continue
        name = str(raw.get("name", "")).strip()
        url = raw.get("browser_download_url")
        if not name or not url:
            continue
        lowered = name.lower()
        if not _asset_is_linux_binary(lowered):
            continue
        is_cuda = _is_cuda_asset_name(lowered)
        if mode == "cpu" and is_cuda:
            continue
        if mode == "gpu" and not is_cuda:
            continue
        assets.append((raw, is_cuda, _extract_cuda_version(lowered)))

    def _rank_key(item: tuple[dict[str, Any], bool, tuple[int, int]]) -> tuple:
        asset, is_cuda, cuda_version = item
        name = str(asset.get("name", "")).lower()
        cpu_hint = 0 if any(tok in name for tok in ("cpu", "nocuda", "no_cuda")) else 1
        size = _asset_size_bytes(asset)
        size_key = size if size > 0 else sys.maxsize
        if mode == "gpu":
            # Prefer higher CUDA version first, then smaller binaries.
            return (-cuda_version[0], -cuda_version[1], size_key, name)
        if mode == "auto_gpu":
            # Prefer CUDA variants first; among CUDA assets prefer newer.
            cuda_rank = 0 if is_cuda else 1
            return (
                cuda_rank,
                -cuda_version[0],
                -cuda_version[1],
                cpu_hint,
                size_key,
                name,
            )
        return (cpu_hint, size_key, name)

    return [item[0] for item in sorted(assets, key=_rank_key)]


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


def _gnina_max_download_bytes(mode: str) -> int:
    """Resolve GNINA auto-install size limit from environment."""
    raw = os.environ.get("HEDGEHOG_GNINA_MAX_DOWNLOAD_BYTES")
    if not raw:
        if mode in ("gpu", "auto_gpu"):
            return _GNINA_DEFAULT_MAX_DOWNLOAD_BYTES_GPU
        return _GNINA_DEFAULT_MAX_DOWNLOAD_BYTES
    try:
        value = int(raw)
        if value > 0:
            return value
    except ValueError:
        pass
    if mode in ("gpu", "auto_gpu"):
        return _GNINA_DEFAULT_MAX_DOWNLOAD_BYTES_GPU
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


def _gnina_variant() -> str:
    """Resolve requested GNINA install variant from environment.

    Supported values:
      - ``cpu``: CPU/non-CUDA assets only
      - ``gpu``: CUDA assets only
      - ``auto``: prefer CUDA assets when a GPU is detected
    """
    raw = str(os.environ.get("HEDGEHOG_GNINA_VARIANT", "auto")).strip().lower()
    if raw in {"cpu", "gpu", "auto"}:
        return raw
    return "auto"


def _select_mode_for_variant(variant: str) -> str:
    """Map user variant selection to internal asset selection mode."""
    if variant == "cpu":
        return "cpu"
    if variant == "gpu":
        return "gpu"
    return "auto_gpu" if _has_nvidia_gpu() else "cpu"


def _has_nvidia_gpu() -> bool:
    """Return True when a visible NVIDIA GPU is detected via nvidia-smi."""
    nvidia_smi = shutil.which("nvidia-smi")
    if not nvidia_smi:
        return False
    try:
        result = subprocess.run(
            [nvidia_smi, "-L"],
            capture_output=True,
            text=True,
            timeout=5,
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return result.returncode == 0 and bool(result.stdout.strip())


def _is_cuda_asset_name(name: str) -> bool:
    """Return True if release asset name indicates a CUDA build."""
    return "cuda" in name.lower()


def _extract_cuda_version(name: str) -> tuple[int, int]:
    """Extract CUDA major/minor version from an asset name.

    Examples:
      - ``gnina.1.3.2.cuda12.8`` -> ``(12, 8)``
      - ``gnina.cuda11`` -> ``(11, 0)``
    """
    match = re.search(r"cuda(\d+)(?:\.(\d+))?", name.lower())
    if not match:
        return (0, 0)
    major = int(match.group(1))
    minor = int(match.group(2) or 0)
    return (major, minor)


def _gnina_env() -> dict[str, str]:
    """Build an environment dict with extended ``LD_LIBRARY_PATH`` for GNINA.

    GNINA links against several shared libraries (libcudnn, libtorch, etc.)
    that may live inside a PyTorch installation or a Conda environment.
    """
    env = os.environ.copy()
    extra_paths = _collect_gnina_library_paths()
    existing_parts = _split_library_path(env.get("LD_LIBRARY_PATH", ""))
    all_parts = _dedupe_library_paths(extra_paths + existing_parts)

    if all_parts:
        env["LD_LIBRARY_PATH"] = os.pathsep.join(all_parts)

    return env


def _collect_gnina_library_paths() -> list[str]:
    """Collect candidate library directories required by GNINA."""
    extra_paths: list[str] = []

    # PyTorch bundled libraries
    try:
        import torch

        torch_lib = Path(torch.__file__).parent / "lib"
        if torch_lib.is_dir():
            extra_paths.append(str(torch_lib))
    except ImportError:
        pass

    # NVIDIA Python package libraries (e.g. nvidia/cudnn/lib, nvidia/cublas/lib)
    for site_packages in _iter_site_packages_dirs():
        nvidia_root = site_packages / "nvidia"
        if not nvidia_root.is_dir():
            continue
        for child in sorted(nvidia_root.iterdir()):
            lib_dir = child / "lib"
            if lib_dir.is_dir():
                extra_paths.append(str(lib_dir))

    # Active conda environment first
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        conda_prefix_lib = Path(conda_prefix) / "lib"
        if conda_prefix_lib.is_dir():
            extra_paths.append(str(conda_prefix_lib))

    # Common conda roots in home directory
    conda_roots = ("miniforge", "miniforge3", "miniconda3", "mambaforge", "anaconda3")
    home = Path.home()
    for root_name in conda_roots:
        lib_dir = home / root_name / "lib"
        if lib_dir.is_dir():
            extra_paths.append(str(lib_dir))

    return _dedupe_library_paths(extra_paths)


def _iter_site_packages_dirs() -> list[Path]:
    """Return existing site-packages directories from interpreter search paths."""
    candidates: list[Path] = []

    for raw in sys.path:
        if not raw:
            continue
        path = Path(raw)
        if path.is_dir() and "site-packages" in path.parts:
            candidates.append(path)

    try:
        for raw in site.getsitepackages():
            path = Path(raw)
            if path.is_dir():
                candidates.append(path)
    except Exception:
        pass

    seen: set[str] = set()
    result: list[Path] = []
    for path in candidates:
        key = str(path.resolve())
        if key in seen:
            continue
        seen.add(key)
        result.append(path)
    return result


def _split_library_path(value: str) -> list[str]:
    """Split and normalize a PATH-like string."""
    if not value:
        return []
    return [part.strip() for part in value.split(os.pathsep) if part.strip()]


def _dedupe_library_paths(paths: list[str]) -> list[str]:
    """Keep only existing unique directories, preserving input order."""
    seen: set[str] = set()
    result: list[str] = []
    for raw in paths:
        path = Path(raw).expanduser()
        if not path.is_dir():
            continue
        normalized = str(path.resolve())
        if normalized in seen:
            continue
        seen.add(normalized)
        result.append(normalized)
    return result
