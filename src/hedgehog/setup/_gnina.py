"""Auto-installation of GNINA docking binary."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

from hedgehog.setup._download import confirm_download, download_with_progress

_GNINA_CACHE_DIR = Path.home() / ".hedgehog" / "bin"


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
    if not confirm_download("GNINA", "~200 MB"):
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
    """Query GitHub releases API for the latest GNINA download URL.

    Returns the ``browser_download_url`` for the first release asset whose
    name does **not** contain "cuda" (case-insensitive).

    Raises:
        RuntimeError: If the API call fails or no suitable asset is found.
    """
    api_url = "https://api.github.com/repos/gnina/gnina/releases/latest"
    req = urllib.request.Request(api_url, headers={"Accept": "application/json"})

    try:
        with urllib.request.urlopen(req, timeout=30) as resp:  # noqa: S310
            data = json.loads(resp.read())
    except Exception as exc:
        raise RuntimeError(
            f"Failed to query GNINA releases from GitHub: {exc}"
        ) from exc

    for asset in data.get("assets", []):
        name = asset.get("name", "")
        if "cuda" not in name.lower():
            return asset["browser_download_url"]

    raise RuntimeError(
        "No suitable GNINA release asset found on GitHub (all assets appear CUDA-specific)."
    )


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
