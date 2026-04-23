"""Auto-installation of AiZynthFinder retrosynthesis tool."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

from hedgehog.configs.logger import logger
from hedgehog.setup._download import confirm_download, resolve_uv_binary

_AIZYNTHFINDER_DIRNAME = "aizynthfinder"
_RETROSYNTHESIS_EXTRA = "retrosynthesis"


def _get_aizynthfinder_root(project_root: Path) -> Path:
    """Return the project-local AiZynthFinder workspace root."""
    return project_root / "modules" / _AIZYNTHFINDER_DIRNAME


def _get_config_path(project_root: Path) -> Path:
    """Return the project-local AiZynthFinder config path."""
    return _get_aizynthfinder_root(project_root) / "public" / "config.yml"


def _ensure_supported_python() -> None:
    """Fail fast when upstream AiZynthFinder does not support this interpreter."""
    version = sys.version_info[:3]
    if (3, 10) <= version < (3, 13):
        return

    version_str = ".".join(str(part) for part in version)
    raise RuntimeError(
        "AiZynthFinder upstream supports Python 3.10-3.12; "
        f"current interpreter is {version_str}."
    )


def _has_aizynthfinder_package(uv_bin: str, project_root: Path) -> bool:
    """Check whether the project environment already provides aizynthfinder."""
    probe = subprocess.run(
        [
            uv_bin,
            "run",
            "python",
            "-c",
            "import importlib.util, sys; "
            "sys.exit(0 if importlib.util.find_spec('aizynthfinder') else 1)",
        ],
        cwd=project_root,
        check=False,
        timeout=120,
    )
    return probe.returncode == 0


def _install_dependencies(
    uv_bin: str, project_root: Path, package_available: bool | None = None
) -> None:
    """Install the optional AiZynthFinder dependency into the project env."""
    if package_available is None:
        package_available = _has_aizynthfinder_package(uv_bin, project_root)
    if package_available:
        return

    logger.info(
        "Installing AiZynthFinder dependency into project environment "
        "(uv sync --extra %s)...",
        _RETROSYNTHESIS_EXTRA,
    )
    subprocess.run(
        [uv_bin, "sync", "--extra", _RETROSYNTHESIS_EXTRA],
        cwd=project_root,
        check=True,
        timeout=1800,
    )
    logger.info("AiZynthFinder dependency installed successfully")


def _ensure_public_data(
    uv_bin: str, project_root: Path, public_dir: Path, config_yml: Path
) -> None:
    """Download public data/models and ensure config.yml exists."""
    public_dir.mkdir(
        parents=True, exist_ok=True
    )  # Standard perms; public data directory

    def _download() -> None:
        logger.info("Downloading AiZynthFinder public data (models)...")
        subprocess.run(
            [
                uv_bin,
                "run",
                "python",
                "-m",
                "aizynthfinder.tools.download_public_data",
                str(public_dir),
            ],
            cwd=project_root,
            check=True,
            timeout=7200,
        )
        logger.info("Public data downloaded successfully")

    if not any(public_dir.iterdir()):
        _download()

    # Re-run if config.yml is missing (partial copy or interrupted download).
    if not config_yml.exists():
        logger.warning(
            "AiZynthFinder config.yml is missing; re-running public data setup."
        )
        _download()

    if not config_yml.exists():
        raise RuntimeError(
            f"AiZynthFinder config.yml not found after setup: {config_yml}"
        )


def _copy_logging_yml(project_root: Path, aizynth_dir: Path) -> None:
    """Copy logging.yml into the aizynthfinder data directory."""
    logging_candidates = [
        project_root / "src" / "hedgehog" / "synthesis" / "logging.yml",
        project_root / "src" / "hedgehog" / "stages" / "synthesis" / "logging.yml",
    ]
    logging_src = next((p for p in logging_candidates if p.exists()), None)
    if logging_src is None:
        return
    data_dir = aizynth_dir / "aizynthfinder" / "data"
    data_dir.mkdir(parents=True, exist_ok=True)  # Standard perms; public data directory
    shutil.copy2(logging_src, data_dir / "logging.yml")
    logger.info("Copied logging.yml to %s", data_dir)


def ensure_aizynthfinder(project_root: Path) -> Path:
    """Ensure AiZynthFinder is installed and return the config path.

    Installs the upstream ``aizynthfinder`` package into the project
    environment and downloads the public data into ``modules/aizynthfinder``.

    Steps:
      1. Return immediately if config already exists.
      2. Verify the current Python version is supported upstream.
      3. Verify ``uv`` is on PATH.
      4. Prompt the user to confirm the download.
      5. Install the optional ``retrosynthesis`` extra if needed.
      6. Download public data (models) if ``public/`` is empty.
      7. Copy ``logging.yml`` into the aizynthfinder data directory.

    Args:
        project_root: Absolute path to the hedgehog project root.

    Returns:
        Path to ``public/config.yml`` inside the AiZynthFinder tree.

    Raises:
        RuntimeError: If prerequisites are missing or the user declines.
    """
    aizynth_dir = _get_aizynthfinder_root(project_root)
    public_dir = aizynth_dir / "public"
    config_yml = _get_config_path(project_root)
    _ensure_supported_python()
    uv_bin = resolve_uv_binary()

    # 1. Already installed data
    if config_yml.exists():
        package_available = _has_aizynthfinder_package(uv_bin, project_root)
        if not package_available:
            if not confirm_download("AiZynthFinder", "package only"):
                raise RuntimeError("AiZynthFinder download declined by user.")
            _install_dependencies(uv_bin, project_root, package_available=False)
        _copy_logging_yml(project_root, aizynth_dir)
        logger.info("AiZynthFinder config found at %s", config_yml)
        return config_yml

    # 2. User confirmation
    if not confirm_download("AiZynthFinder", "800 MB (package + public data)"):
        raise RuntimeError("AiZynthFinder download declined by user.")

    aizynth_dir.mkdir(parents=True, exist_ok=True)

    # 3. Install dependencies
    _install_dependencies(uv_bin, project_root)

    # 4. Download public data
    _ensure_public_data(uv_bin, project_root, public_dir, config_yml)

    # 5. Copy logging.yml
    _copy_logging_yml(project_root, aizynth_dir)

    return config_yml
