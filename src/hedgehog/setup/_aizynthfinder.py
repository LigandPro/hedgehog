"""Auto-installation of AiZynthFinder retrosynthesis tool."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from hedgehog.configs.logger import logger
from hedgehog.setup._download import confirm_download, resolve_uv_binary


def _ensure_retro_repo(project_root: Path, retro_dir: Path) -> None:
    """Clone or update the retrosynthesis repository."""
    if retro_dir.exists() and not (retro_dir / ".git").exists():
        logger.warning(
            "Found incomplete retrosynthesis directory at %s; recreating it.",
            retro_dir,
        )
        shutil.rmtree(retro_dir)

    if not retro_dir.exists():
        logger.info("Cloning retrosynthesis repository...")
        modules_dir = project_root / "modules"
        modules_dir.mkdir(
            parents=True, exist_ok=True
        )  # Standard perms; project data directory
        subprocess.run(
            ["git", "clone", "https://github.com/LigandPro/retrosynthesis.git"],
            cwd=modules_dir,
            check=True,
            timeout=1800,
        )
        logger.info("Repository cloned successfully")
        return

    if (retro_dir / ".git").exists():
        logger.info("Updating retrosynthesis repository...")
        try:
            subprocess.run(
                ["git", "pull", "--ff-only"],
                cwd=retro_dir,
                check=True,
                timeout=1800,
            )
        except subprocess.TimeoutExpired:
            logger.warning(
                "git pull timed out after 1800s; continuing with existing checkout."
            )


def _install_dependencies(uv_bin: str, aizynth_dir: Path) -> None:
    """Install AiZynthFinder Python dependencies via uv sync."""
    if (aizynth_dir / ".venv").exists():
        return
    logger.info("Installing AiZynthFinder dependencies (uv sync)...")
    subprocess.run(
        [uv_bin, "sync"],
        cwd=aizynth_dir,
        check=True,
        timeout=600,
    )
    logger.info("Dependencies installed successfully")


def _ensure_public_data(
    uv_bin: str, aizynth_dir: Path, public_dir: Path, config_yml: Path
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
            cwd=aizynth_dir,
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

    Replicates the logic of ``modules/install_aizynthfinder.sh`` as a
    Python function so the synthesis stage can auto-install when needed.

    Steps:
      1. Return immediately if config already exists.
      2. Verify ``git`` and ``uv`` are on PATH.
      3. Prompt the user to confirm the download.
      4. Clone the retrosynthesis repo if missing.
      5. Install Python dependencies via ``uv sync`` if ``.venv`` is missing.
      6. Download public data (models) if ``public/`` is empty.
      7. Copy ``logging.yml`` into the aizynthfinder data directory.

    Args:
        project_root: Absolute path to the hedgehog project root.

    Returns:
        Path to ``public/config.yml`` inside the AiZynthFinder tree.

    Raises:
        RuntimeError: If prerequisites are missing or the user declines.
    """
    retro_dir = project_root / "modules" / "retrosynthesis"
    aizynth_dir = retro_dir / "aizynthfinder"
    public_dir = aizynth_dir / "public"
    config_yml = public_dir / "config.yml"

    # 1. Already installed
    if config_yml.exists():
        _copy_logging_yml(project_root, aizynth_dir)
        logger.info("AiZynthFinder config found at %s", config_yml)
        return config_yml

    # 2. Prerequisites
    if not shutil.which("git"):
        raise RuntimeError(
            "git is not installed. Please install git to set up AiZynthFinder."
        )
    uv_bin = resolve_uv_binary()

    # 3. User confirmation
    if not confirm_download("AiZynthFinder", "~800 MB (repo + models)"):
        raise RuntimeError("AiZynthFinder download declined by user.")

    # 4. Ensure repository layout exists
    if not aizynth_dir.exists():
        _ensure_retro_repo(project_root, retro_dir)

    if not aizynth_dir.exists():
        raise RuntimeError(
            f"AiZynthFinder directory not found after setup: {aizynth_dir}"
        )

    # 5. Install dependencies
    _install_dependencies(uv_bin, aizynth_dir)

    # 6. Download public data
    _ensure_public_data(uv_bin, aizynth_dir, public_dir, config_yml)

    # 7. Copy logging.yml
    _copy_logging_yml(project_root, aizynth_dir)

    return config_yml
