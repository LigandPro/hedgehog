"""Auto-installation of AiZynthFinder retrosynthesis tool."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from hedgehog.configs.logger import logger
from hedgehog.setup._download import confirm_download


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
        logger.info("AiZynthFinder config found at %s", config_yml)
        return config_yml

    # 2. Prerequisites
    if not shutil.which("git"):
        raise RuntimeError(
            "git is not installed. Please install git to set up AiZynthFinder."
        )
    if not shutil.which("uv"):
        raise RuntimeError(
            "uv is not installed. Please install uv to set up AiZynthFinder."
        )

    # 3. User confirmation
    if not confirm_download("AiZynthFinder", "~800 MB (repo + models)"):
        raise RuntimeError("AiZynthFinder download declined by user.")

    # 4. Clone repository
    if not retro_dir.exists():
        logger.info("Cloning retrosynthesis repository...")
        modules_dir = project_root / "modules"
        modules_dir.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            ["git", "clone", "https://github.com/LigandPro/retrosynthesis.git"],
            cwd=modules_dir,
            check=True,
            timeout=600,
        )
        logger.info("Repository cloned successfully")

    # 5. Install dependencies
    if not (aizynth_dir / ".venv").exists():
        logger.info("Installing AiZynthFinder dependencies (uv sync)...")
        subprocess.run(
            ["uv", "sync"],
            cwd=aizynth_dir,
            check=True,
            timeout=600,
        )
        logger.info("Dependencies installed successfully")

    # 6. Download public data
    public_dir.mkdir(parents=True, exist_ok=True)
    if not any(public_dir.iterdir()):
        logger.info("Downloading AiZynthFinder public data (models)...")
        subprocess.run(
            [
                "uv",
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

    # 7. Copy logging.yml
    logging_src = (
        project_root / "src" / "hedgehog" / "stages" / "synthesis" / "logging.yml"
    )
    data_dir = aizynth_dir / "aizynthfinder" / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    if logging_src.exists():
        shutil.copy2(logging_src, data_dir / "logging.yml")
        logger.info("Copied logging.yml to %s", data_dir)

    return config_yml
