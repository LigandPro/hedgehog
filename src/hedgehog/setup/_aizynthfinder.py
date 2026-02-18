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

    # 4. Ensure repository layout exists
    if not aizynth_dir.exists():
        if retro_dir.exists() and not (retro_dir / ".git").exists():
            logger.warning(
                "Found incomplete retrosynthesis directory at %s; recreating it.",
                retro_dir,
            )
            shutil.rmtree(retro_dir)

        if not retro_dir.exists():
            logger.info("Cloning retrosynthesis repository...")
            modules_dir = project_root / "modules"
            modules_dir.mkdir(parents=True, exist_ok=True)
            subprocess.run(
                ["git", "clone", "https://github.com/LigandPro/retrosynthesis.git"],
                cwd=modules_dir,
                check=True,
                timeout=1800,
            )
            logger.info("Repository cloned successfully")
        elif (retro_dir / ".git").exists():
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

    if not aizynth_dir.exists():
        raise RuntimeError(
            f"AiZynthFinder directory not found after setup: {aizynth_dir}"
        )

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
    def _download_public_data() -> None:
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

    if not any(public_dir.iterdir()):
        _download_public_data()

    # The public data download is also responsible for writing config.yml.
    # If the directory already had content but config.yml is missing (partial
    # manual copy, interrupted download), try running it again.
    if not config_yml.exists():
        logger.warning("AiZynthFinder config.yml is missing; re-running public data setup.")
        _download_public_data()

    if not config_yml.exists():
        raise RuntimeError(f"AiZynthFinder config.yml not found after setup: {config_yml}")

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
