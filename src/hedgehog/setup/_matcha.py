"""Clone and update the official Matcha repository for hedgehog docking."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from hedgehog.configs.logger import logger
from hedgehog.setup._download import confirm_download

_MATCHA_REPO_URL = "https://github.com/LigandPro/Matcha.git"
_DEFAULT_CHECKOUT_DIR = "modules/matcha_remote"


def _run_git(cmd: list[str], cwd: Path) -> None:
    subprocess.run(cmd, cwd=cwd, check=True, timeout=1800)


def ensure_matcha_checkout(
    project_root: Path,
    *,
    repo_url: str = _MATCHA_REPO_URL,
    checkout_dir: str | Path | None = None,
    pr_number: int | None = None,  # noqa: ARG001 — kept for config compatibility
) -> Path:
    """Clone LigandPro/Matcha (or update it) and return the checkout path."""
    del pr_number  # only origin/main is supported

    target_dir = (
        Path(checkout_dir).expanduser()
        if checkout_dir is not None
        else (project_root / _DEFAULT_CHECKOUT_DIR)
    )
    if not target_dir.is_absolute():
        target_dir = (project_root / target_dir).resolve()

    if not shutil.which("git"):
        raise RuntimeError("git is required to fetch Matcha from GitHub.")

    if not target_dir.exists() and not confirm_download(
        "Matcha checkout",
        "~repository checkout from LigandPro/Matcha plus uv environment",
    ):
        raise RuntimeError("Matcha download declined by user.")

    if not target_dir.exists() or not (target_dir / ".git").is_dir():
        if target_dir.exists():
            logger.warning(
                "Removing incomplete Matcha checkout at %s before recloning.",
                target_dir,
            )
            shutil.rmtree(target_dir)
        target_dir.parent.mkdir(parents=True, exist_ok=True)
        logger.info("Cloning Matcha from %s into %s", repo_url, target_dir)
        _run_git(["git", "clone", repo_url, str(target_dir)], cwd=target_dir.parent)
    else:
        logger.info("Updating Matcha checkout at %s", target_dir)
        _run_git(["git", "fetch", "origin", "main"], cwd=target_dir)

    _run_git(["git", "checkout", "--detach", "origin/main"], cwd=target_dir)

    if not (target_dir / "matcha" / "cli.py").is_file():
        raise RuntimeError(
            f"Matcha checkout is missing matcha/cli.py after cloning {repo_url}."
        )

    return target_dir
