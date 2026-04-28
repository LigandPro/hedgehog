"""Setup helpers for FSScore checkout and isolated worker runtime."""

from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass
from pathlib import Path

from hedgehog.configs.logger import logger
from hedgehog.setup._download import confirm_download, resolve_uv_binary

_FSSCORE_REPO_URL = "https://github.com/schwallergroup/fsscore.git"
_FSSCORE_SCORE_SCRIPT = Path("src") / "fsscore" / "score.py"
_FSSCORE_MODEL_PATH = Path("models") / "pretrain_graph_GGLGGL_ep242_best_valloss.ckpt"
_FSSCORE_OPTIONAL_ENV_ROOT_ENV = "HEDGEHOG_OPTIONAL_ENV_ROOT"
_FSSCORE_DEFAULT_ENV_DIR = ".venv-fsscore-worker"


@dataclass(frozen=True)
class FSScoreRuntime:
    """Resolved FSScore runtime paths for HEDGEHOG workers."""

    checkout_path: Path
    worker_python: Path
    model_path: Path
    env_path: Path


def _looks_like_fsscore_checkout(path: Path) -> bool:
    """Return True when required FSScore files are present in checkout."""
    return (path / _FSSCORE_SCORE_SCRIPT).exists() and (
        path / _FSSCORE_MODEL_PATH
    ).exists()


def _resolve_fsscore_env_path(project_root: Path) -> Path:
    """Resolve isolated FSScore worker env path."""
    optional_env_root = os.environ.get(_FSSCORE_OPTIONAL_ENV_ROOT_ENV, "").strip()
    if optional_env_root:
        return Path(optional_env_root).expanduser() / "fsscore"
    return project_root / _FSSCORE_DEFAULT_ENV_DIR


def _check_fsscore_runtime(worker_python: Path, model_path: Path) -> bool:
    """Return True when FSScore worker runtime is usable."""
    if not worker_python.exists() or not model_path.exists():
        return False

    probe_script = (
        "import importlib.util\n"
        "import pathlib\n"
        "import sys\n"
        "model_path = pathlib.Path(sys.argv[1])\n"
        "if not model_path.exists():\n"
        "    raise SystemExit(2)\n"
        "for module_name in ('fsscore', 'torch', 'pytorch_lightning'):\n"
        "    if importlib.util.find_spec(module_name) is None:\n"
        "        raise SystemExit(3)\n"
        "import pkg_resources\n"
        "print('ok')\n"
    )

    completed = subprocess.run(
        [str(worker_python), "-c", probe_script, str(model_path)],
        shell=False,
        capture_output=True,
        text=True,
        check=False,
    )
    return completed.returncode == 0 and completed.stdout.strip().endswith("ok")


def _install_fsscore_runtime(checkout_path: Path, env_path: Path) -> Path:
    """Create/refresh isolated FSScore env and install runtime dependencies."""
    env_path.parent.mkdir(parents=True, exist_ok=True)
    uv_bin = resolve_uv_binary()

    subprocess.run(
        [uv_bin, "venv", "--python", "3.10", str(env_path)],
        shell=False,
        check=True,
        timeout=900,
    )

    worker_python = env_path / "bin" / "python"
    subprocess.run(
        [
            uv_bin,
            "pip",
            "install",
            "--python",
            str(worker_python),
            "-e",
            str(checkout_path),
        ],
        shell=False,
        check=True,
        timeout=7200,
    )
    subprocess.run(
        [uv_bin, "pip", "install", "--python", str(worker_python), "setuptools<81"],
        shell=False,
        check=True,
        timeout=600,
    )
    return worker_python


def ensure_fsscore_checkout(project_root: Path) -> Path:
    """Ensure FSScore source checkout exists in ``modules/fsscore``."""
    checkout_path = project_root / "modules" / "fsscore"
    if _looks_like_fsscore_checkout(checkout_path):
        return checkout_path

    if checkout_path.exists():
        raise RuntimeError(
            "Existing modules/fsscore checkout is missing required files "
            "(src/fsscore/score.py or pretrained checkpoint)."
        )

    if not confirm_download("FSScore source checkout", "~500 MB (repo + checkpoint)"):
        raise RuntimeError("FSScore setup declined by user.")

    checkout_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Cloning FSScore repository to %s", checkout_path)
    subprocess.run(
        [
            "git",
            "clone",
            "--depth",
            "1",
            _FSSCORE_REPO_URL,
            str(checkout_path),
        ],
        shell=False,
        check=True,
        timeout=900,
    )

    if not _looks_like_fsscore_checkout(checkout_path):
        raise RuntimeError(
            "FSScore checkout completed but required files are missing "
            "(src/fsscore/score.py or pretrained checkpoint)."
        )

    logger.info("FSScore checkout is ready: %s", checkout_path)
    return checkout_path


def ensure_fsscore_runtime(project_root: Path) -> FSScoreRuntime:
    """Ensure FSScore checkout and isolated worker runtime are ready."""
    checkout_path = ensure_fsscore_checkout(project_root)
    model_path = checkout_path / _FSSCORE_MODEL_PATH
    env_path = _resolve_fsscore_env_path(project_root)
    worker_python = env_path / "bin" / "python"

    if _check_fsscore_runtime(worker_python, model_path):
        logger.info("FSScore runtime is ready: %s", worker_python)
        return FSScoreRuntime(
            checkout_path=checkout_path,
            worker_python=worker_python,
            model_path=model_path,
            env_path=env_path,
        )

    if not confirm_download(
        "FSScore isolated worker environment",
        "~7 GB (PyTorch + FSScore runtime dependencies)",
    ):
        raise RuntimeError("FSScore runtime setup declined by user.")

    logger.info("Installing FSScore runtime environment to %s", env_path)
    worker_python = _install_fsscore_runtime(checkout_path, env_path)

    if not _check_fsscore_runtime(worker_python, model_path):
        raise RuntimeError(
            "FSScore runtime installation completed but runtime probe failed "
            "(check torch/lightning/setuptools compatibility)."
        )

    logger.info("FSScore runtime is ready: %s", worker_python)
    return FSScoreRuntime(
        checkout_path=checkout_path,
        worker_python=worker_python,
        model_path=model_path,
        env_path=env_path,
    )
