"""Setup helper for optional GASA scorer checkout and isolated worker env."""

from __future__ import annotations

import os
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from hedgehog.setup._download import confirm_download, resolve_uv_binary

_GASA_REPO_URL = "https://github.com/cadd-synthetic/GASA.git"
_GASA_REQUIRED_FILES = (
    Path("gasa.py"),
    Path("model") / "gasa.pth",
    Path("model") / "gasa.json",
)
_GASA_DEPS = (
    "setuptools<81",
    "wheel",
    "numpy==1.26.4",
    "pandas==2.2.3",
    "scikit-learn==1.5.2",
    "hyperopt==0.2.7",
    "packaging==25.0",
    "rdkit-pypi==2022.9.5",
    "torch==2.2.2",
    "dgl==1.1.3",
    "dgllife==0.3.2",
)


@dataclass(frozen=True)
class GasaSetupResult:
    """Resolved checkout and worker Python path for GASA runtime."""

    repo_path: Path
    worker_python: Path


def _run(cmd: list[str], cwd: Path, timeout: int = 1800) -> None:
    subprocess.run(cmd, cwd=cwd, check=True, timeout=timeout)


def _venv_python(venv_dir: Path) -> Path:
    if os.name == "nt":
        return venv_dir / "Scripts" / "python.exe"
    return venv_dir / "bin" / "python"


def _resolve_python_binary(python_bin: str | None) -> str:
    if python_bin:
        explicit = Path(python_bin).expanduser()
        if explicit.exists():
            return str(explicit)
        found = shutil.which(python_bin)
        if found:
            return found
        if any(sep in python_bin for sep in ("/", "\\")):
            raise RuntimeError(
                f"Requested Python interpreter was not found: {python_bin}"
            )
        return python_bin

    for candidate in ("python3.10", "python3.11"):
        found = shutil.which(candidate)
        if found:
            return found

    # Allow uv to provision a compatible interpreter when none is preinstalled.
    return "3.10"


def _looks_like_gasa_checkout(path: Path) -> bool:
    return all((path / rel).exists() for rel in _GASA_REQUIRED_FILES)


def _ensure_gasa_checkout(project_root: Path) -> Path:
    checkout_path = project_root / "modules" / "gasa"
    if _looks_like_gasa_checkout(checkout_path):
        return checkout_path

    if checkout_path.exists():
        raise RuntimeError(
            "Existing modules/gasa checkout is missing required files "
            "(gasa.py, model/gasa.pth, or model/gasa.json)."
        )

    if not confirm_download("GASA source checkout", "~20 MB (repo + pretrained model)"):
        raise RuntimeError("GASA setup declined by user.")

    checkout_path.parent.mkdir(parents=True, exist_ok=True)
    _run(
        [
            "git",
            "clone",
            "--depth",
            "1",
            _GASA_REPO_URL,
            str(checkout_path),
        ],
        cwd=project_root,
        timeout=900,
    )

    if not _looks_like_gasa_checkout(checkout_path):
        raise RuntimeError(
            "GASA checkout completed but required files are missing "
            "(gasa.py, model/gasa.pth, or model/gasa.json)."
        )
    return checkout_path


def _resolve_worker_venv_dir(project_root: Path) -> Path:
    optional_env_root = os.environ.get("HEDGEHOG_OPTIONAL_ENV_ROOT")
    if optional_env_root:
        return Path(optional_env_root).expanduser() / "gasa"
    return project_root / ".venv-gasa-worker"


def _verify_gasa_dependencies(venv_python: Path) -> None:
    _run(
        [
            str(venv_python),
            "-c",
            "import numpy, pandas, sklearn, hyperopt, rdkit, torch, dgl, dgllife",
        ],
        cwd=venv_python.parent.parent,
        timeout=120,
    )


def ensure_gasa_worker(
    project_root: Path, python_bin: str | None = None
) -> GasaSetupResult:
    """Ensure GASA checkout/env exist and return worker execution paths."""
    repo_path = _ensure_gasa_checkout(project_root)
    venv_dir = _resolve_worker_venv_dir(project_root)
    venv_python = _venv_python(venv_dir)
    selected_python = _resolve_python_binary(python_bin)

    if venv_python.exists():
        try:
            _verify_gasa_dependencies(venv_python)
        except Exception:
            # Recreate stale/broken env below.
            pass
        else:
            return GasaSetupResult(repo_path=repo_path, worker_python=venv_python)

    if not confirm_download(
        "GASA worker dependencies",
        "~4 GB (PyTorch + DGL + RDKit stack)",
    ):
        raise RuntimeError("GASA worker setup declined by user.")

    uv_bin = resolve_uv_binary()
    if venv_dir.exists():
        if venv_dir.is_dir():
            shutil.rmtree(venv_dir)
        else:
            venv_dir.unlink()
    _run(
        [uv_bin, "venv", "--python", selected_python, str(venv_dir)],
        cwd=project_root,
        timeout=600,
    )

    if not venv_python.exists():
        raise RuntimeError(f"Failed to create virtualenv at {venv_dir}")

    _run(
        [
            uv_bin,
            "pip",
            "install",
            "--python",
            str(venv_python),
            *_GASA_DEPS,
        ],
        cwd=repo_path,
        timeout=3600,
    )

    try:
        _verify_gasa_dependencies(venv_python)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(
            "GASA worker dependencies were installed but verification failed. "
            f"Command: {shlex.join([str(venv_python), '-c', 'import torch, dgl, dgllife'])}. "
            f"Error: {exc}"
        ) from exc

    return GasaSetupResult(repo_path=repo_path, worker_python=venv_python)
