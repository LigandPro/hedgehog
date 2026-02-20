"""Set up isolated virtual environment for nvMolKit conformer worker."""

from __future__ import annotations

import os
import shlex
import shutil
import subprocess
from pathlib import Path

from hedgehog.setup._download import confirm_download, resolve_uv_binary

_MATCHA_REPO_URL = "https://github.com/LigandPro/Matcha.git"


def _venv_python(venv_dir: Path) -> Path:
    if os.name == "nt":
        return venv_dir / "Scripts" / "python.exe"
    return venv_dir / "bin" / "python"


def _venv_worker_entry(venv_dir: Path) -> Path:
    if os.name == "nt":
        return venv_dir / "Scripts" / "matcha-nvmolkit-worker.exe"
    return venv_dir / "bin" / "matcha-nvmolkit-worker"


def _run(cmd: list[str], cwd: Path, timeout: int = 1800) -> None:
    subprocess.run(cmd, cwd=cwd, check=True, timeout=timeout)


def _resolve_python_binary(python_bin: str | None) -> str:
    if python_bin:
        explicit = Path(python_bin)
        if explicit.exists():
            return str(explicit)
        found = shutil.which(python_bin)
        if found:
            return found
        raise RuntimeError(f"Requested Python interpreter was not found: {python_bin}")

    for candidate in ("python3.12", "python3.11", "python3.10"):
        found = shutil.which(candidate)
        if found:
            return found

    raise RuntimeError(
        "No supported Python interpreter found. Install one of: "
        "python3.12, python3.11, python3.10"
    )


def _verify_worker(worker_entry: Path, venv_python: Path, cwd: Path) -> None:
    if worker_entry.exists():
        _run([str(worker_entry), "--help"], cwd=cwd, timeout=60)
        return
    _run(
        [str(venv_python), "-m", "matcha_nvmolkit_worker.cli", "--help"],
        cwd=cwd,
        timeout=60,
    )


def _ensure_matcha_setup_script(project_root: Path) -> Path:
    cache_root = project_root / ".cache" / "nvmolkit-worker"
    repo_dir = cache_root / "Matcha"
    script_path = repo_dir / "scripts" / "setup_nvmolkit_worker.py"
    cache_root.mkdir(parents=True, exist_ok=True)

    if not repo_dir.exists():
        _run(
            ["git", "clone", _MATCHA_REPO_URL, str(repo_dir)],
            cwd=cache_root,
            timeout=1800,
        )
    else:
        try:
            _run(["git", "fetch", "--all", "--tags"], cwd=repo_dir, timeout=1800)
            _run(["git", "pull", "--ff-only"], cwd=repo_dir, timeout=1800)
        except (OSError, subprocess.SubprocessError):
            shutil.rmtree(repo_dir, ignore_errors=True)
            _run(
                ["git", "clone", _MATCHA_REPO_URL, str(repo_dir)],
                cwd=cache_root,
                timeout=1800,
            )

    if not script_path.exists():
        raise RuntimeError(f"setup_nvmolkit_worker.py not found in {repo_dir}")
    return script_path


def ensure_nvmolkit_worker(project_root: Path, python_bin: str | None = None) -> Path:
    """Ensure nvmolkit worker virtualenv exists and return worker entry path."""
    resolve_uv_binary()  # Pre-flight check for consistent error semantics.
    selected_python = _resolve_python_binary(python_bin)
    venv_dir = project_root / ".venv-nvmolkit-worker"
    venv_python = _venv_python(venv_dir)
    worker_entry = _venv_worker_entry(venv_dir)

    if venv_python.exists():
        try:
            _verify_worker(worker_entry, venv_python, project_root)
            return worker_entry if worker_entry.exists() else venv_python
        except (OSError, subprocess.SubprocessError, RuntimeError):
            pass

    if not confirm_download(
        "nvMolKit worker dependencies", "~2 GB (CUDA + nvMolKit build stack)"
    ):
        raise RuntimeError("nvMolKit worker setup declined by user.")

    script_path = _ensure_matcha_setup_script(project_root)

    cmd = [
        selected_python,
        str(script_path),
        "--venv",
        str(venv_dir),
    ]
    _run(cmd, cwd=script_path.parent.parent, timeout=7200)

    try:
        _verify_worker(worker_entry, venv_python, project_root)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(
            "nvMolKit worker was installed but failed verification. "
            f"Command: {shlex.join([str(venv_python), '-m', 'matcha_nvmolkit_worker.cli', '--help'])}. "
            f"Error: {exc}"
        ) from exc

    return worker_entry if worker_entry.exists() else venv_python
