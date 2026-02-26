"""Set up isolated virtual environment for nvMolKit worker."""

from __future__ import annotations

import os
import shlex
import shutil
import subprocess
from pathlib import Path

from hedgehog.setup._download import confirm_download, resolve_uv_binary


def _venv_python(venv_dir: Path) -> Path:
    if os.name == "nt":
        return venv_dir / "Scripts" / "python.exe"
    return venv_dir / "bin" / "python"


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
        "No supported Python interpreter found. Install one of: python3.12, python3.11, python3.10"
    )


def _verify_worker(venv_python: Path, cwd: Path) -> None:
    _run(
        [str(venv_python), "-c", "import nvmolkit"],
        cwd=cwd,
        timeout=60,
    )


def ensure_nvmolkit_worker(project_root: Path, python_bin: str | None = None) -> Path:
    """Ensure nvMolKit worker virtualenv exists and return worker entry path."""
    uv_bin = resolve_uv_binary()

    selected_python = _resolve_python_binary(python_bin)
    venv_dir = project_root / ".venv-nvmolkit-worker"
    venv_python = _venv_python(venv_dir)

    if venv_python.exists():
        try:
            _verify_worker(venv_python, project_root)
            return venv_python
        except Exception:
            pass

    if not confirm_download("nvMolKit worker dependencies", "~1 GB (nvMolKit stack)"):
        raise RuntimeError("nvMolKit worker setup declined by user.")

    _run([selected_python, "-m", "venv", str(venv_dir)], cwd=project_root, timeout=600)

    if not venv_python.exists():
        raise RuntimeError(f"Failed to create virtualenv at {venv_dir}")

    install_cmd = [
        uv_bin,
        "pip",
        "install",
        "--python",
        str(venv_python),
        "nvmolkit",
    ]
    _run(install_cmd, cwd=project_root, timeout=3600)

    try:
        _verify_worker(venv_python, project_root)
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(
            "nvMolKit worker was installed but failed verification. "
            f"Command: {shlex.join([str(venv_python), '-c', 'import nvmolkit'])}. "
            f"Error: {exc}"
        ) from exc

    return venv_python
