"""Set up isolated virtual environment for Shepherd-Score worker."""

from __future__ import annotations

import shlex
import subprocess
from pathlib import Path

from hedgehog.setup._common import (
    resolve_python_binary as _resolve_python_binary,
)
from hedgehog.setup._common import (
    run as _run,
)
from hedgehog.setup._common import (
    venv_python as _venv_python,
)
from hedgehog.setup._common import (
    venv_worker_entry as _venv_worker_entry_raw,
)
from hedgehog.setup._common import (
    verify_worker as _verify_worker_raw,
)
from hedgehog.setup._download import confirm_download, resolve_uv_binary

_WORKER_SCRIPT_NAME = "hedgehog-shepherd-worker"
_FALLBACK_MODULE = "hedgehog.workers.shepherd_worker"


def _venv_worker_entry(venv_dir: Path) -> Path:
    return _venv_worker_entry_raw(venv_dir, _WORKER_SCRIPT_NAME)


def _verify_worker(worker_entry: Path, venv_python: Path, cwd: Path) -> None:
    _verify_worker_raw(worker_entry, venv_python, cwd, _FALLBACK_MODULE)


def ensure_shepherd_worker(project_root: Path, python_bin: str | None = None) -> Path:
    """Ensure shepherd worker virtualenv exists and return worker entry path."""
    uv_bin = resolve_uv_binary()

    selected_python = _resolve_python_binary(python_bin)
    venv_dir = project_root / ".venv-shepherd-worker"
    venv_python = _venv_python(venv_dir)
    worker_entry = _venv_worker_entry(venv_dir)

    if venv_python.exists():
        try:
            _verify_worker(worker_entry, venv_python, project_root)
            return worker_entry if worker_entry.exists() else venv_python
        except (OSError, subprocess.SubprocessError, RuntimeError):
            pass

    if not confirm_download(
        "Shepherd worker dependencies", "~1 GB (PyTorch/Open3D stack)"
    ):
        raise RuntimeError("Shepherd worker setup declined by user.")

    _run([selected_python, "-m", "venv", str(venv_dir)], cwd=project_root, timeout=600)

    if not venv_python.exists():
        raise RuntimeError(f"Failed to create virtualenv at {venv_dir}")

    install_cmd = [
        uv_bin,
        "pip",
        "install",
        "--python",
        str(venv_python),
        "-e",
        ".[shepherd]",
    ]
    _run(install_cmd, cwd=project_root, timeout=3600)

    try:
        _verify_worker(worker_entry, venv_python, project_root)
    except Exception as exc:  # noqa: BLE001 — intentional: Shepherd verification may raise diverse errors
        raise RuntimeError(
            "Shepherd worker was installed but failed verification. "
            f"Command: {shlex.join([str(venv_python), '-m', _FALLBACK_MODULE, '--help'])}. "
            f"Error: {exc}"
        ) from exc

    return worker_entry if worker_entry.exists() else venv_python
