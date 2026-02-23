"""Shared helper utilities for isolated worker virtual-environment setup."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path


def venv_python(venv_dir: Path) -> Path:
    """Return the path to the Python interpreter inside a virtual environment."""
    if os.name == "nt":
        return venv_dir / "Scripts" / "python.exe"
    return venv_dir / "bin" / "python"


def venv_worker_entry(venv_dir: Path, script_name: str) -> Path:
    """Return the path to a console-script entry point inside a virtual environment.

    Args:
        venv_dir: Root of the virtual environment.
        script_name: Base name of the entry-point script (without ``.exe`` suffix).
    """
    if os.name == "nt":
        return venv_dir / "Scripts" / f"{script_name}.exe"
    return venv_dir / "bin" / script_name


def run(cmd: list[str], cwd: Path, timeout: int = 1800) -> None:
    """Run a subprocess command with a timeout."""
    subprocess.run(cmd, cwd=cwd, check=True, timeout=timeout)


def resolve_python_binary(python_bin: str | None) -> str:
    """Resolve a Python interpreter path, falling back to well-known candidates.

    Args:
        python_bin: Explicit interpreter path/name requested by the user, or *None*
            to auto-detect.

    Raises:
        RuntimeError: If no suitable interpreter is found.
    """
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


def verify_worker(
    worker_entry: Path,
    venv_python_path: Path,
    cwd: Path,
    fallback_module: str,
) -> None:
    """Verify that a worker entry-point is functional.

    Tries the console-script first; falls back to ``python -m <fallback_module>``.

    Args:
        worker_entry: Path to the console-script entry-point.
        venv_python_path: Path to the venv Python interpreter.
        cwd: Working directory for the subprocess.
        fallback_module: Dotted module path used as ``-m`` fallback.
    """
    if worker_entry.exists():
        run([str(worker_entry), "--help"], cwd=cwd, timeout=60)
        return
    run(
        [str(venv_python_path), "-m", fallback_module, "--help"],
        cwd=cwd,
        timeout=60,
    )
