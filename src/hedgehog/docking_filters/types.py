"""Shared types, constants, and helper utilities for docking filters."""

from __future__ import annotations

import os
import shlex
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from hedgehog.configs.logger import logger

# Suppress common warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _project_root() -> Path:
    # src/hedgehog/docking_filters/types.py -> project root
    return Path(__file__).resolve().parent.parent.parent.parent


def _resolve_existing_path(base_folder: Path, path: str | Path) -> Path:
    """Resolve a path (possibly relative) to an existing absolute path."""
    p = Path(path)
    if p.is_absolute():
        return p

    candidates = [
        (base_folder / p).resolve(),
        (_project_root() / p).resolve(),
        (Path.cwd() / p).resolve(),
    ]
    for c in candidates:
        if c.exists():
            return c
    return (base_folder / p).resolve()


def _make_passthrough_df(n_mols: int, columns: dict[str, Any]) -> pd.DataFrame:
    """Build a 'pass everyone' DataFrame with mol_idx and the given constant columns.

    Args:
        n_mols: Number of molecules (rows).
        columns: Mapping of column name -> constant value for every row.

    Returns:
        DataFrame with ``mol_idx`` in ``range(n_mols)`` plus the requested columns.
    """
    data: dict[str, Any] = {"mol_idx": range(n_mols)}
    for col, val in columns.items():
        data[col] = [val] * n_mols if isinstance(val, float) and np.isnan(val) else val
    return pd.DataFrame(data)


# Sentinel for conformer-RMSD error results.
_CONFORMER_ERROR: dict[str, Any] = {
    "min_rmsd": float("inf"),
    "n_conformers_generated": 0,
    "passed": False,
}

# Sentinel for posebusters error results.
_POSEBUSTERS_ERROR: dict[str, Any] = {
    "no_clashes": False,
    "no_volume_clash": False,
    "not_too_far_away": False,
    "no_internal_clash": False,
    "passed": False,
}


def _error_result(**defaults: Any) -> dict[str, Any]:
    """Build a default error-result dict.

    Returns a dict with ``error: None`` merged with *defaults*.
    """
    out = dict(defaults)
    out.setdefault("error", None)
    return out


def _fail_result(error: str, **defaults: Any) -> dict[str, Any]:
    """Build a failure-result dict with a specific error message."""
    out = dict(defaults)
    out["error"] = error
    return out


def _aggregate_filter_results(
    raw_results: list[dict[str, Any]],
    col_mapping: dict[str, str],
    pass_col: str,
    filter_name: str,
) -> pd.DataFrame:
    """Aggregate per-molecule filter results into a DataFrame.

    Args:
        raw_results: List of dicts returned by per-molecule worker functions.
        col_mapping: Mapping of ``{result_key: dataframe_column}``.
        pass_col: Name of the boolean pass column in the DataFrame.
        filter_name: Human-readable filter name for log messages.

    Returns:
        DataFrame with ``mol_idx``, mapped columns, and ``pass_col``.
    """
    results: list[dict[str, Any]] = []
    for i, res in enumerate(raw_results):
        if res.get("error"):
            logger.warning("%s failed for mol %d: %s", filter_name, i, res["error"])
        row: dict[str, Any] = {"mol_idx": i}
        for src_key, dst_col in col_mapping.items():
            row[dst_col] = res[src_key]
        row[pass_col] = res["passed"]
        results.append(row)

    df = pd.DataFrame(results)
    logger.info(
        "%s: %d/%d passed",
        filter_name,
        int(df[pass_col].sum()),
        len(df),
    )
    return df


def _resolve_worker_command(
    env_var: str,
    venv_name: str,
    entry_name: str,
    module_path: str,
    install_hint: str,
) -> list[str]:
    """Resolve a worker command from environment, venv entry-point, or venv python.

    Resolution order:
    1. Environment variable (shell-split).
    2. Platform-specific entry-point script inside the venv.
    3. ``python -m <module_path>`` using the venv python.

    Args:
        env_var: Environment variable name (e.g. ``HEDGEHOG_SHEPHERD_WORKER_CMD``).
        venv_name: Venv directory name under project root (e.g. ``.venv-shepherd-worker``).
        entry_name: Entry-point script base name (e.g. ``hedgehog-shepherd-worker``).
        module_path: Dotted module path for ``python -m`` fallback.
        install_hint: Human-readable install instruction for the error message.

    Returns:
        Command tokens ready for ``subprocess.run()``.

    Raises:
        RuntimeError: If no command could be resolved.
    """
    raw_cmd = os.environ.get(env_var, "").strip()
    if raw_cmd:
        parts = shlex.split(raw_cmd)
        if not parts:
            raise RuntimeError(f"{env_var} is empty")
        return parts

    project_root = _project_root()
    if os.name == "nt":
        worker_entry = project_root / venv_name / "Scripts" / f"{entry_name}.exe"
        venv_python = project_root / venv_name / "Scripts" / "python.exe"
    else:
        worker_entry = project_root / venv_name / "bin" / entry_name
        venv_python = project_root / venv_name / "bin" / "python"

    if worker_entry.exists():
        return [str(worker_entry)]
    if venv_python.exists():
        return [str(venv_python), "-m", module_path]

    raise RuntimeError(
        f"No {entry_name} command found. Set {env_var} or run: {install_hint}"
    )
