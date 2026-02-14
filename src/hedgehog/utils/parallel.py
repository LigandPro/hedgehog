"""Unified parallelization utilities for all pipeline stages."""

from __future__ import annotations

import multiprocessing
import os
from collections.abc import Callable, Sequence
from typing import TypeVar

T = TypeVar("T")
R = TypeVar("R")


def resolve_n_jobs(
    stage_config: dict | None = None,
    global_config: dict | None = None,
    default: int = -1,
) -> int:
    """Resolve the number of parallel workers consistently across all stages.

    Priority:
      1. ``stage_config["n_jobs"]`` (per-stage override)
      2. ``global_config["n_jobs"]`` (master config)
      3. *default* argument

    Special values:
      - ``-1`` or ``0``: use all available cores (``SLURM_CPUS_PER_TASK``
        environment variable if set, otherwise ``os.cpu_count()``).

    Returns:
        Positive integer >= 1.
    """
    n_jobs: int | None = None

    if stage_config and "n_jobs" in stage_config:
        n_jobs = int(stage_config["n_jobs"])
    elif global_config and "n_jobs" in global_config:
        n_jobs = int(global_config["n_jobs"])
    else:
        n_jobs = default

    if n_jobs <= 0:
        slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
        if slurm_cpus:
            n_jobs = int(slurm_cpus)
        else:
            n_jobs = os.cpu_count() or 1

    return max(1, n_jobs)


def parallel_map(
    func: Callable,
    items: Sequence,
    n_jobs: int,
    chunksize: int | None = None,
    progress: Callable[[int, int], None] | None = None,
) -> list:
    """Apply *func* to every element of *items*, optionally in parallel.

    When ``n_jobs == 1`` (or only one item), processing is sequential — useful
    for debugging and environments where forking is problematic.

    Args:
        func: A **picklable** callable (top-level function, not a lambda).
        items: Sequence of arguments to map over.
        n_jobs: Number of worker processes.  ``1`` → sequential.
        chunksize: Items sent to each worker at a time.  ``None`` →
            auto-calculated as ``max(1, len(items) // (n_jobs * 4))``.

    Returns:
        List of results in the same order as *items*.
    """
    length = len(items)
    if length == 0:
        return []

    if n_jobs == 1 or length == 1:
        out = []
        for idx, item in enumerate(items, start=1):
            out.append(func(item))
            if progress:
                progress(idx, length)
        return out

    if chunksize is None:
        chunksize = max(1, length // (n_jobs * 4))

    ctx = multiprocessing.get_context("fork")
    with ctx.Pool(processes=n_jobs) as pool:
        out = []
        for idx, result in enumerate(
            pool.imap(func, items, chunksize=chunksize), start=1
        ):
            out.append(result)
            if progress:
                progress(idx, length)
        return out
