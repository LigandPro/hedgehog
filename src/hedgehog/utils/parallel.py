"""Unified parallelization utilities for all pipeline stages."""

from __future__ import annotations

import multiprocessing
import os
import sys
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
    stage_n_jobs: int | None = None
    global_n_jobs: int | None = None

    if stage_config and "n_jobs" in stage_config:
        stage_n_jobs = int(stage_config["n_jobs"])
    if global_config and "n_jobs" in global_config:
        global_n_jobs = int(global_config["n_jobs"])

    if stage_n_jobs is not None:
        if stage_n_jobs > 0:
            n_jobs = stage_n_jobs
        elif global_n_jobs is not None and global_n_jobs > 0:
            n_jobs = global_n_jobs
        else:
            n_jobs = default
    elif global_n_jobs is not None:
        n_jobs = global_n_jobs
    else:
        n_jobs = default

    if n_jobs <= 0:
        env_n_jobs = os.environ.get("MOLSCORE_NJOBS")
        if env_n_jobs:
            try:
                parsed = int(env_n_jobs)
                if parsed > 0:
                    n_jobs = parsed
            except ValueError:
                pass

    if n_jobs <= 0:
        slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
        if slurm_cpus:
            n_jobs = int(slurm_cpus)
        else:
            n_jobs = os.cpu_count() or 1

    return max(1, n_jobs)


def _get_mp_context() -> multiprocessing.context.BaseContext:
    """Select a safe multiprocessing start method.

    Notes:
    - ``fork`` is fast but can deadlock in multi-threaded parents (e.g., Rich
      progress rendering). Prefer ``forkserver`` on Linux when available.
    - ``spawn`` is the most portable fallback.
    - Override via ``HEDGEHOG_MP_START_METHOD`` if needed.
    """
    forced = os.environ.get("HEDGEHOG_MP_START_METHOD")
    if forced:
        return multiprocessing.get_context(forced)

    if sys.platform.startswith("linux"):
        try:
            return multiprocessing.get_context("forkserver")
        except ValueError:
            return multiprocessing.get_context("spawn")

    return multiprocessing.get_context("spawn")


def parallel_map(
    func: Callable,
    items: Sequence,
    n_jobs: int,
    chunksize: int | None = None,
    progress: Callable[[int, int], None] | None = None,
    initializer: Callable[..., None] | None = None,
    initargs: tuple = (),
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
        initializer: Optional initializer called once per worker (and once in
            sequential mode) to set up per-process globals.
        initargs: Arguments passed to *initializer*.

    Returns:
        List of results in the same order as *items*.
    """
    length = len(items)
    if length == 0:
        return []

    if n_jobs == 1 or length == 1:
        if initializer:
            initializer(*initargs)
        out = []
        for idx, item in enumerate(items, start=1):
            out.append(func(item))
            if progress:
                progress(idx, length)
        return out

    if chunksize is None:
        chunksize = max(1, length // (n_jobs * 4))

    ctx = _get_mp_context()
    with ctx.Pool(
        processes=n_jobs,
        initializer=initializer,
        initargs=initargs,
    ) as pool:
        out = []
        for idx, result in enumerate(
            pool.imap(func, items, chunksize=chunksize), start=1
        ):
            out.append(result)
            if progress:
                progress(idx, length)
        return out
