from __future__ import annotations

import importlib
import logging
import os
from collections.abc import Callable
from functools import lru_cache, wraps
from typing import Any

_SUPPRESS_ENV = "HEDGEHOG_SHOW_PANDASTOOLS_WARNINGS"
_PANDASTOOLS_LOGGERS = (
    "rdkit.Chem.PandasTools",
    "rdkit.Chem.PandasPatcher",
)
_PARALLEL_PATCH_FLAG = "__hedgehog_parallel_patch_installed"
_WRAPPED_FN_FLAG = "__hedgehog_parallel_wrapper_installed"


def _should_suppress_pandastools_warning() -> bool:
    """Return True when PandasTools warning suppression is enabled."""
    return os.environ.get(_SUPPRESS_ENV, "").strip() != "1"


def suppress_pandastools_warning() -> None:
    """Silence noisy RDKit PandasTools warnings in this process."""
    if not _should_suppress_pandastools_warning():
        return

    for logger_name in _PANDASTOOLS_LOGGERS:
        ext_logger = logging.getLogger(logger_name)
        ext_logger.setLevel(logging.ERROR)
        if not any(
            isinstance(handler, logging.NullHandler) for handler in ext_logger.handlers
        ):
            ext_logger.addHandler(logging.NullHandler())


def _parallel_worker_initializer(
    user_initializer: Callable[..., Any] | None,
    user_initargs: tuple[Any, ...],
) -> None:
    """Apply warning suppression in every joblib worker process."""
    suppress_pandastools_warning()
    if user_initializer is not None:
        user_initializer(*user_initargs)


def _wrap_parallel_runner(func: Callable[..., Any]) -> Callable[..., Any]:
    """Wrap datamol parallel helpers to inject worker initializer."""
    if getattr(func, _WRAPPED_FN_FLAG, False):
        return func

    @wraps(func)
    def wrapped(*args: Any, **kwargs: Any) -> Any:
        if not _should_suppress_pandastools_warning():
            return func(*args, **kwargs)

        user_initializer = kwargs.pop("initializer", None)
        user_initargs = kwargs.pop("initargs", ())
        if user_initargs is None:
            user_initargs = ()
        elif not isinstance(user_initargs, tuple):
            user_initargs = tuple(user_initargs)

        kwargs["initializer"] = _parallel_worker_initializer
        kwargs["initargs"] = (user_initializer, user_initargs)
        return func(*args, **kwargs)

    setattr(wrapped, _WRAPPED_FN_FLAG, True)
    return wrapped


def _patch_datamol_parallel_helpers(dm_module: Any) -> None:
    """Patch datamol parallel helpers once per process."""
    if getattr(dm_module, _PARALLEL_PATCH_FLAG, False):
        return

    for fn_name in ("parallelized", "parallelized_with_batches"):
        current = getattr(dm_module, fn_name, None)
        if callable(current):
            setattr(dm_module, fn_name, _wrap_parallel_runner(current))

    setattr(dm_module, _PARALLEL_PATCH_FLAG, True)


@lru_cache(maxsize=1)
def import_datamol_quietly() -> Any:
    """Import datamol while muting native stdout/stderr noise."""
    suppress_pandastools_warning()

    stdout_fd = os.dup(1)
    stderr_fd = os.dup(2)
    devnull_fd = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull_fd, 1)
        os.dup2(devnull_fd, 2)
        dm_module = importlib.import_module("datamol")
        _patch_datamol_parallel_helpers(dm_module)
        # Force-load convert once while warning suppression is active.
        importlib.import_module("datamol.convert")
        return dm_module
    finally:
        os.dup2(stdout_fd, 1)
        os.dup2(stderr_fd, 2)
        os.close(devnull_fd)
        os.close(stdout_fd)
        os.close(stderr_fd)
