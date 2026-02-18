from __future__ import annotations

import importlib
import os
from functools import lru_cache
from typing import Any


@lru_cache(maxsize=1)
def import_datamol_quietly() -> Any:
    """Import datamol while muting native stdout/stderr noise."""
    stdout_fd = os.dup(1)
    stderr_fd = os.dup(2)
    devnull_fd = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull_fd, 1)
        os.dup2(devnull_fd, 2)
        return importlib.import_module("datamol")
    finally:
        os.dup2(stdout_fd, 1)
        os.dup2(stderr_fd, 2)
        os.close(devnull_fd)
        os.close(stdout_fd)
        os.close(stderr_fd)
