"""Backwards-compatibility shim. Use hedgehog.docking instead."""

import importlib
import warnings

warnings.warn(
    "hedgehog.stages.docking is deprecated. Use hedgehog.docking instead.",
    DeprecationWarning,
    stacklevel=2,
)

from hedgehog.docking import *  # noqa: E402, F401, F403
from hedgehog.docking import main, run_docking  # noqa: E402, F401


def __getattr__(name):
    return importlib.import_module(f"hedgehog.docking.{name}")
