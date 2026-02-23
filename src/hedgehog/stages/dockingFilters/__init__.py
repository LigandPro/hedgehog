"""Backwards-compatibility shim. Use hedgehog.docking_filters instead."""

import importlib
import warnings

warnings.warn(
    "hedgehog.stages.dockingFilters is deprecated. Use hedgehog.docking_filters instead.",
    DeprecationWarning,
    stacklevel=2,
)

from hedgehog.docking_filters import *  # noqa: E402, F401, F403
from hedgehog.docking_filters import docking_filters_main  # noqa: E402, F401


def __getattr__(name):
    return importlib.import_module(f"hedgehog.docking_filters.{name}")
