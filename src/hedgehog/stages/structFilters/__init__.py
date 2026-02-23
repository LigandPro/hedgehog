"""Backwards-compatibility shim. Use hedgehog.struct_filters instead."""

import importlib
import warnings

warnings.warn(
    "hedgehog.stages.structFilters is deprecated. Use hedgehog.struct_filters instead.",
    DeprecationWarning,
    stacklevel=2,
)

from hedgehog.struct_filters import *  # noqa: E402, F401, F403
from hedgehog.struct_filters import main  # noqa: E402, F401


def __getattr__(name):
    return importlib.import_module(f"hedgehog.struct_filters.{name}")
