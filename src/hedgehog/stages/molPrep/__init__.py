"""Backwards-compatibility shim. Use hedgehog.molprep instead."""

import importlib
import warnings

warnings.warn(
    "hedgehog.stages.molPrep is deprecated. Use hedgehog.molprep instead.",
    DeprecationWarning,
    stacklevel=2,
)

from hedgehog.molprep import *  # noqa: E402, F401, F403
from hedgehog.molprep import main  # noqa: E402, F401


def __getattr__(name):
    return importlib.import_module(f"hedgehog.molprep.{name}")
