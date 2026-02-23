"""Backwards-compatibility shim. Use hedgehog.synthesis instead."""

import importlib
import warnings

warnings.warn(
    "hedgehog.stages.synthesis is deprecated. Use hedgehog.synthesis instead.",
    DeprecationWarning,
    stacklevel=2,
)

from hedgehog.synthesis import *  # noqa: E402, F401, F403
from hedgehog.synthesis import main  # noqa: E402, F401


def __getattr__(name):
    return importlib.import_module(f"hedgehog.synthesis.{name}")
