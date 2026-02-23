"""Backwards-compatibility shim. Use hedgehog.descriptors instead."""

import importlib
import warnings

warnings.warn(
    "hedgehog.stages.descriptors is deprecated. Use hedgehog.descriptors instead.",
    DeprecationWarning,
    stacklevel=2,
)

from hedgehog.descriptors import *  # noqa: E402, F401, F403
from hedgehog.descriptors import main  # noqa: E402, F401


def __getattr__(name):
    return importlib.import_module(f"hedgehog.descriptors.{name}")
