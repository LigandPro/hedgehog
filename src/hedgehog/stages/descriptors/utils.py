"""Backwards-compatibility shim — use hedgehog.descriptors.utils instead."""

import importlib as _il
import sys as _sys
import warnings as _w

_w.warn(
    "hedgehog.stages.descriptors.utils is deprecated, use hedgehog.descriptors.utils instead",
    DeprecationWarning,
    stacklevel=2,
)
_sys.modules[__name__] = _il.import_module("hedgehog.descriptors.utils")
