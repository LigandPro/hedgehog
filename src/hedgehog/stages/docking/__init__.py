"""Backwards-compatibility shim — use hedgehog.docking instead."""

import importlib as _il
import sys as _sys
import warnings as _w

_w.warn(
    "hedgehog.stages.docking is deprecated, use hedgehog.docking instead",
    DeprecationWarning,
    stacklevel=2,
)
_sys.modules[__name__] = _il.import_module("hedgehog.docking")
