"""Backwards-compatibility shim — use hedgehog.docking.main instead."""

import importlib as _il
import sys as _sys
import warnings as _w

_w.warn(
    "hedgehog.stages.docking.main is deprecated, use hedgehog.docking.main instead",
    DeprecationWarning,
    stacklevel=2,
)
_new = _il.import_module("hedgehog.docking.main")
_sys.modules[__name__] = _new
