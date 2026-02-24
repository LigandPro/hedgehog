"""Backwards-compatibility shim — use hedgehog.docking.utils instead."""

import importlib as _il
import sys as _sys
import warnings as _w

_w.warn(
    "hedgehog.stages.docking.utils is deprecated, use hedgehog.docking.utils instead",
    DeprecationWarning,
    stacklevel=2,
)
_new = _il.import_module("hedgehog.docking.utils")
_sys.modules[__name__] = _new
