"""Backwards-compatibility shim — use hedgehog.struct_filters instead."""

import importlib as _il
import sys as _sys
import warnings as _w

_w.warn(
    "hedgehog.stages.structFilters is deprecated, use hedgehog.struct_filters instead",
    DeprecationWarning,
    stacklevel=2,
)
_new = _il.import_module("hedgehog.struct_filters")
_sys.modules[__name__] = _new
