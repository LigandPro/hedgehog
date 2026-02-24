"""Backwards-compatibility shim — use hedgehog.synthesis instead."""

import importlib as _il
import sys as _sys
import warnings as _w

_w.warn(
    "hedgehog.stages.synthesis is deprecated, use hedgehog.synthesis instead",
    DeprecationWarning,
    stacklevel=2,
)
_new = _il.import_module("hedgehog.synthesis")
_sys.modules[__name__] = _new
