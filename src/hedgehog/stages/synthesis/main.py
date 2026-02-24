"""Backwards-compatibility shim — use hedgehog.synthesis.main instead."""

import importlib as _il
import sys as _sys
import warnings as _w

_w.warn(
    "hedgehog.stages.synthesis.main is deprecated, use hedgehog.synthesis.main instead",
    DeprecationWarning,
    stacklevel=2,
)
_new = _il.import_module("hedgehog.synthesis.main")
_sys.modules[__name__] = _new
