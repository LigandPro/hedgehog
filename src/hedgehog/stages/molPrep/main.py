"""Backwards-compatibility shim — use hedgehog.molprep.main instead."""

import importlib as _il
import sys as _sys
import warnings as _w

_w.warn(
    "hedgehog.stages.molPrep.main is deprecated, use hedgehog.molprep.main instead",
    DeprecationWarning,
    stacklevel=2,
)
_new = _il.import_module("hedgehog.molprep.main")
_sys.modules[__name__] = _new
