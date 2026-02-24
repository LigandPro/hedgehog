"""Backwards-compatibility shim — use hedgehog.molprep instead."""

import importlib as _il
import sys as _sys
import warnings as _w

_w.warn(
    "hedgehog.stages.molPrep is deprecated, use hedgehog.molprep instead",
    DeprecationWarning,
    stacklevel=2,
)
_new = _il.import_module("hedgehog.molprep")
_sys.modules[__name__] = _new
