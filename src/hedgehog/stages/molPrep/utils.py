"""Backwards-compatibility shim — use hedgehog.molprep.utils instead."""

import importlib as _il
import sys as _sys
import warnings as _w

_w.warn(
    "hedgehog.stages.molPrep.utils is deprecated, use hedgehog.molprep.utils instead",
    DeprecationWarning,
    stacklevel=2,
)
_new = _il.import_module("hedgehog.molprep.utils")
_sys.modules[__name__] = _new
