"""Backwards-compatibility shim — use hedgehog.descriptors instead."""

import importlib as _il
import sys as _sys
import warnings as _w

_w.warn(
    "hedgehog.stages.descriptors is deprecated, use hedgehog.descriptors instead",
    DeprecationWarning,
    stacklevel=2,
)
_sys.modules[__name__] = _il.import_module("hedgehog.descriptors")
