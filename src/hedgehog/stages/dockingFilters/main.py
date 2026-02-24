"""Backwards-compatibility shim — use hedgehog.docking_filters.main instead."""

import warnings as _w

_w.warn(
    "hedgehog.stages.dockingFilters.main is deprecated, use hedgehog.docking_filters.main instead",
    DeprecationWarning,
    stacklevel=2,
)
from hedgehog.docking_filters.main import *  # noqa: F401,F403,E402
