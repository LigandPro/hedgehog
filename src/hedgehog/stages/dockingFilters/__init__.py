"""Backwards-compatibility shim — use hedgehog.docking_filters instead."""

import warnings as _w

_w.warn(
    "hedgehog.stages.dockingFilters is deprecated, use hedgehog.docking_filters instead",
    DeprecationWarning,
    stacklevel=2,
)
from hedgehog.docking_filters import *  # noqa: F401,F403,E402
