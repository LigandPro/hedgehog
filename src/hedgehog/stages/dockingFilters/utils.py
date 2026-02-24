"""Backwards-compatibility shim — use hedgehog.docking_filters.utils instead."""

import warnings as _w

_w.warn(
    "hedgehog.stages.dockingFilters.utils is deprecated, use hedgehog.docking_filters.utils instead",
    DeprecationWarning,
    stacklevel=2,
)
from hedgehog.docking_filters.utils import *  # noqa: F401,F403,E402
