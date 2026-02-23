"""Backwards-compatibility shims for hedgehog.stages.* imports.

The stage packages have been moved to top-level hedgehog.* modules:
  hedgehog.stages.molPrep       -> hedgehog.molprep
  hedgehog.stages.descriptors   -> hedgehog.descriptors
  hedgehog.stages.structFilters -> hedgehog.struct_filters
  hedgehog.stages.synthesis     -> hedgehog.synthesis
  hedgehog.stages.docking       -> hedgehog.docking
  hedgehog.stages.dockingFilters -> hedgehog.docking_filters

These shims emit a DeprecationWarning and will be removed in a future release.
"""

from __future__ import annotations

import importlib
import warnings


def __getattr__(name: str):
    """Lazy re-export old stage names to new flat packages."""
    _STAGE_MAP = {
        "molPrep": "hedgehog.molprep",
        "descriptors": "hedgehog.descriptors",
        "structFilters": "hedgehog.struct_filters",
        "synthesis": "hedgehog.synthesis",
        "docking": "hedgehog.docking",
        "dockingFilters": "hedgehog.docking_filters",
    }
    if name in _STAGE_MAP:
        new_path = _STAGE_MAP[name]
        warnings.warn(
            f"hedgehog.stages.{name} is deprecated, use {new_path} instead",
            DeprecationWarning,
            stacklevel=2,
        )
        return importlib.import_module(new_path)
    raise AttributeError(f"module 'hedgehog.stages' has no attribute {name!r}")
