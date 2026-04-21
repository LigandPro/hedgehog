"""Docking filters module for post-docking pose quality assessment."""

from __future__ import annotations


def __getattr__(name: str):
    """Lazily load public symbols to avoid eager heavy imports."""
    if name == "docking_filters_main":
        from .main import docking_filters_main

        return docking_filters_main
    if name == "apply_pose_quality_filter":
        from .utils import apply_pose_quality_filter

        return apply_pose_quality_filter
    if name == "apply_interaction_filter":
        from .utils import apply_interaction_filter

        return apply_interaction_filter
    if name == "apply_shepherd_score_filter":
        from .utils import apply_shepherd_score_filter

        return apply_shepherd_score_filter
    if name == "apply_conformer_deviation_filter":
        from .utils import apply_conformer_deviation_filter

        return apply_conformer_deviation_filter
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    """Expose lazy exports in module introspection."""
    return sorted(set(globals()) | set(__all__))


__all__ = [
    "docking_filters_main",
    "apply_pose_quality_filter",
    "apply_interaction_filter",
    "apply_shepherd_score_filter",
    "apply_conformer_deviation_filter",
]
