"""Auto-installation helpers for external tools."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass


def __getattr__(name: str):
    """Lazily import public API symbols so the package can be imported
    before all sub-modules exist (they may be created by parallel agents)."""
    if name == "ensure_gnina":
        from hedgehog.setup._gnina import ensure_gnina

        return ensure_gnina
    if name == "ensure_aizynthfinder":
        from hedgehog.setup._aizynthfinder import ensure_aizynthfinder

        return ensure_aizynthfinder
    if name == "ensure_rascore_model":
        from hedgehog.setup._rascore import ensure_rascore_model

        return ensure_rascore_model
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = ["ensure_gnina", "ensure_aizynthfinder", "ensure_rascore_model"]
