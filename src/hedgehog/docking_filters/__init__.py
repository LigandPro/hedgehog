"""Docking filters module for post-docking pose quality assessment."""

from __future__ import annotations


def __getattr__(name: str):
    """Lazily load public symbols to avoid eager heavy imports."""
    if name == "docking_filters_main":
        from .main import docking_filters_main

        return docking_filters_main
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    """Expose lazy exports in module introspection."""
    return sorted(set(globals()) | set(__all__))


__all__ = ["docking_filters_main"]
