"""Descriptors stage: physicochemical descriptor computation and filtering."""

from __future__ import annotations


def __getattr__(name: str):
    """Lazily load public symbols to avoid eager heavy imports."""
    if name == "run":
        from hedgehog.descriptors.stage import run

        return run
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    """Expose lazy exports in module introspection."""
    return sorted(set(globals()) | set(__all__))


__all__ = ["run"]
