"""HEDGEHOG reporting module for generating HTML reports with visualizations."""

from __future__ import annotations


def __getattr__(name: str):
    """Lazily load public symbols to avoid eager heavy imports."""
    if name == "ReportGenerator":
        from hedgehog.reporting.report_generator import ReportGenerator

        return ReportGenerator
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    """Expose lazy exports in module introspection."""
    return sorted(set(globals()) | set(__all__))


__all__ = ["ReportGenerator"]
