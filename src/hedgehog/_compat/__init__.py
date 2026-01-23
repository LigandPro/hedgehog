"""Compatibility helpers for third-party integrations."""

from __future__ import annotations

from importlib import import_module
from types import ModuleType


def _get_module(name: str) -> ModuleType | None:
    try:
        return import_module(name)
    except Exception:  # noqa: BLE001
        return None


def ensure_pandas_adjustment() -> None:
    """Expose pandas text adjustment helpers where legacy callers expect them."""
    pandas_format = _get_module("pandas.io.formats.format")
    pandas_printing = _get_module("pandas.io.formats.printing")

    if pandas_format is None or pandas_printing is None:
        return

    if not hasattr(pandas_format, "get_adjustment"):
        pandas_format.get_adjustment = pandas_printing.get_adjustment  # type: ignore[attr-defined]

    if not hasattr(pandas_format, "_get_adjustment") and hasattr(
        pandas_printing, "_get_adjustment"
    ):
        pandas_format._get_adjustment = pandas_printing._get_adjustment  # type: ignore[attr-defined]
