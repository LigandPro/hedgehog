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
    """Expose pandas helpers where legacy third-party callers still expect them."""
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

    # pandas 3 removed DataFrame.applymap in favor of DataFrame.map.
    # Older RDKit PandasPatcher still accesses applymap at import time.
    pandas_module = _get_module("pandas")
    if pandas_module is None or not hasattr(pandas_module, "DataFrame"):
        return

    dataframe_class = pandas_module.DataFrame
    if not hasattr(dataframe_class, "applymap") and hasattr(dataframe_class, "map"):
        dataframe_class.applymap = dataframe_class.map  # type: ignore[attr-defined]
