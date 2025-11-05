"""Environment-wide compatibility shims loaded before any other modules."""

from __future__ import annotations

from hedge._compat import ensure_pandas_adjustment

ensure_pandas_adjustment()

