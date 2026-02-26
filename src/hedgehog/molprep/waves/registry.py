from __future__ import annotations

from collections.abc import Callable
from typing import Any

_waves: list[Callable] = []


def register(fn: Callable) -> Callable:
    _waves.append(fn)
    return fn


def run_waves(context: Any) -> Any:
    for wave_fn in _waves:
        context = wave_fn(context)
    return context
