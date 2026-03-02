"""Wave registry for descriptors stage extensibility."""

from collections.abc import Callable
from typing import Any

_waves: list[Callable] = []


def register(fn: Callable) -> Callable:
    """Register a wave function for the descriptors stage.

    Wave functions receive a context dict and return a (possibly modified) context dict.
    They are executed in registration order after the main stage logic.
    """
    _waves.append(fn)
    return fn


def run_waves(context: Any) -> Any:
    """Execute all registered wave functions in order.

    Args:
        context: Stage context dict passed through each wave.

    Returns:
        The context after all waves have been applied.
    """
    for wave_fn in _waves:
        context = wave_fn(context)
    return context
