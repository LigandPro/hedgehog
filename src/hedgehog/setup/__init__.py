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
    if name == "ensure_shepherd_worker":
        from hedgehog.setup._shepherd_worker import ensure_shepherd_worker

        return ensure_shepherd_worker
    if name == "ensure_nvmolkit_worker":
        from hedgehog.setup._nvmolkit_worker import ensure_nvmolkit_worker

        return ensure_nvmolkit_worker
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    """Expose lazy exports in module introspection without eager binding."""
    return sorted(set(globals()) | set(__all__))


def ensure_sync_model(*args, **kwargs):
    """Expose SYNC model setup as a real module attribute."""
    from hedgehog.setup._sync import ensure_sync_model as _impl

    return _impl(*args, **kwargs)


def ensure_matcha_checkout(*args, **kwargs):
    """Expose Matcha checkout setup as a real module attribute."""
    from hedgehog.setup._matcha import ensure_matcha_checkout as _impl

    return _impl(*args, **kwargs)


__all__ = [
    "ensure_gnina",
    "ensure_aizynthfinder",
    "ensure_rascore_model",
    "ensure_sync_model",
    "ensure_shepherd_worker",
    "ensure_nvmolkit_worker",
    "ensure_matcha_checkout",
]
