"""Auto-installation helpers for external tools."""

from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from hedgehog.setup._fsscore import FSScoreRuntime
    from hedgehog.setup._gasa import GasaSetupResult
    from hedgehog.setup._nonpher import NonpherCheckResult, NonpherEnsureResult


_LAZY_EXPORTS: dict[str, tuple[str, str]] = {
    "ensure_gnina": ("hedgehog.setup._gnina", "ensure_gnina"),
    "ensure_aizynthfinder": ("hedgehog.setup._aizynthfinder", "ensure_aizynthfinder"),
    "ensure_rascore_model": ("hedgehog.setup._rascore", "ensure_rascore_model"),
    "ensure_sync_model": ("hedgehog.setup._sync", "ensure_sync_model"),
    "ensure_scscore_model": ("hedgehog.setup._scscore", "ensure_scscore_model"),
    "ensure_fsscore_checkout": ("hedgehog.setup._fsscore", "ensure_fsscore_checkout"),
    "ensure_fsscore_runtime": ("hedgehog.setup._fsscore", "ensure_fsscore_runtime"),
    "FSScoreRuntime": ("hedgehog.setup._fsscore", "FSScoreRuntime"),
    "ensure_gasa_worker": ("hedgehog.setup._gasa", "ensure_gasa_worker"),
    "GasaSetupResult": ("hedgehog.setup._gasa", "GasaSetupResult"),
    "ensure_shepherd_worker": (
        "hedgehog.setup._shepherd_worker",
        "ensure_shepherd_worker",
    ),
    "ensure_nvmolkit_worker": (
        "hedgehog.setup._nvmolkit_worker",
        "ensure_nvmolkit_worker",
    ),
    "check_nonpher_runtime": ("hedgehog.setup._nonpher", "check_nonpher_runtime"),
    "create_nonpher_complexity_filter": (
        "hedgehog.setup._nonpher",
        "create_nonpher_complexity_filter",
    ),
    "nonpher_lobachevsky_setup_commands": (
        "hedgehog.setup._nonpher",
        "nonpher_lobachevsky_setup_commands",
    ),
    "NonpherCheckResult": ("hedgehog.setup._nonpher", "NonpherCheckResult"),
    "ensure_nonpher_external_runtime": (
        "hedgehog.setup._nonpher",
        "ensure_nonpher_external_runtime",
    ),
    "ensure_nonpher_uv_runtime": (
        "hedgehog.setup._nonpher",
        "ensure_nonpher_uv_runtime",
    ),
    "resolve_nonpher_python": ("hedgehog.setup._nonpher", "resolve_nonpher_python"),
    "run_nonpher_batch_external": (
        "hedgehog.setup._nonpher",
        "run_nonpher_batch_external",
    ),
    "NONPHER_PYTHON_ENV_VAR": ("hedgehog.setup._nonpher", "NONPHER_PYTHON_ENV_VAR"),
    "NonpherEnsureResult": ("hedgehog.setup._nonpher", "NonpherEnsureResult"),
}


def __getattr__(name: str):
    """Lazily import public API symbols so the package can be imported
    before all sub-modules exist (they may be created by parallel agents)."""
    target = _LAZY_EXPORTS.get(name)
    if target is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    module_name, attr_name = target
    module = import_module(module_name)
    return getattr(module, attr_name)


def __dir__() -> list[str]:
    """Expose lazy exports in module introspection without eager binding."""
    return sorted(set(globals()) | set(__all__))


def ensure_matcha_checkout(*args, **kwargs):
    """Expose Matcha checkout setup as a real module attribute."""
    from hedgehog.setup._matcha import ensure_matcha_checkout as _impl

    return _impl(*args, **kwargs)


__all__ = [
    "ensure_gnina",
    "ensure_aizynthfinder",
    "ensure_rascore_model",
    "ensure_sync_model",
    "ensure_scscore_model",
    "ensure_fsscore_checkout",
    "ensure_fsscore_runtime",
    "FSScoreRuntime",
    "ensure_gasa_worker",
    "GasaSetupResult",
    "ensure_shepherd_worker",
    "ensure_nvmolkit_worker",
    "check_nonpher_runtime",
    "create_nonpher_complexity_filter",
    "nonpher_lobachevsky_setup_commands",
    "NonpherCheckResult",
    "ensure_nonpher_external_runtime",
    "ensure_nonpher_uv_runtime",
    "resolve_nonpher_python",
    "run_nonpher_batch_external",
    "NONPHER_PYTHON_ENV_VAR",
    "NonpherEnsureResult",
    "ensure_matcha_checkout",
]
