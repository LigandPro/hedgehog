"""Auto-installation helpers for external tools."""

from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from hedgehog.setup._fsscore import FSScoreRuntime
    from hedgehog.setup._gasa import GasaSetupResult
    from hedgehog.setup._nonpher import NonpherCheckResult, NonpherEnsureResult


NONPHER_PYTHON_ENV_VAR = "HEDGEHOG_NONPHER_PYTHON"


_LAZY_EXPORTS: dict[str, tuple[str, str]] = {
    "ensure_gnina": ("hedgehog.setup._gnina", "ensure_gnina"),
    "ensure_aizynthfinder": ("hedgehog.setup._aizynthfinder", "ensure_aizynthfinder"),
    "ensure_rascore_model": ("hedgehog.setup._rascore", "ensure_rascore_model"),
    "FSScoreRuntime": ("hedgehog.setup._fsscore", "FSScoreRuntime"),
    "GasaSetupResult": ("hedgehog.setup._gasa", "GasaSetupResult"),
    "ensure_shepherd_worker": (
        "hedgehog.setup._shepherd_worker",
        "ensure_shepherd_worker",
    ),
    "ensure_nvmolkit_worker": (
        "hedgehog.setup._nvmolkit_worker",
        "ensure_nvmolkit_worker",
    ),
    "NonpherCheckResult": ("hedgehog.setup._nonpher", "NonpherCheckResult"),
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


def ensure_sync_model(*args, **kwargs):
    """Expose SYNC model setup as a real module attribute."""
    from hedgehog.setup._sync import ensure_sync_model as _impl

    return _impl(*args, **kwargs)


def ensure_scscore_model(*args, **kwargs):
    """Expose SCScore model setup as a real module attribute."""
    from hedgehog.setup._scscore import ensure_scscore_model as _impl

    return _impl(*args, **kwargs)


def ensure_fsscore_checkout(*args, **kwargs):
    """Expose FSScore checkout setup as a real module attribute."""
    from hedgehog.setup._fsscore import ensure_fsscore_checkout as _impl

    return _impl(*args, **kwargs)


def ensure_fsscore_runtime(*args, **kwargs):
    """Expose FSScore runtime setup as a real module attribute."""
    from hedgehog.setup._fsscore import ensure_fsscore_runtime as _impl

    return _impl(*args, **kwargs)


def ensure_gasa_worker(*args, **kwargs):
    """Expose GASA worker setup as a real module attribute."""
    from hedgehog.setup._gasa import ensure_gasa_worker as _impl

    return _impl(*args, **kwargs)


def check_nonpher_runtime(*args, **kwargs):
    """Expose Nonpher runtime checks as a real module attribute."""
    from hedgehog.setup._nonpher import check_nonpher_runtime as _impl

    return _impl(*args, **kwargs)


def create_nonpher_complexity_filter(*args, **kwargs):
    """Expose Nonpher filter setup as a real module attribute."""
    from hedgehog.setup._nonpher import create_nonpher_complexity_filter as _impl

    return _impl(*args, **kwargs)


def nonpher_lobachevsky_setup_commands(*args, **kwargs):
    """Expose Nonpher Lobachevsky setup commands as a real module attribute."""
    from hedgehog.setup._nonpher import nonpher_lobachevsky_setup_commands as _impl

    return _impl(*args, **kwargs)


def ensure_nonpher_external_runtime(*args, **kwargs):
    """Expose Nonpher external runtime setup as a real module attribute."""
    from hedgehog.setup._nonpher import ensure_nonpher_external_runtime as _impl

    return _impl(*args, **kwargs)


def ensure_nonpher_uv_runtime(*args, **kwargs):
    """Expose Nonpher uv runtime setup as a real module attribute."""
    from hedgehog.setup._nonpher import ensure_nonpher_uv_runtime as _impl

    return _impl(*args, **kwargs)


def resolve_nonpher_python(*args, **kwargs):
    """Expose Nonpher Python resolution as a real module attribute."""
    from hedgehog.setup._nonpher import resolve_nonpher_python as _impl

    return _impl(*args, **kwargs)


def run_nonpher_batch_external(*args, **kwargs):
    """Expose Nonpher external batch execution as a real module attribute."""
    from hedgehog.setup._nonpher import run_nonpher_batch_external as _impl

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
