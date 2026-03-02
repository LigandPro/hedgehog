"""Best-effort nvMolKit runtime enablement helpers."""

from __future__ import annotations

import importlib
import os
import sys
from functools import lru_cache
from pathlib import Path

from hedgehog.configs.logger import logger

_LOGGED_CONTEXTS: set[tuple[str, str, str]] = set()


def _iter_worker_site_packages(project_root: Path):
    venv_dir = project_root / ".venv-nvmolkit-worker"
    if os.name == "nt":
        candidate = venv_dir / "Lib" / "site-packages"
        if candidate.exists():
            yield candidate
        return

    version_tag = f"python{sys.version_info.major}.{sys.version_info.minor}"
    candidate = venv_dir / "lib" / version_tag / "site-packages"
    if candidate.exists():
        yield candidate


def _try_import_nvmolkit() -> tuple[bool, str]:
    try:
        importlib.import_module("nvmolkit")
        return True, "current-environment"
    except ModuleNotFoundError:
        return False, ""
    except Exception as exc:  # noqa: BLE001
        logger.debug("nvMolKit import failed: %s", exc)
        return False, ""


@lru_cache(maxsize=8)
def _resolve_nvmolkit_source(project_root_str: str) -> str:
    imported, source = _try_import_nvmolkit()
    if imported:
        return source

    if not project_root_str:
        return ""

    project_root = Path(project_root_str)
    for site_packages in _iter_worker_site_packages(project_root):
        site_str = str(site_packages)
        if site_str not in sys.path:
            sys.path.insert(0, site_str)

        imported, _ = _try_import_nvmolkit()
        if imported:
            return f"worker-site-packages:{site_str}"

    return ""


def maybe_enable_nvmolkit(
    *,
    project_root: Path | None = None,
    context: str = "runtime",
) -> bool:
    """Enable nvMolKit if it is importable from current env or worker venv."""
    root = ""
    if project_root is not None:
        root = str(project_root.resolve())

    source = _resolve_nvmolkit_source(root)
    if not source:
        return False

    key = (root, context, source)
    if key not in _LOGGED_CONTEXTS:
        logger.info("nvMolKit enabled for %s (%s)", context, source)
        _LOGGED_CONTEXTS.add(key)
    return True
