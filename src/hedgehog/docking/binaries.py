import os
import shutil
from pathlib import Path

from hedgehog._constants import TOOL_GNINA
from hedgehog.configs.logger import logger


def _validate_optional_tool_path(tool_path, tool_label):
    """Validate optional external tool path and return usable value or None."""
    if not tool_path:
        return None

    path = Path(str(tool_path))
    if path.exists():
        if not path.is_file() or not os.access(path, os.X_OK):
            logger.warning(
                "%s is not executable: %s. Falling back to built-in behavior.",
                tool_label,
                tool_path,
            )
            return None
        return str(path)

    resolved = shutil.which(str(tool_path))
    if resolved:
        return resolved

    logger.warning(
        "%s not found: %s. Falling back to built-in behavior.",
        tool_label,
        tool_path,
    )
    return None


def _is_executable_file(path: str) -> bool:
    """Return True when path points to an executable file."""
    return os.path.isfile(path) and os.access(path, os.X_OK)


def _is_path_like(path_value: str) -> bool:
    """Return True when value should be treated as a filesystem path."""
    return (
        os.path.isabs(path_value)
        or path_value.startswith(".")
        or path_value.startswith("~")
        or "/" in path_value
        or "\\" in path_value
    )


def _resolve_configured_binary_path(
    binary_value: str, tool_name: str, config_dir: Path | None = None
) -> str:
    """Resolve and validate configured binary path."""
    candidate = Path(binary_value).expanduser()
    if not candidate.is_absolute() and config_dir is not None:
        candidate = config_dir / candidate
    candidate = candidate.resolve()

    if not candidate.exists():
        raise FileNotFoundError(
            f"Configured {tool_name} binary does not exist: {candidate}"
        )
    if not candidate.is_file():
        raise FileNotFoundError(
            f"Configured {tool_name} binary is not a file: {candidate}"
        )
    if not os.access(candidate, os.X_OK):
        raise PermissionError(
            f"Configured {tool_name} binary is not executable: {candidate}"
        )
    return str(candidate)


def _resolve_docking_binary(
    config_path: str | None,
    tool_name: str,
    *,
    config_dir: str | Path | None = None,
) -> str:
    """Resolve a docking binary path from config or PATH.

    Args:
        config_path: Optional value from config (path or command name).
        tool_name: Tool name for PATH lookup (e.g. 'smina', 'gnina').
        config_dir: Directory of docking config file (for relative paths).

    Returns:
        Resolved absolute path to the binary.

    Raises:
        FileNotFoundError: If the binary cannot be found.
    """
    config_candidate = str(config_path or "").strip()
    if config_candidate:
        expanded = os.path.expanduser(config_candidate)
        if _is_executable_file(expanded):
            return expanded

        resolved_cfg = shutil.which(config_candidate)
        if resolved_cfg and _is_executable_file(resolved_cfg):
            return resolved_cfg

    found = shutil.which(tool_name)
    if found and _is_executable_file(found):
        return found

    if tool_name == TOOL_GNINA:
        from hedgehog.setup import ensure_gnina

        try:
            return ensure_gnina()
        except RuntimeError as exc:
            raise FileNotFoundError(str(exc)) from exc

    raise FileNotFoundError(
        f"Docking binary '{tool_name}' not found. "
        f"Provide absolute path in config or ensure it's on PATH."
    )
