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


def _is_real_binary(path: str) -> bool:
    """Check if a file is a real compiled binary (ELF), not a script wrapper."""
    try:
        with open(path, "rb") as f:
            header = f.read(4)
        return header == b"\x7fELF"
    except Exception:
        return False


def _resolve_docking_binary(config_path: str, tool_name: str) -> str:
    """Resolve a docking binary path from config or PATH.

    Args:
        config_path: Path from config (absolute path or bare tool name).
        tool_name: Tool name for PATH lookup (e.g. 'smina', 'gnina').

    Returns:
        Resolved absolute path to the binary.

    Raises:
        FileNotFoundError: If the binary cannot be found.
    """
    if os.path.isabs(config_path) and os.path.isfile(config_path):
        return config_path

    found = shutil.which(tool_name)
    if found and _is_real_binary(found):
        return found
    elif found:
        logger.debug(
            "%s found at %s but is a script wrapper, not a real binary — skipping",
            tool_name,
            found,
        )

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
