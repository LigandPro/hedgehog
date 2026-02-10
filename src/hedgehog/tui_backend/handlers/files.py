"""File handling for TUI backend."""

from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from ..server import JsonRpcServer


def _validate_path(resolved_path: Path, allowed_base: Path | None = None) -> None:
    """Validate that a resolved path is within the allowed base directory.

    Raises ValueError if the path escapes the allowed directory tree.
    """
    if allowed_base is None:
        allowed_base = Path.home()
    try:
        resolved_path.relative_to(allowed_base)
    except ValueError as err:
        raise ValueError(
            f"Access denied: path '{resolved_path}' is outside '{allowed_base}'"
        ) from err


def _normalize_extensions(extensions: list[str] | None) -> set[str] | None:
    """Normalize file extensions to lowercase without leading dots."""
    if extensions is None:
        return None
    return {ext.lstrip(".").lower() for ext in extensions}


def _resolve_and_validate_dir(path: str) -> Path:
    """Resolve a path, validate it, and ensure it is an existing directory."""
    dir_path = Path(path).expanduser().resolve()
    _validate_path(dir_path)
    if not dir_path.exists():
        raise FileNotFoundError(f"Directory not found: {path}")
    if not dir_path.is_dir():
        raise NotADirectoryError(f"Not a directory: {path}")
    return dir_path


class FilesHandler:
    """Handler for file-related RPC methods."""

    def __init__(self, server: "JsonRpcServer"):
        self.server = server

    def list_files(self, path: str, extensions: list[str] | None = None) -> list[str]:
        """List files in a directory, optionally filtered by extension."""
        dir_path = _resolve_and_validate_dir(path)
        ext_set = _normalize_extensions(extensions)

        files = []
        for item in sorted(dir_path.iterdir()):
            if item.is_file():
                if ext_set is None or item.suffix.lstrip(".").lower() in ext_set:
                    files.append(str(item))

        return files

    def list_directory(
        self, path: str, extensions: list[str] | None = None
    ) -> list[dict[str, Any]]:
        """List directory contents with metadata."""
        dir_path = _resolve_and_validate_dir(path)
        ext_set = _normalize_extensions(extensions)

        entries = []

        try:
            for item in sorted(
                dir_path.iterdir(), key=lambda x: (not x.is_dir(), x.name.lower())
            ):
                # Skip hidden files
                if item.name.startswith("."):
                    continue

                is_dir = item.is_dir()

                # Filter by extension for files
                if not is_dir and ext_set is not None:
                    file_ext = item.suffix.lstrip(".").lower()
                    if file_ext not in ext_set:
                        continue

                entry = {
                    "name": item.name,
                    "path": str(item),
                    "isDirectory": is_dir,
                }

                if not is_dir:
                    try:
                        stat = item.stat()
                        entry["size"] = stat.st_size
                        entry["modified"] = stat.st_mtime
                    except OSError:
                        pass

                entries.append(entry)
        except PermissionError as e:
            raise PermissionError(f"Permission denied: {path}") from e

        return entries

    def count_molecules(self, path: str) -> dict[str, Any]:
        """Count molecules in a file (CSV, SMI, SDF)."""
        file_path = Path(path).expanduser().resolve()
        _validate_path(file_path)

        if not file_path.exists():
            raise FileNotFoundError(f"File not found: {path}")

        if not file_path.is_file():
            raise ValueError(f"Not a file: {path}")

        ext = file_path.suffix.lower()
        count = 0

        try:
            if ext in (".csv",):
                # CSV: count lines minus header
                with open(file_path) as f:
                    count = sum(1 for _ in f) - 1  # subtract header
                    count = max(0, count)
            elif ext in (".smi", ".smiles"):
                # SMI: count non-empty lines
                with open(file_path) as f:
                    count = sum(1 for line in f if line.strip())
            elif ext in (".sdf", ".mol2"):
                # SDF: count $$$$ delimiters using streaming to avoid OOM on large files
                with open(file_path) as f:
                    count = sum(1 for line in f if "$$$$" in line)
                    if count == 0:
                        # Fallback: might be single molecule without delimiter
                        count = 1 if file_path.stat().st_size > 0 else 0
            else:
                # Unknown format - try counting lines
                with open(file_path) as f:
                    count = sum(1 for line in f if line.strip())
        except Exception as e:
            return {"count": 0, "error": str(e)}

        return {"count": count, "path": str(file_path)}
