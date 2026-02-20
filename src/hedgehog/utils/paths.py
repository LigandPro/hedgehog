"""Shared path utilities."""

from pathlib import Path


def process_path(folder_to_save, key_word=None):
    """Ensure path ends with '/' and create directory if needed."""
    folder_to_save = str(folder_to_save)
    if not folder_to_save.endswith("/"):
        folder_to_save += "/"

    if key_word:
        folder_to_save += f"{key_word}/"

    Path(folder_to_save).mkdir(parents=True, exist_ok=True)
    return folder_to_save


def resolve_existing_path(
    base: Path, path: str | Path, project_root: Path | None = None
) -> Path:
    """Resolve a possibly-relative path to an existing absolute path."""
    p = Path(path)
    if p.is_absolute():
        return p

    candidates = [(base / p).resolve()]
    if project_root is not None:
        candidates.append((project_root / p).resolve())
    candidates.append((Path.cwd() / p).resolve())

    for c in candidates:
        if c.exists():
            return c

    return (base / p).resolve()
