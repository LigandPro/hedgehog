"""Shared path helpers for stage outputs."""

from __future__ import annotations

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
