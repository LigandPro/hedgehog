"""Shared download utilities for auto-installing external tools."""

from __future__ import annotations

import os
import urllib.request
from pathlib import Path

from rich.progress import (
    BarColumn,
    DownloadColumn,
    Progress,
    TimeRemainingColumn,
    TransferSpeedColumn,
)


def confirm_download(tool_name: str, size_hint: str) -> bool:
    """Prompt user to confirm downloading a missing tool.

    Respects environment variable overrides:
      - HEDGEHOG_AUTO_INSTALL=1  → auto-accept (no prompt)
      - HEDGEHOG_NON_INTERACTIVE=1 → auto-decline (no prompt)

    Returns True if the user accepts, False otherwise.
    """
    if os.environ.get("HEDGEHOG_AUTO_INSTALL") == "1":
        return True
    if os.environ.get("HEDGEHOG_NON_INTERACTIVE") == "1":
        return False

    try:
        answer = input(f"{tool_name} not found. Download (~{size_hint})? [Y/n] ")
    except EOFError:
        return False

    return answer.strip().lower() in ("", "y", "yes")


def download_with_progress(url: str, dest: Path, description: str) -> None:
    """Download a file with a Rich progress bar and atomic write.

    The file is first written to a temporary path (``dest.downloading``)
    and renamed on success.  On any error the temporary file is removed.

    Args:
        url: URL to download.
        dest: Final destination path (parent dirs created automatically).
        description: Label shown in the progress bar.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".downloading")

    try:
        response = urllib.request.urlopen(url)  # noqa: S310
        total = int(response.headers.get("Content-Length", 0))

        with (
            Progress(
                BarColumn(),
                DownloadColumn(),
                TransferSpeedColumn(),
                TimeRemainingColumn(),
            ) as progress,
            open(tmp, "wb") as fh,
        ):
            task = progress.add_task(description, total=total or None)
            while True:
                chunk = response.read(1024 * 64)
                if not chunk:
                    break
                fh.write(chunk)
                progress.advance(task, len(chunk))

        os.rename(tmp, dest)
    except BaseException:
        if tmp.exists():
            tmp.unlink()
        raise
