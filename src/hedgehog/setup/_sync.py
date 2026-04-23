"""Setup helper for the SYNC synthesizability classifier checkpoint."""

from __future__ import annotations

import hashlib
from pathlib import Path

from hedgehog.setup._download import confirm_download, download_with_progress

SYNC_MODEL_URL = (
    "https://raw.githubusercontent.com/XYxiyang/SYNC/main/"
    "targetdiff/pretrained_models/classifier_emb.ckpt"
)
SYNC_MODEL_SHA256 = "dae9fc427ea1df0ba8f1d2cf3699bab526ea40a1126cb3d3bb7b3245b9b277dc"
SYNC_MODEL_RELATIVE_PATH = Path("modules") / "sync" / "classifier_emb.ckpt"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def ensure_sync_model(project_root: Path) -> Path:
    """Ensure the SYNC classifier checkpoint exists locally."""
    model_path = project_root / SYNC_MODEL_RELATIVE_PATH
    if model_path.exists():
        if _sha256(model_path) == SYNC_MODEL_SHA256:
            return model_path
        model_path.unlink()

    if not confirm_download("SYNC classifier checkpoint", "2.7 MB"):
        raise RuntimeError("SYNC classifier checkpoint download declined by user.")

    download_with_progress(
        SYNC_MODEL_URL,
        model_path,
        "Downloading SYNC classifier checkpoint",
    )
    if _sha256(model_path) != SYNC_MODEL_SHA256:
        model_path.unlink(missing_ok=True)
        raise RuntimeError("Downloaded SYNC classifier checkpoint checksum mismatch.")
    return model_path
