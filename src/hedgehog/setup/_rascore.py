"""Auto-download of RAScore model used by synthesis stage."""

from __future__ import annotations

from pathlib import Path

from hedgehog.configs.logger import logger
from hedgehog.setup._download import download_with_progress

_RASCORE_MODEL_URL = (
    "https://raw.githubusercontent.com/MorganCThomas/MolScore/main/"
    "molscore/data/models/RAScore/XGB_chembl_ecfp_counts/model.pkl"
)
_MIN_MODEL_SIZE_BYTES = 1_000_000


def ensure_rascore_model(project_root: Path) -> Path:
    """Ensure RAScore model exists locally and return its path.

    The model is downloaded from the official MolScore repository if missing.
    """
    model_path = (
        project_root
        / "modules"
        / "MolScore"
        / "molscore"
        / "data"
        / "models"
        / "RAScore"
        / "XGB_chembl_ecfp_counts"
        / "model.pkl"
    )

    if model_path.exists() and model_path.stat().st_size >= _MIN_MODEL_SIZE_BYTES:
        return model_path

    logger.info("Downloading RAScore model to %s", model_path)
    download_with_progress(_RASCORE_MODEL_URL, model_path, "RAScore model")

    if not model_path.exists() or model_path.stat().st_size < _MIN_MODEL_SIZE_BYTES:
        model_path.unlink(missing_ok=True)
        raise RuntimeError(
            "Downloaded RAScore model is missing or invalid (too small)."
        )

    logger.info("RAScore model downloaded successfully: %s", model_path)
    return model_path
