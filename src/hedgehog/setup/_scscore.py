"""Auto-download of the SCScore model used by the synthesis stage."""

from __future__ import annotations

from pathlib import Path

from hedgehog.configs.logger import logger
from hedgehog.setup._download import download_with_progress

_SCSCORE_MODEL_URL = (
    "https://raw.githubusercontent.com/connorcoley/scscore/master/"
    "models/full_reaxys_model_1024bool/model.ckpt-10654.as_numpy.json.gz"
)
_MIN_MODEL_SIZE_BYTES = 1_000_000


def ensure_scscore_model(project_root: Path) -> Path:
    """Ensure SCScore standalone numpy model weights exist locally."""
    model_path = (
        project_root
        / "modules"
        / "scscore"
        / "models"
        / "full_reaxys_model_1024bool"
        / "model.ckpt-10654.as_numpy.json.gz"
    )

    if model_path.exists() and model_path.stat().st_size >= _MIN_MODEL_SIZE_BYTES:
        return model_path

    logger.info("Downloading SCScore model to %s", model_path)
    download_with_progress(_SCSCORE_MODEL_URL, model_path, "SCScore model")

    if not model_path.exists() or model_path.stat().st_size < _MIN_MODEL_SIZE_BYTES:
        model_path.unlink(missing_ok=True)
        raise RuntimeError(
            "Downloaded SCScore model is missing or invalid (too small)."
        )

    logger.info("SCScore model downloaded successfully: %s", model_path)
    return model_path
