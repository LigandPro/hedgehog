from pathlib import Path
from typing import Any

import pandas as pd

from hedge.configs.logger import load_config, logger
from hedge.stages.descriptors.utils import (
    compute_metrics,
    draw_filtered_mols,
    filter_molecules,
    process_path,
)


def main(
    data: pd.DataFrame | None,
    config: dict[str, Any],
    subfolder: str | None,
) -> pd.DataFrame | None:
    """Compute and filter molecular descriptors.

    Computes default set of 22 physicochemical descriptors per molecule
    using RDKit, filters molecules based on configurable thresholds, and
    generates distribution plots.

    Args:
        data: DataFrame with molecules (must have 'smiles' column)
        config: Configuration file
        subfolder: Optional subfolder for output

    Returns
    -------
        pd.DataFrame: DataFrame with computed descriptors or None
    """
    folder_to_save = process_path(config["folder_to_save"])

    if subfolder is not None:
        if not subfolder.endswith("/"):
            subfolder = subfolder + "/"
        folder_to_save = folder_to_save + subfolder
        Path(folder_to_save).mkdir(parents=True, exist_ok=True)

    config_descriptors = load_config(config["config_descriptors"])
    borders = config_descriptors["borders"]

    descriptors_folder = folder_to_save + "Descriptors/"
    Path(descriptors_folder).mkdir(parents=True, exist_ok=True)

    if data is None or len(data) == 0:
        logger.warning(
            "No molecules provided for descriptor calculation. Skipping."
        )
        return None

    metrics_df = compute_metrics(data, descriptors_folder)

    if config_descriptors["filter_data"]:
        filter_molecules(metrics_df, borders, descriptors_folder)
        draw_filtered_mols(metrics_df, descriptors_folder, config)

    return metrics_df
