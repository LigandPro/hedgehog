from pathlib import Path

from hedgehog.configs.logger import load_config, logger
from hedgehog.stages.descriptors.utils import (
    compute_metrics,
    draw_filtered_mols,
    filter_molecules,
    process_path,
)


def main(data, config, subfolder):
    """
    Computes default set of 22 physicochemical descriptors per molecule using RDKit,
    filters molecules based on configurable thresholds, and generates
    distribution plots.

    Args:
        data: DataFrame with molecules (must have 'smiles' column)
        config: Configuration file
        subfolder: Optional subfolder for output (e.g., 'stages/01_descriptors_initial' or 'stages/06_descriptors_final')
    """
    if data is None or len(data) == 0:
        logger.warning("No molecules provided for descriptor calculation. Skipping.")
        return None

    folder_to_save = Path(process_path(config["folder_to_save"]))
    subfolder = subfolder or str(Path("stages") / "01_descriptors_initial")
    descriptors_folder = folder_to_save / subfolder

    metrics_folder = descriptors_folder / "metrics"
    filtered_folder = descriptors_folder / "filtered"
    plots_folder = descriptors_folder / "plots"

    for folder in [descriptors_folder, metrics_folder, filtered_folder, plots_folder]:
        Path(folder).mkdir(parents=True, exist_ok=True)

    config_descriptors = load_config(config["config_descriptors"])
    metrics_df = compute_metrics(data, metrics_folder, config=config)

    if config_descriptors["filter_data"]:
        filter_molecules(metrics_df, config_descriptors["borders"], filtered_folder)
        draw_filtered_mols(metrics_df, plots_folder, config)

    return metrics_df
