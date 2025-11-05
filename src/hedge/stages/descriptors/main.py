import os 

from hedge.configs.logger import logger, load_config
from hedge.stages.descriptors.utils import compute_metrics, filter_molecules, draw_filtered_mols, process_path


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
    folder_to_save = process_path(config['folder_to_save'])

    # Determine the descriptors folder path
    if subfolder is not None:
        # subfolder is already the full path like 'stages/01_descriptors_initial'
        descriptors_folder = os.path.join(folder_to_save, subfolder)
    else:
        descriptors_folder = os.path.join(folder_to_save, 'stages', '01_descriptors_initial')

    os.makedirs(descriptors_folder, exist_ok=True)

    # Create subdirectories for organized output
    metrics_folder = os.path.join(descriptors_folder, 'metrics')
    filtered_folder = os.path.join(descriptors_folder, 'filtered')
    plots_folder = os.path.join(descriptors_folder, 'plots')
    os.makedirs(metrics_folder, exist_ok=True)
    os.makedirs(filtered_folder, exist_ok=True)
    os.makedirs(plots_folder, exist_ok=True)

    config_descriptors = load_config(config['config_descriptors'])
    borders = config_descriptors['borders']

    if data is None or len(data) == 0:
        logger.warning('No molecules provided for descriptor calculation. Skipping.')
        return None

    # Compute metrics and save to metrics/ subfolder
    metrics_df = compute_metrics(data, metrics_folder, config=config)

    if config_descriptors['filter_data']:
        # Filter molecules and save to filtered/ subfolder
        filter_molecules(metrics_df, borders, filtered_folder)
        # Draw plots and save to plots/ subfolder
        draw_filtered_mols(metrics_df, plots_folder, config)

    return metrics_df
    