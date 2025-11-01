import os 

from configs.config_utils import load_config
from logger_config import logger
from src.hedge.stages.descriptors.utils import compute_metrics, filter_molecules, draw_filtered_mols, process_path


def main(data, config, subfolder):
    """
    Computes default set of 22 physicochemical descriptors per molecule using RDKit,
    filters molecules based on configurable thresholds, and generates
    distribution plots.
    
    Args:
        data: DataFrame with molecules (must have 'smiles' column)
        config: Configuration file
        subfolder: Optional subfolder for output
    """
    folder_to_save = process_path(config['folder_to_save'])
    
    if subfolder is not None:
        folder_to_save = folder_to_save + process_path(subfolder)
        folder_to_save = process_path(folder_to_save)

    config_descriptors = load_config(config['config_descriptors'])
    borders = config_descriptors['borders']
    
    os.makedirs(folder_to_save + 'Descriptors', exist_ok=True)
    if not os.path.exists(folder_to_save):
        os.makedirs(folder_to_save)

    if data is None or len(data) == 0:
        logger.warning('No molecules provided for descriptor calculation. Skipping.')
        return None

    metrics_df = compute_metrics(data, folder_to_save, config=config)  

    if config_descriptors['filter_data']:
        filter_molecules(metrics_df, borders, folder_to_save)
        draw_filtered_mols(metrics_df, folder_to_save, config)
    
    return metrics_df