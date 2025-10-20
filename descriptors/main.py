import os

from configs.config_utils import load_config
from descriptors.utils import *
from pathlib import Path
from utils.mol_index import assign_mol_idx


def main(data, config, subfolder):
    folder_to_save = process_path(config['folder_to_save'])
    
    if subfolder != None:
        subfolder = process_path(subfolder)
        folder_to_save = folder_to_save + subfolder
        folder_to_save = process_path(folder_to_save)

    config_descriptors = load_config(config['config_descriptors'])
    borders = config_descriptors['borders']
    
    os.makedirs(folder_to_save + 'Descriptors', exist_ok=True)
    if not os.path.exists(folder_to_save):
        os.makedirs(folder_to_save)

    if 'mol_idx' not in data.columns:
        raise ValueError('mol_idx missing in descriptors input; IDs must be assigned once during preprocessing.')

    metrics_df = compute_metrics(data, folder_to_save, config=config)  
    
    if config_descriptors['filter_data']:
        filter_molecules(metrics_df, borders, folder_to_save)
        draw_filtered_mols(metrics_df, folder_to_save, config)
    return