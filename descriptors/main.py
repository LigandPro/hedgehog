import os

from configs.config_utils import load_config
from descriptors.utils import *


def main(data, config, mode):
    folder_to_save = config['folder_to_save']

    config_descriptors = load_config(config['config_descriptors'])
    borders = config_descriptors['borders']
    
    os.makedirs(folder_to_save + '/Descriptors', exist_ok=True)
    if not os.path.exists(folder_to_save):
        os.makedirs(folder_to_save)


    metrics_df = compute_metrics(data, folder_to_save, mode, config)  

    if config_descriptors['filter_data']:
        filter_molecules(metrics_df, borders, folder_to_save, mode)
        draw_filtered_mols(metrics_df, folder_to_save, config)
    
    return