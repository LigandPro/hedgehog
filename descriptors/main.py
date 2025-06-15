import os

from configs.config_utils import load_config
from descriptors.utils import *


def main(data, config):
    folder_to_save = config['folder_to_save']

    config_descriptors = load_config(config['config_descriptors'])
    borders = config_descriptors['borders']
    syba_model = load_syba() 
    
    if not os.path.exists(folder_to_save):
        os.makedirs(folder_to_save)

    metrics_df = compute_metrics(data, syba_model, folder_to_save)  

    if config_descriptors['filter_data']:
        filter_molecules(metrics_df, borders, folder_to_save)
        draw_filtered_mols(metrics_df, folder_to_save, config)
    
    return