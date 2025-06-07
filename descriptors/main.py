import os
import pandas as pd

from configs.config_utils import load_config
from descriptors.utils import *
from logger_config import logger

from syba import syba


def main(data, config):
    model_name = config['model_name']
    folder_to_save = config['folder_to_save']
    target_mols_path = config['target_mols_path']

    config_descriptors = load_config(config['config_descriptors'])
    n_jobs = config_descriptors['n_jobs']
    batch_size = config_descriptors['batch_size']
    allowed_chars = config_descriptors['borders']['allowed_chars']
    borders = config_descriptors['borders']
    syba_model = read_syba() 
    
    if not os.path.exists(folder_to_save):
        os.makedirs(folder_to_save)
    
    target_data = load_inhibitors(target_mols_path)
    dict_generated = loading_generated_mols(data, model_name)
    # dict_generated = dict(sorted(dict_generated.items(), key=lambda x: int(''.join(filter(str.isdigit, x[0])))))

    check_intersection(dict_generated, target_data)

    tanimoto_similarity_claculation(dict_generated, target_data, folder_to_save)
    descriptors_dict = compute_descriptors(dict_generated)
    dict_generated = compute_mce18_score(dict_generated, n_jobs)
    dict_generated = compute_syba_score(dict_generated, syba_model, n_jobs)
    target_descriptors_df = computing_target_descriptors(target_data, syba_model, n_jobs)
    dict_generated["Target"] = target_descriptors_df

    metrics_dict = collect_metrics_dict(dict_generated, descriptors_dict, target_descriptors_df)
    save_plots(metrics_dict, dict_generated, folder_to_save, target_data, n_jobs, batch_size) 
    metrics_df = compute_metrics(data, syba_model, model_name, folder_to_save)  

    if config_descriptors['filter_data']:
        filter_molecules(metrics_df, borders, folder_to_save)
        draw_filtered_mols(metrics_df, folder_to_save, config)