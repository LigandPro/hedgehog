import os

from structFilters.utils import *

from logger_config import logger 
from configs.config_utils import load_config


def main(config, mode, prefix):
    model_name = get_model_name(config, mode=mode)

    sample_size = config['sample_size']
    folder_to_save = config['folder_to_save']

    os.makedirs(folder_to_save + f'/{prefix}StructFilters', exist_ok=True)
    if prefix == 'before_descriptors':
        if mode == 'single_comparison':
            path = config['generated_mols_path']
    else:
        path = folder_to_save + '/Descriptors/' + f'leftMolsAfterDescriptorsSMILES.csv'

    config_structFilters = load_config(config['config_structFilters'])
    
    filters_to_calculate = {}
    for k, v in config_structFilters.items():
        if 'calculate_' in k:
            k = k.replace('calculate_', '')
            filters_to_calculate[k] = v

    for k, v in filters_to_calculate.items():
        if v:     
            apply_func = filter_function_applier(k)
            filter_results = process_one_file(config, path, apply_func, sample_size, mode)
            if filter_results is not None:
                final_res, final_extended = get_basic_stats(config_structFilters, filter_results, model_name, filter_name=k, mode=mode)
                
                path_to_save = folder_to_save + f'/{prefix}StructFilters/' + f'{camelcase(k)}'

                final_res.to_csv(f'{path_to_save}metrics.csv', index=False)
                final_extended.to_csv(f"{path_to_save}extended.csv", index=False)

                filtered_mols = dropFalse(final_extended)
                filtered_mols.to_csv(f'{path_to_save}filteredMols.csv', index=False)
            else:
                logger.warning(f"No molecules to process for {model_name}")
    
    plot_calculated_stats(config, prefix)
    plot_restriction_ratios(config, prefix)

    if config_structFilters['filter_data']:
        filter_data(config, prefix)

    return
