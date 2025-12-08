import os
import pandas as pd

from hedgehog.configs.logger import logger, load_config
from hedgehog.stages.structFilters.utils import *

def _order_identity_columns(df):
    """Order dataframe columns with identity columns first."""
    id_cols = ['smiles', 'model_name', 'mol_idx']
    existing_id_cols = [c for c in id_cols if c in df.columns]
    ordered_cols = existing_id_cols + [c for c in df.columns if c not in id_cols]
    return df[ordered_cols]


def _get_input_path(config, stage_dir, folder_to_save):
    """Determine input path based on stage directory.

    Args:
        config: Configuration dictionary
        stage_dir: Stage directory path (e.g., 'stages/02_structural_filters_pre' or 'stages/03_structural_filters_post')
        folder_to_save: Base output folder

    Returns:
        Path to input file
    """
    descriptors_enabled = False
    if 'config_descriptors' in config:
        try:
            desc_config = load_config(config['config_descriptors'])
            descriptors_enabled = desc_config.get('run', False)
        except:
            pass
    
    struct_filters_config = load_config(config.get('config_structFilters', ''))
    run_before_descriptors = struct_filters_config.get('run_before_descriptors', False)
    
    is_pre_descriptors = '02_structural_filters_pre' in stage_dir or 'pre' in stage_dir.lower()
    is_post_descriptors = '03_structural_filters_post' in stage_dir or stage_dir == 'StructFilters'

    if is_post_descriptors:
        descriptors_path = os.path.join(folder_to_save, 'stages', '01_descriptors_initial', 'filtered', 'filtered_molecules.csv')
        if os.path.exists(descriptors_path) and os.path.getsize(descriptors_path) > 0:
            return descriptors_path
        legacy_path = folder_to_save + 'Descriptors/passDescriptorsSMILES.csv'
        if os.path.exists(legacy_path) and os.path.getsize(legacy_path) > 0:
            return legacy_path
        
        if descriptors_enabled:
            logger.error("Descriptors are enabled but output not found. Cannot proceed with post-descriptors filters.")
            raise FileNotFoundError(f"Descriptors output not found at {descriptors_path} or {legacy_path}")
        
        pre_filters_path = os.path.join(folder_to_save, 'stages', '02_structural_filters_pre', 'filtered_molecules.csv')
        if os.path.exists(pre_filters_path) and os.path.getsize(pre_filters_path) > 0:
            return pre_filters_path
        
        sampled_path = os.path.join(folder_to_save, 'input', 'sampled_molecules.csv')
        if os.path.exists(sampled_path):
            return sampled_path
        legacy_sampled = os.path.join(folder_to_save, 'sampled_molecules.csv')
        if os.path.exists(legacy_sampled):
            return legacy_sampled
    
    if is_pre_descriptors:
        sampled_path = os.path.join(folder_to_save, 'input', 'sampled_molecules.csv')
        if os.path.exists(sampled_path):
            return sampled_path
        legacy_sampled = os.path.join(folder_to_save, 'sampled_molecules.csv')
        if os.path.exists(legacy_sampled):
            return legacy_sampled
    
    return config['generated_mols_path']


def main(config, stage_dir):
    """Main entry point for structural filters stage.

    Args:
        config: Configuration dictionary
        stage_dir: Stage directory path (e.g., 'stages/02_structural_filters_pre' or 'stages/03_structural_filters_post')
    """
    sample_size = config['sample_size']
    folder_to_save = process_path(config['folder_to_save'])

    output_dir = os.path.join(folder_to_save, stage_dir)
    os.makedirs(output_dir, exist_ok=True)
    input_path = _get_input_path(config, stage_dir, folder_to_save)
    
    try:
        input_df = pd.read_csv(input_path)
        if 'model_name' not in input_df.columns:
            inferred_model = os.path.splitext(os.path.basename(input_path))[0]
            input_df['model_name'] = inferred_model
            logger.info("model_name column missing; using '%s'", inferred_model)
        if 'mol_idx' not in input_df.columns:
            input_df['mol_idx'] = range(len(input_df))
        model_names = sorted(input_df['model_name'].dropna().unique().tolist())
        model_name = model_names[0] if len(model_names) == 1 else model_names
    except Exception as e:
        logger.error(f"Could not load input data from {input_path}: {e}")
        raise
    
    process_input_path = input_path
    
    config_structFilters = load_config(config['config_structFilters'])
    filters_to_calculate = {k.replace('calculate_', ''): v 
                            for k, v in config_structFilters.items() 
                            if 'calculate_' in k and v
                          }
 
    for filter_name in filters_to_calculate:
        apply_func = filter_function_applier(filter_name)
        filter_results = process_one_file(config, process_input_path, apply_func, sample_size)
        if filter_results is None:
            logger.warning(f"No molecules to process")
            continue
        final_res, final_extended = get_basic_stats(config_structFilters, filter_results, model_name, filter_name=filter_name)

        filter_subdir = os.path.join(output_dir, filter_name)
        os.makedirs(filter_subdir, exist_ok=True)

        final_res = _order_identity_columns(final_res)
        final_res.to_csv(os.path.join(filter_subdir, 'metrics.csv'), index=False)

        final_extended = _order_identity_columns(final_extended)
        final_extended.to_csv(os.path.join(filter_subdir, 'extended.csv'), index=False)

        if 'pass' not in final_extended.columns:
            if 'pass_filter' in final_extended.columns:
                final_extended['pass'] = final_extended['pass_filter']
            else:
                logger.warning(f"'pass' column not found in {filter_name} extended results. Assuming all molecules pass.")
                final_extended['pass'] = True
        filtered_mols = final_extended[final_extended['pass'] == True].copy()

        filtered_mols = _order_identity_columns(filtered_mols)
        filtered_mols.to_csv(os.path.join(filter_subdir, 'filtered_molecules.csv'), index=False)

    plot_calculated_stats(config, stage_dir)
    plot_restriction_ratios(config, stage_dir)

    is_post_descriptors = '03_structural_filters_post' in stage_dir or stage_dir == 'StructFilters'
    if is_post_descriptors:
        plot_filter_failures_analysis(config, stage_dir)

    is_single_stage = config.get('_run_single_stage_override') == 'struct_filters'
    if config_structFilters.get('filter_data', False) or is_single_stage:
        filter_data(config, stage_dir)

    inject_identity_columns_to_all_csvs(config, stage_dir)
