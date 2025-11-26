import os
import glob
import pandas as pd

from hedge.configs.logger import logger, load_config
from hedge.stages.structFilters.utils import *

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
    # For post-descriptors filters, look for descriptors output
    if '03_structural_filters_post' in stage_dir or stage_dir == 'StructFilters':
        # New structure
        descriptors_path = os.path.join(folder_to_save, 'stages', '01_descriptors_initial', 'filtered', 'filtered_molecules.csv')
        if os.path.exists(descriptors_path):
            return descriptors_path
        # Legacy structure
        legacy_path = folder_to_save + 'Descriptors/passDescriptorsSMILES.csv'
        if os.path.exists(legacy_path):
            return legacy_path
        # Fallback to sampled molecules
        sampled_path = os.path.join(folder_to_save, 'input', 'sampled_molecules.csv')
        if os.path.exists(sampled_path):
            logger.info("Descriptors output not found, using sampled_molecules.csv")
            return sampled_path
        # Legacy sampled molecules
        legacy_sampled = folder_to_save + 'sampledMols.csv'
        if os.path.exists(legacy_sampled):
            logger.info("Using legacy sampledMols.csv")
            return legacy_sampled
        logger.info("No processed data found, using molecules from config")
        return config['generated_mols_path']

    # For pre-descriptors filters, use original input
    matched = glob.glob(config['generated_mols_path'])
    if len(matched) > 1:
        sampled_path = os.path.join(folder_to_save, 'input', 'sampled_molecules.csv')
        if os.path.exists(sampled_path):
            return sampled_path
        # Legacy fallback
        return folder_to_save + 'sampledMols.csv'

    single_path = matched[0]
    try:
        df_check = pd.read_csv(single_path)
        lower_cols = {c.lower(): c for c in df_check.columns}
        candidate = lower_cols.get('model_name') or lower_cols.get('name')
        if candidate and df_check[candidate].nunique(dropna=True) > 1:
            sampled_path = os.path.join(folder_to_save, 'input', 'sampled_molecules.csv')
            if os.path.exists(sampled_path):
                return sampled_path
            return folder_to_save + 'sampledMols.csv'
    except Exception:
        pass

    return config['generated_mols_path']


def main(config, stage_dir):
    """Main entry point for structural filters stage.

    Args:
        config: Configuration dictionary
        stage_dir: Stage directory path (e.g., 'stages/02_structural_filters_pre' or 'stages/03_structural_filters_post')
    """
    sample_size = config['sample_size']
    folder_to_save = process_path(config['folder_to_save'])

    # Determine output directory
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
        logger.error(f"Could not load input data: {e}")
        raise
    
    config_structFilters = load_config(config['config_structFilters'])
    filters_to_calculate = {k.replace('calculate_', ''): v 
                            for k, v in config_structFilters.items() 
                            if 'calculate_' in k and v
                          }
 
    for filter_name in filters_to_calculate:
        apply_func = filter_function_applier(filter_name)
        filter_results = process_one_file(config, input_path, apply_func, sample_size)
        
        if filter_results is None:
            logger.warning(f"No molecules to process")
            continue
        
        final_res, final_extended = get_basic_stats(config_structFilters, filter_results, model_name, filter_name=filter_name)

        # Create filter-specific subdirectory
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

    # Only run failure analysis for post-descriptors filters
    is_post_descriptors = '03_structural_filters_post' in stage_dir or stage_dir == 'StructFilters'
    if is_post_descriptors:
        plot_filter_failures_analysis(config, stage_dir)

    is_single_stage = config.get('_run_single_stage_override') == 'struct_filters'
    if config_structFilters.get('filter_data', False) or is_single_stage:
        filter_data(config, stage_dir)

    inject_identity_columns_to_all_csvs(config, stage_dir)
