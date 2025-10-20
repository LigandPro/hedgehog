import os
import json
import pandas as pd
from rdkit import Chem

from structFilters.utils import *

from logger_config import logger 
from configs.config_utils import load_config
from utils.mol_index import assign_mol_idx
from structFilters.utils import inject_identity_columns_to_all_csvs


def main(config, prefix):
    model_name = get_model_name(config)

    sample_size = config['sample_size']
    folder_to_save = process_path(config['folder_to_save'])
    if prefix == 'beforeDescriptors':
        os.makedirs(folder_to_save + f'{prefix}_StructFilters/', exist_ok=True)
    else:
        os.makedirs(folder_to_save + f'StructFilters/', exist_ok=True)

    if prefix == 'beforeDescriptors':
        matched = glob.glob(config['generated_mols_path'])
        if len(matched) > 1:
            path = folder_to_save + 'sampledMols.csv'
        else:
            single_path = matched[0]
            use_sampled = False
            try:
                df_check = pd.read_csv(single_path)
                lower_cols = {c.lower(): c for c in df_check.columns}
                candidate = lower_cols.get('model_name') or lower_cols.get('name')
                if candidate is not None and df_check[candidate].nunique(dropna=True) > 1:
                    use_sampled = True
            except Exception:
                pass
            path = folder_to_save + 'sampledMols.csv' if use_sampled else config['generated_mols_path']
    else:
        path = folder_to_save + 'Descriptors/' + f'passDescriptorsSMILES.csv'
    config_structFilters = load_config(config['config_structFilters'])
    
    filters_to_calculate = {}
    for k, v in config_structFilters.items():
        if 'calculate_' in k:
            k = k.replace('calculate_', '')
            filters_to_calculate[k] = v

    for k, v in filters_to_calculate.items():
        if v:     
            apply_func = filter_function_applier(k)
            filter_results = process_one_file(config, path, apply_func, sample_size)
            if filter_results is not None:
                final_res, final_extended = get_basic_stats(config_structFilters, filter_results, model_name, filter_name=k)
                
                if prefix == 'beforeDescriptors':
                    path_to_save = folder_to_save + f'{prefix}_StructFilters/' + f'{camelcase(k)}'
                else:
                    path_to_save = folder_to_save + f'StructFilters/' + f'{camelcase(k)}'

                if 'smiles' in final_res.columns:
                    id_cols = ['smiles', 'model_name', 'mol_idx']
                    if 'model_name' not in final_res.columns:
                        final_res['model_name'] = 'single'
                    ordered_cols = [c for c in id_cols if c in final_res.columns] + [c for c in final_res.columns if c not in id_cols]
                    final_res = final_res[ordered_cols]
                final_res.to_csv(f'{path_to_save}_metrics.csv', index=False)
                try:
                    base_folder = process_path(config['folder_to_save'])
                    id_path = base_folder + 'Descriptors/passDescriptorsSMILES.csv'
                    src_df = None
                    if os.path.exists(id_path):
                        src_df = pd.read_csv(id_path)
                    else:
                        sample_path = base_folder + 'sampledMols.csv'
                        if os.path.exists(sample_path):
                            src_df = pd.read_csv(sample_path)
                    if src_df is not None and 'smiles' in final_extended.columns:
                        if 'smiles' not in src_df.columns:
                            smiles_col = next((c for c in src_df.columns if c.lower() == 'smiles'), None)
                            if smiles_col:
                                src_df = src_df.rename(columns={smiles_col: 'smiles'})
                        if 'model_name' not in src_df.columns:
                            src_df['model_name'] = 'single'
                        if 'model_name' not in final_extended.columns:
                            final_extended['model_name'] = 'single'
                        if 'mol_idx' in src_df.columns:
                            exact_lookup = {(row['smiles'], row['model_name']): row['mol_idx'] for _, row in src_df.iterrows()}
                            for col in ['mol_idx_x', 'mol_idx_y', 'mol_idx']:
                                if col in final_extended.columns:
                                    final_extended.drop(columns=[col], inplace=True)
                            assigned = [exact_lookup.get((s, m)) for s, m in zip(final_extended['smiles'], final_extended['model_name'])]
                            final_extended['mol_idx'] = assigned
                except Exception:
                    pass
                if 'smiles' in final_extended.columns:
                    id_cols = ['smiles', 'model_name', 'mol_idx']
                    if 'model_name' not in final_extended.columns:
                        final_extended['model_name'] = 'single'
                    ordered_cols = [c for c in id_cols if c in final_extended.columns] + [c for c in final_extended.columns if c not in id_cols]
                    final_extended = final_extended[ordered_cols]
                final_extended.to_csv(f"{path_to_save}_extended.csv", index=False)

                filtered_mols = final_extended[final_extended['full_pass'] == True]
                if 'mol_idx_x' in filtered_mols.columns or 'mol_idx_y' in filtered_mols.columns:
                    filtered_mols['mol_idx'] = filtered_mols.get('mol_idx', filtered_mols.get('mol_idx_x', None))
                    drop_cols = [c for c in ['mol_idx_x', 'mol_idx_y'] if c in filtered_mols.columns]
                    if drop_cols:
                        filtered_mols.drop(columns=drop_cols, inplace=True)
                if 'smiles' in filtered_mols.columns:
                    id_cols = ['smiles', 'model_name', 'mol_idx']
                    if 'model_name' not in filtered_mols.columns:
                        filtered_mols['model_name'] = 'single'
                    ordered_cols = [c for c in id_cols if c in filtered_mols.columns] + [c for c in filtered_mols.columns if c not in id_cols]
                    filtered_mols = filtered_mols[ordered_cols]
                filtered_mols.to_csv(f'{path_to_save}_filteredMols.csv', index=False)
            else:
                logger.warning(f"No molecules to process for {model_name}")
    
    plot_calculated_stats(config, prefix)
    plot_restriction_ratios(config, prefix)

    if config_structFilters['filter_data']:
        filter_data(config, prefix)

    inject_identity_columns_to_all_csvs(config, prefix)

    return
