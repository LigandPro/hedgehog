
import os 
import glob 
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


from rdkit import Chem
from tqdm.auto import tqdm
from functools import reduce

import datamol as dm
import medchem as mc

from logger_config import logger
from configs.config_utils import load_config


def camelcase(any_str):
    camel_str = ''.join(word.capitalize() for word in any_str.split('_'))
    return camel_str


def get_model_name(config, mode):
    assert mode in ['single_comparison', 'multi_comparison'], "Mode must be either 'single_comparison' or 'multi_comparison'"
    
    if mode == 'single_comparison':
        return config['generated_mols_path'].split('/')[-1].split('.')[0]
    else:
        paths = glob.glob(config['generated_mols_path'])
        model_names = [path.split('/')[-1].split('.')[0] for path in paths]
        return model_names


def dropFalse(df):
    df = df[df['full_pass'] == True]
    return df


def sdf_to_mols(sdf_file, SUBSAMPLE):
    mols_list = []
    smiles_list = []
    
    molecules = dm.read_sdf(sdf_file)
    i = 0
    for mol in molecules:
        if i >= SUBSAMPLE:
            break
        i += 1
        if mol is not None: 
            mols_list.append(mol)
            smiles = dm.to_smiles(mol)
            smiles_list.append(smiles)
    
    return mols_list, smiles_list


def dropna(mols, smiles):
    df = pd.DataFrame({
        'mols' : mols,
        'smiles' : smiles
    })

    df.dropna(inplace=True)
    return df.mols.tolist(), df.smiles.tolist()


def format_number(x, p=None):
    if x >= 1e6:
        return f'{x/1e6:.1f}M'
    elif x >= 1e3:
        return f'{x/1e3:.1f}K'
    else:
        return f'{x:.0f}'


def get_model_colors(model_names, cmap=None):
   return dict(zip(model_names, plt.cm.YlOrRd(np.linspace(1, 0, len(model_names) + 1)) 
                   if cmap is None 
                   else plt.colormaps.get_cmap(cmap)(np.linspace(1, 0, len(model_names) + 1))))


def clean_name(name):
    to_replace = ['metrics', '_', '.csv']
    for ban in to_replace:
        name = name.replace(ban, '')
    return name.strip()
    


def process_path(folder_to_save, key_word=None):
    if not folder_to_save.endswith('/'):
        folder_to_save = folder_to_save + '/'

    if key_word:
        folder_to_save = folder_to_save + f'{key_word}/'

    return folder_to_save


def common_postprocessing_statistics(filter_results, res_df, stat, extend):          
    if stat is not None:
        res_df = pd.concat([stat, res_df])

    to_smi = lambda mol: Chem.MolToSmiles(mol, canonical=True)
    filter_results["SMILES"] = filter_results.mol.apply(to_smi)
    filter_results = filter_results.drop(columns='mol')

    if extend is not None:
        filter_extended = pd.concat([extend, filter_results]).copy()
    else:
        filter_extended = filter_results.copy()

    return res_df,filter_extended


def process_one_file(config, input_path, apply_filter, subsample):
    mode = config['mode']
    
    input_type = input_path[input_path.rfind(".")+1:]
    assert input_type in {"csv", "smi", 'sdf', 'txt'}

    if input_type == 'csv':
        data = pd.read_csv(input_path)
        smiles_col = None
        if 'smiles' in data.columns: smiles_col = 'smiles'
        elif 'SMILES' in data.columns: smiles_col = 'SMILES'

        if ('SMILES' not in data.columns) and ('smiles' not in data.columns):
            for col in data.columns:
                if "smiles" in col.lower():
                    smiles_col = col
                    break

        assert smiles_col is not None

        if mode == 'single_comparison':
            smiles = data[smiles_col].tolist()
            mols = [dm.to_mol(x) for x in smiles]
        else:
            assert 'model_name' in data.columns, "CSV must contain 'model_name' column in multi_comparison mode"
            smiles_str = data[smiles_col].tolist()
            model_names = data['model_name'].tolist()
            mols = [dm.to_mol(x) for x in smiles_str]
            smiles = list(zip(smiles_str, model_names, mols))

    elif input_type == "smi" or input_type == 'txt':
        with open(input_path, 'r') as file:
            lines = [line.rstrip('\n') for line in file]
        if subsample <= len(lines):
            lines = np.random.permutation(lines)[:subsample].tolist()

        if mode == 'single_comparison':
            smiles = lines
            mols = [dm.to_mol(x) for x in smiles]
        else:
            smiles = []
            model_names = []
            for line in lines:
                parts = line.split()
                if len(parts) == 2:
                    smi, model = parts
                else:
                    smi, model = parts[:1], "unknown"
                smiles.append(smi)
                model_names.append(model)
            mols = [dm.to_mol(x) for x in smiles]
            smiles = list(zip(smiles, model_names, mols))

    elif input_type == 'sdf':
        assert mode == 'single_comparison', "SDF must be in single_comparison mode"
        mols, smiles = sdf_to_mols(input_path, subsample)

    if mode == 'single_comparison':
        mols, smiles = dropna(mols, smiles)
    else:
        smiles = [(smi, model, mol) for (smi, model, mol) in smiles if mol is not None]
        mols = [mol for (_, _, mol) in smiles]

    assert len(mols) == len(smiles), f"{len(mols)}, {len(smiles)}"
    if mode == 'single_comparison':
        assert len(mols) <= subsample
    
    for mol, smi in zip(mols, smiles):
        if mode == 'multi_comparison':
            smi_val = smi[0]
        else:
            smi_val = smi
        assert mol is not None, f"{smi_val}"

    final_result = None
    if len(mols) > 0:
        if mode == 'multi_comparison':
            final_result = apply_filter(config, mols, smiles)
        else:
            final_result = apply_filter(config, mols)

    return final_result


def add_model_name_col(final_result, smiles_with_model, mode):
    if len(final_result) == len(smiles_with_model):
        if mode == 'single_comparison':
            smiles = smiles_with_model
            final_result['smiles'] = smiles
            final_result['model_name'] = model_names 
        else:
            smiles, model_names, mols = zip(*smiles_with_model)
            final_result['smiles'] = smiles
            final_result['model_name'] = model_names
    else:
        logger.warning("Length mismatch between results and smiles_with_model! Attempting to align valid mols.")
        valid_pairs = []
        mols_in_result = list(final_result['mol'])
        for mol in mols_in_result:
            for smi, model, orig_mol in smiles_with_model:
                if mol == orig_mol:
                    valid_pairs.append((smi, model))
                    break
        if len(valid_pairs) == len(final_result):
            smiles, model_names = zip(*valid_pairs)
            final_result['smiles'] = smiles
            final_result['model_name'] = model_names
        else:
            logger.error("Could not align all mols with smiles/model_name. Check your data integrity.")

    return final_result


# list of functions to apply filters
def filter_function_applier(filter_name):
    if filter_name == 'common_alerts':
        return apply_structural_alerts
    elif filter_name == 'NIBR':
        return apply_nibr_filter
    elif filter_name == 'bredt':
        return apply_bredt_filter
    elif filter_name == 'molgraph_stats':
        return apply_molgraph_stats
    elif filter_name == 'molcomplexity':
        return apply_molcomplexity_filters
    else:
        raise ValueError(f"Filter {filter_name} not found")


def apply_structural_alerts(config, mols, smiles_modelName_mols=None):
    mode = config['mode']
    n_jobs = config['n_jobs']

    config_structFilters = load_config(config['config_structFilters'])
    alert_names = config_structFilters['common_alerts_names']
    itoalertname = {i: n for i, n in enumerate(alert_names.keys())}
    
    logger.info(f"Processing {len(mols)} filtered molecules")
    result_frames = []

    for name in tqdm(alert_names):
        alert = mc.structural.CommonAlertsFilters(alerts_set=[name])
        results = alert(mols=mols,
                        n_jobs=n_jobs,
                        )
        result_frames.append(results)

    final_result = result_frames[0]
    final_result['full_pass'] = final_result["pass_filter"].copy()
    final_result['any_pass'] = final_result["pass_filter"].copy()
    final_result[f'pass_filter_{itoalertname[0]}'] = final_result["pass_filter"].copy()
    final_result[f'reasons_{itoalertname[0]}'] = final_result["reasons"].copy()
    
    for i in tqdm(range(1, len(result_frames))):  
        res = result_frames[i]
        name = itoalertname[i]

        final_result[f"pass_filter_{name}"] = res['pass_filter'].copy()
        final_result[f"reasons_{name}"] = res['reasons'].copy()
        
        final_result['full_pass'] = final_result['full_pass'] * final_result[f"pass_filter_{name}"]
        final_result['any_pass'] = final_result['any_pass'] + final_result[f"pass_filter_{name}"]

    if smiles_modelName_mols is not None:
        final_result = add_model_name_col(final_result, smiles_modelName_mols, mode)
    return final_result


def apply_bredt_filter(config, mols: list[Chem.Mol], smiles_modelName_mols=None):
    logger.info(f"Processing {len(mols)} molecules to calculate Bredt filter")
    mode = config['mode']
    out = mc.functional.bredt_filter(mols=mols,
                                     n_jobs=-1,
                                     progress=False,
                                     return_idx=False,
                                    )
    results = pd.DataFrame({'mol' : mols,
                            'full_pass' : out 
                            })    
    if smiles_modelName_mols is not None:
        results = add_model_name_col(results, smiles_modelName_mols, mode) 
    return results


def apply_nibr_filter(config, mols: list[Chem.Mol], smiles_modelName_mols=None):
    logger.info(f"Processing {len(mols)} molecules to calculate NIBR filter")
    mode = config['mode']
    n_jobs = config['n_jobs']

    config_structFilters = load_config(config['config_structFilters'])
    scheduler = config_structFilters['nibr_scheduler']

    nibr_filters = mc.structural.NIBRFilters()
    results = nibr_filters(mols=mols,
                           n_jobs=n_jobs,
                           scheduler=scheduler,
                           keep_details=True,
                        ) 

    if smiles_modelName_mols is not None:
        results = add_model_name_col(results, smiles_modelName_mols, mode)      
    return results


def apply_molgraph_stats(config, mols: list[Chem.rdchem.Mol], smiles_modelName_mols=None):
    logger.info(f"Processing {len(mols)} molecules to calculate Molecular Graph statistics")
    mode = config['mode']
    severities = list(range(1, 12))

    results = {'mol' : mols}
    for s in severities:
        out = mc.functional.molecular_graph_filter(mols=mols,
                                                   max_severity=s,
                                                   n_jobs=-1,
                                                   progress=False,
                                                   return_idx=False,
                                                  )
        results[f'pass_filter_{s}'] = out
    results = pd.DataFrame(results)

    if smiles_modelName_mols is not None:
        results = add_model_name_col(results, smiles_modelName_mols, mode)
    return results


def apply_molcomplexity_filters(config, mols: list[Chem.Mol], smiles_modelName_mols=None):
    logger.info(f"Processing {len(mols)} molecules to calculate Complexity filters")
    mode = config['mode']
    final_result = pd.DataFrame({'mol' : mols,
                                 'full_pass' : True,
                                 'pass_any' : False
                               })
    alert_names = mc.complexity.ComplexityFilter.list_default_available_filters()
    for name in tqdm(alert_names):
        cfilter = mc.complexity.ComplexityFilter(complexity_metric=name)
        final_result[f"pass_{name}"] = final_result["mol"].apply(cfilter)
        final_result['full_pass'] = final_result['full_pass'] * final_result[f"pass_{name}"]
        final_result['pass_any'] = final_result['pass_any'] + final_result[f"pass_{name}"]
    
    if smiles_modelName_mols is not None:
        final_result = add_model_name_col(final_result, smiles_modelName_mols, mode)
    return final_result


def get_basic_stats(config, filter_results: pd.DataFrame, model_name:str, filter_name:str, mode, stat=None, extend=None):
    if mode == 'multi_comparison':
        all_res = []
        all_extended = []
        for model, group in filter_results.groupby('model_name'):
            res_df, filter_extended = get_basic_stats(config, group.copy(), model, filter_name, 'single_comparison', stat, extend,)
            all_res.append(res_df)
            all_extended.append(filter_extended)

        res_df = pd.concat(all_res, ignore_index=True)
        filter_extended = pd.concat(all_extended, ignore_index=True)
        return res_df, filter_extended

    
    num_mol = len(filter_results)

    filter_results.dropna(subset='mol', inplace=True)
    filter_results['model_name'] = model_name 

    if filter_name == 'common_alerts':
        alert_names = config['common_alerts_names']
        
        any_banned_percent = (filter_results.full_pass == False).mean()
        all_banned_percent = (filter_results.any_pass  == False).mean()
        
        res_df = pd.DataFrame({"model_name" : [model_name],
                               "num_mol" : [num_mol],
                               "all_banned_ratio" : [all_banned_percent],
                               "any_banned_ratio" : [any_banned_percent]
                             })
        for name in alert_names:
            res_df[f"{name}_banned_ratio"] = 1 - filter_results[f"pass_filter_{name}"].mean()
        
        res_df, filter_extended = common_postprocessing_statistics(filter_results, res_df, stat, extend)
        return res_df, filter_extended


    elif filter_name == 'bredt':
        res_df = pd.DataFrame({"model_name" : [model_name],
                               "num_mol" : [num_mol],
                              })
        res_df[f"banned_ratio"] = 1 - filter_results[f"full_pass"].mean()

        res_df, filter_extended = common_postprocessing_statistics(filter_results, res_df, stat, extend)
        return res_df, filter_extended
    

    elif filter_name == 'NIBR':
        mean_severity = filter_results.severity.mean()
        max_severity = filter_results.severity.max()

        mean_n_covalent_motif = filter_results.n_covalent_motif.mean()
        mean_nonzero_special_mol = (filter_results.special_mol > 0).mean()

        res_df = pd.DataFrame({"model_name" : [model_name],
                               "num_mol" : [num_mol],
                               "mean_severity" : [mean_severity],
                               "max_severity" : [max_severity],
                               "mean_n_covalent_motif" : [mean_n_covalent_motif],
                               "mean_nonzero_special_mol" : [mean_nonzero_special_mol]
                             })
        res_df[f"banned_ratio"] = 1 - filter_results[f"pass_filter"].mean()
        
        res_df, filter_extended = common_postprocessing_statistics(filter_results, res_df, stat, extend)
        filter_extended.rename(columns={'pass_filter' : 'full_pass'}, inplace=True)
        return res_df, filter_extended


    elif filter_name == 'molgraph_stats':
        res_df = pd.DataFrame({"model_name" : [model_name],
                               "num_mol" : [num_mol],
                              })
        for i in range(1, 12):
            res_df[f"banned_ratio_s_{i}"] = 1 - filter_results[f"pass_filter_{i}"].mean()

        res_df, filter_extended = common_postprocessing_statistics(filter_results, res_df, stat, extend)
        pass_filter_cols = [col for col in filter_extended.columns if 'pass_filter' in col]
        filter_extended['full_pass'] = filter_extended[pass_filter_cols].all(axis=1)
        return res_df, filter_extended


    elif filter_name == 'molcomplexity':
        any_banned_percent = 1 - filter_results.full_pass.mean()
        all_banned_percent = 1 - filter_results.pass_any.mean()

        res_df = pd.DataFrame({"model_name" : [model_name],
                                "num_mol" : [num_mol],
                                "all_banned_ratio" : [all_banned_percent],
                                "any_banned_ratio" : [any_banned_percent]
                             })
        alert_names = mc.complexity.ComplexityFilter.list_default_available_filters()
        for name in alert_names:
            res_df[f"{name}_banned_ratio"] = 1 - filter_results[f"pass_{name}"].mean()
        
        res_df, filter_extended = common_postprocessing_statistics(filter_results, res_df, stat, extend)
        return res_df, filter_extended
        
    else:
        raise ValueError(f"Filter {filter_name} not found")
    

def check_paths(paths):
    required_patterns = ['nibr', 'bredt', 'molgraphstats', 'molcomplexity', 'commonalerts']
    missing_patterns = [pattern for pattern in required_patterns if not any(pattern in path.lower() for path in paths)]
    if len(missing_patterns) > 0:
        raise AssertionError(f"Invalid filter name(s) missing: {', '.join(missing_patterns)}")
    return True


def plot_calculated_stats(config, prefix):
    mode = config['mode']
    folder_to_save = process_path(config['folder_to_save'])

    if prefix == 'beforeDescriptors':
        paths = glob.glob(folder_to_save + f'{prefix}_StructFilters/' + '*metrics.csv')
    else:
        paths = glob.glob(folder_to_save + f'StructFilters/' + '*metrics.csv')
    check_paths(paths)

    model_name_set = get_model_name(config, mode)
    datas=[]
    filter_names = []

    for path in paths:
        data = pd.read_csv(path)
        data.set_index('model_name', inplace=True)

        banned_cols = [col for col in data.columns if 'banned_ratio' in col]
        data_filtered = data[banned_cols + ['num_mol']].copy()
        for banned_col in banned_cols:
            data_filtered.loc[:, f'num_banned_{banned_col}'] = data_filtered[banned_col] * data_filtered['num_mol']
        datas.append(data_filtered)
        filter_name = path.split(f'{model_name_set}/')[-1].split('_metrics.csv')[0]
        filter_names.append(filter_name)
    
    filter_results = {}
    if prefix == 'beforeDescriptors':
        filters_to_find = glob.glob(folder_to_save + f'{prefix}_StructFilters/*filteredMols.csv')
    else:
        filters_to_find = glob.glob(folder_to_save + f'StructFilters/*filteredMols.csv')
    
    for path in filters_to_find:
        try:
            filter_data = pd.read_csv(path)
            model_name = filter_data['model_name']
            filter_name = path.split('/')[-1].split('filteredMols.csv')[0].strip('_')
            if 'full_pass' in filter_data.columns:
                passed = filter_data[filter_data['full_pass'] == True]
                num_passed_by_model = passed.groupby('model_name').size().to_dict()

            if num_passed_by_model is not None:
                filter_results[filter_name] = num_passed_by_model
            elif num_passed_by_model is None:
                filter_results[filter_name] = {model_name: 0}

        except (IndexError, FileNotFoundError):
            filter_results[filter_name] = {model_name: 0}
    for filter, values in filter_results.items():
        if len(values) != len(model_name):
            for model in model_name:
                if model not in values.keys():
                    filter_results[filter][model] = 0        
    for filter_name, models in filter_results.items():
        filter_results[filter_name] = dict(sorted(models.items()))

    n_plots = len(datas)
    n_cols = 2
    n_rows = (n_plots + n_cols - 1) // n_cols  
    
    fig = plt.figure(figsize=(40, 5*n_rows)) 
    for idx, (data, filter_name) in enumerate(zip(datas, filter_names)):
        ax = plt.subplot(n_rows, n_cols, idx + 1)
        models = data.index
        x = np.arange(len(models))
        width = 0.8 
        total_mols = data['num_mol'].sum()
        total = ax.barh(x, data.loc[models, 'num_mol'], width, label=f'Total Molecules ({format_number(total_mols)})', color='#E5E5E5', alpha=0.5)

        clean_filter_name = filter_name.split('/')[-1].lower()
        for known_filter in filter_results.keys():
            if known_filter.lower() in clean_filter_name:
                if mode == 'single_comparison':
                    passed = list(filter_results[known_filter].values())[0]
                    bar_center_x = data.loc[models, 'num_mol'].values[0] / 2
                    bar_center_y = x[0]
                    model_total = data.loc[model, 'num_mol']
                    text = f'Passed molecules: {passed} ({(passed/model_total*100):.1f}%)'
                    ax.annotate(text, (bar_center_x, bar_center_y), ha='center', va='center', fontsize=12, color='black', fontweight='bold',
                                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=3), zorder=1000)
                else:
                    for i, (model, passed) in enumerate(filter_results[known_filter].items()):
                        bar_center_x = data.loc[models, 'num_mol'].values[0] / 2
                        bar_center_y = x[i]
                        model_total = data.loc[model, 'num_mol']
                        if model_total != 0:
                            text = f'Passed molecules: {passed} ({(passed/model_total*100):.1f}%)'
                        else:
                            text = f'Passed molecules: {passed} (0%)'
                        ax.annotate(text, (bar_center_x, bar_center_y), ha='center', va='center', fontsize=12, color='black', fontweight='bold',
                                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=3), zorder=1000)


        for i, model in enumerate(models):
            count = data.loc[model, 'num_mol']
            max_bar_width = data['num_mol'].max()
            text_x_position = max_bar_width 
            ax.text(text_x_position, i, int(count),  va='center', ha='left', fontsize=12, color='black', fontweight='bold')     
        
        banned_bars = []
        banned_percentages = [] 
        ratio_cols = [col for col in data.columns if 'banned_ratio' in col and 'num_banned' not in col]
        colors = get_model_colors(model_names=ratio_cols, cmap='Paired')
        for col, color in zip(ratio_cols, colors.values()):
            num_banned_col = f'num_banned_{col}'
            ratio_name = col.replace('banned_ratio', '').strip('_')
            ratio_name = clean_name(ratio_name)
            
            banned_count = data[num_banned_col]
            total_banned = banned_count.sum()
            banned_percent = (total_banned / total_mols) * 100 if total_mols > 0 else 0
            banned_percentages.append(banned_percent)
            
            if banned_percent == 0.0:
                label = f'{ratio_name} (0%)'
            else:
                label = f'{ratio_name} ({format_number(total_banned)}, {banned_percent:.1f}%)'
            
            bar = ax.barh(x, banned_count, width, label=label, color=color, alpha=0.8)
            banned_bars.append(bar)
        
        clean_filter_name = filter_name.split('/')[-1] 
        ax.set_title(clean_name(clean_filter_name), fontsize=14, pad=20, fontweight='bold')
        
        ax.set_yticks(x)
        ax.set_yticklabels(models, fontsize=12)
        
        ax.xaxis.set_major_formatter(plt.FuncFormatter(format_number))

        # max_mols = int(data['num_mol'].max())
        # step = max(1, max_mols // 5) 
        # ticks = list(range(0, max_mols + 1, step))
        # ax.set_xticks(ticks)
        # ax.set_xticklabels(ticks)
        ax.set_xlim(left=0)

        if mode == 'single_comparison':
            ax.set_xlabel(f'Number of Molecules (Total: {format_number(total_mols)})', fontsize=12, labelpad=10)
        ax.set_ylabel('Models', fontsize=12, labelpad=10)
        
        ax.grid(True, axis='x', alpha=0.2, linestyle='--')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        sorted_indices = np.argsort(banned_percentages)[::-1] 
        sorted_handles = [total] + [banned_bars[i] for i in sorted_indices]
        
        legend = ax.legend(handles=sorted_handles, loc='center left',  bbox_to_anchor=(1.02, 0.5), fontsize=11, ncol=1)
        legend.get_frame().set_alpha(0.9)
        legend.get_frame().set_edgecolor('lightgray')

    plt.subplots_adjust(right=0.85, hspace=0.6, wspace=0.5) 
    
    if prefix == 'beforeDescriptors':
        plt.savefig(folder_to_save + f'{prefix}_StructFilters/' + f'MoleculeCountsComparison.png', dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    else:
        plt.savefig(folder_to_save + f'StructFilters/' + f'MoleculeCountsComparison.png', dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()

    return


def plot_restriction_ratios(config, prefix):
    mode = config['mode']
    folder_to_save = process_path(config['folder_to_save'])
    folder_name = config['folder_to_save'].split('/')[-1]

    if prefix == 'beforeDescriptors':
        paths = glob.glob(folder_to_save + f'{prefix}_StructFilters/' + '*metrics.csv')
    else:
        paths = glob.glob(folder_to_save + f'StructFilters/' + '*metrics.csv')
    check_paths(paths)
    
    filter_data = {}
    
    if not paths:
        logger.error(f"No data files found in {folder_to_save}")
        return
    
    for path in paths:
        filter_name = path.split(f'{folder_name}/')[-1].split('_metrics.csv')[0]
        data = pd.read_csv(path)
        ratio_cols = [col for col in data.columns if 'banned_ratio' in col]
        
        if not ratio_cols:
            continue
            
        clean_cols = {col: col.replace('_banned_ratio', '').replace('banned_ratio', '').replace('_s', 's') for col in ratio_cols}
        ratios = data[ratio_cols].rename(columns=clean_cols)

        actual_model_names = data['model_name'].tolist()
        ratios.index = actual_model_names
        
        row = ratios.iloc[0]
        if row.isna().all():
            continue
            
        all_value = None
        if 'all' in row.index:
            all_value = row['all']
            row = row.drop('all')
        
        sorted_values = row.sort_values(ascending=False)
        
        if all_value is not None:
            if all_value >= sorted_values.iloc[0]:
                sorted_index = pd.Index(['all']).append(sorted_values.index)
            else:
                sorted_index = sorted_values.index.append(pd.Index(['all']))
            ratios = ratios[sorted_index]
        else:
            ratios = ratios[sorted_values.index]
        
        filter_data[filter_name] = ratios

    if not filter_data:
        logger.error("No valid data to plot")
        return

    plt.style.use('default')
    sns.set_style("white")
    sns.set_context("talk")

    n_filters = len(filter_data)
    n_cols = min(2, n_filters) 
    n_rows = (n_filters + n_cols - 1) // n_cols

    fig = plt.figure(figsize=(16, 7*n_rows))
    fig.suptitle('Comparison of Restriction Ratios Across Different Filters', fontsize=16, y=0.98, fontweight='bold')


    for idx, (filter_name, data) in enumerate(filter_data.items()):
        ax = plt.subplot(n_rows, n_cols, idx + 1)
        
        if not data.empty and data.notna().any().any():
            sns.heatmap(data.T, cmap='Blues', 
                       cbar_kws={'label': 'Restriction Ratio'}, ax=ax, vmin=0, vmax=1,
                       fmt='.3f', annot=True, annot_kws={'size': 10, 'rotation': 0}, cbar=True)
            
            ax.set_title(f'{clean_name(filter_name)} Filter', pad=10, fontsize=11, fontweight='bold')
            plt.setp(ax.get_yticklabels(), rotation=0, ha='right')
            plt.setp(ax.get_xticklabels(), ha='center')
            ax.set_xlabel('Model')
            ax.set_ylabel('Filter Type')

            actual_model_names = data.index.tolist()
            if len(actual_model_names) == len(ax.get_xticklabels()):
                ax.set_xticklabels(actual_model_names)
            
        else:
            ax.text(0.5, 0.5, 'No data available', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
            ax.set_title(f'{clean_name(filter_name)} Filter', pad=10, fontsize=11, fontweight='bold')

    plt.tight_layout()
    if prefix == 'beforeDescriptors':
        plt.savefig(folder_to_save + f'{prefix}_StructFilters/' + 'RestrictionRatiosComparison.png', dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    else:
        plt.savefig(folder_to_save + f'StructFilters/' + 'RestrictionRatiosComparison.png', dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()
    return



def filter_data(config, prefix):
    mode = config['mode']
    if prefix == 'beforeDescriptors':
        folder_to_save = process_path(config['folder_to_save'], key_word=f'{prefix}_StructFilters/')
    else:
        folder_to_save = process_path(config['folder_to_save'], key_word='StructFilters')
    paths = glob.glob(folder_to_save + '*filteredMols.csv')
    check_paths(paths)

    columns_to_drop = ['full_pass', 'any_pass', 'name', 'pass_any']
    datas = []
    for path in paths:
        data = pd.read_csv(path)
        if mode == 'single_comparison':
            model_name = get_model_name(config, mode)
            filter_name = path.split(f'{model_name}/')[-1].split('_filteredMols.csv')[0].strip('_')
        else:
            filter_name = path.split(f'{folder_to_save}')[-1].split('_filteredMols.csv')[0].strip('_')

        for col in columns_to_drop:
            if col in data.columns:
                data.drop(columns=[col], inplace=True)
        
        data = data.rename(columns={col: f"{filter_name}_{col}" for col in data.columns if col != 'SMILES'})
        datas.append(data)

    filtered_data = reduce(lambda x, y: pd.merge(x, y, on='SMILES', how='inner'), datas)
    filtered_data['model_name'] = filtered_data['Nibr_model_name']

    filtered_data.to_csv(folder_to_save + 'leftMolsAfterStructFiltersMetrics.csv', index=False)
    if mode == 'single_comparison':
        filtered_data['SMILES'].to_csv(folder_to_save + 'leftMolsAfterStructFiltersSMILES.csv', index=False)
    else:
        filtered_data[['SMILES', 'model_name']].to_csv(folder_to_save + 'leftMolsAfterStructFiltersSMILES.csv', index_label='SMILES', index=False)

    return filtered_data
