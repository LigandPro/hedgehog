import os 
import glob 
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from rdkit import Chem
from tqdm.auto import tqdm
from functools import reduce
from pandarallel import pandarallel

import datamol as dm
import medchem as mc
from medchem.structural.lilly_demerits import LillyDemeritsFilters

from logger_config import logger
from configs.config_utils import load_config


def camelcase(any_str):
    camel_str = ''.join(word.capitalize() for word in any_str.split('_'))
    return camel_str


def build_identity_map_from_descriptors(config):
    base_folder = process_path(config['folder_to_save'])
    id_path = base_folder + 'Descriptors/passDescriptorsSMILES.csv'
    identity_map = {}
    id_df = None
    try:
        if os.path.exists(id_path):
            id_df = pd.read_csv(id_path)
            if 'smiles' not in id_df.columns:
                smiles_col = next((c for c in id_df.columns if c.lower() == 'smiles'), None)
                if smiles_col:
                    id_df = id_df.rename(columns={smiles_col: 'smiles'})
            if 'model_name' not in id_df.columns:
                id_df['model_name'] = 'single'
            if 'mol_idx' in id_df.columns:
                for _, row in id_df.iterrows():
                    key_raw = (row['smiles'], row['model_name'])
                    identity_map[key_raw] = row['mol_idx']
    except Exception:
        pass
    return identity_map, id_df



def get_model_name(config):
    paths = glob.glob(config['generated_mols_path'])
    if len(paths) > 1:
        return [path.split('/')[-1].split('.')[0] for path in paths]
    single = paths[0]
    try:
        tmp = pd.read_csv(single)
        lower_cols = {c.lower(): c for c in tmp.columns}
        candidate = lower_cols.get('model_name') or lower_cols.get('name')
        if candidate is not None and tmp[candidate].nunique(dropna=True) > 1:
            return sorted(tmp[candidate].dropna().unique().tolist())
    except Exception:
        pass
    return single.split('/')[-1].split('.')[0]


def process_path(folder_to_save, key_word=None):
    if not folder_to_save.endswith('/'):
        folder_to_save = folder_to_save + '/'

    if key_word:
        folder_to_save = folder_to_save + f'{key_word}/'

    os.makedirs(folder_to_save, exist_ok=True)
    return folder_to_save


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
    

def filter_alerts(config):
    df = pd.read_csv(config["alerts_data_path"])
    mask = df.rule_set_name.isin(config["include_rulesets"])
    for ruleset in config["exclude_descriptions"].keys():
        local_mask_1 = ~(df.rule_set_name == ruleset)
        local_mask_2 = ~df.description.isin(config["exclude_descriptions"][ruleset])
        local_mask = local_mask_1 | local_mask_2
        mask &= local_mask
    return df[mask]
    

def process_path(folder_to_save, key_word=None):
    if not folder_to_save.endswith('/'):
        folder_to_save = folder_to_save + '/'

    if key_word:
        folder_to_save = folder_to_save + f'{key_word}/'
    os.makedirs(folder_to_save, exist_ok=True)
    return folder_to_save


def common_postprocessing_statistics(filter_results, res_df, stat, extend):          
    if stat is not None:
        res_df = pd.concat([stat, res_df])

    # Preserve existing SMILES if present; only derive from mol when missing
    if 'smiles' not in filter_results.columns:
        filter_results["smiles"] = filter_results.mol.apply(lambda m: dm.to_smiles(m) if m is not None else None)
    filter_results = filter_results.drop(columns='mol')

    if extend is not None:
        filter_extended = pd.concat([extend, filter_results]).copy()
    else:
        filter_extended = filter_results.copy()

    return res_df,filter_extended


def process_one_file(config, input_path, apply_filter, subsample):
    input_type = input_path[input_path.rfind(".")+1:]
    assert input_type in {"csv", "smi", 'sdf', 'txt'}

    if input_type == 'csv':
        data = pd.read_csv(input_path)
        smiles_col = None
        if 'smiles' in data.columns: smiles_col = 'smiles'
        elif 'smiles' in data.columns: smiles_col = 'smiles'

        if 'smiles' not in data.columns:
            for col in data.columns:
                if "smiles" in col.lower():
                    smiles_col = col
                    break

        assert smiles_col is not None

        has_mol_idx = 'mol_idx' in data.columns
        is_multi = ('model_name' in data.columns and data['model_name'].nunique(dropna=True) > 1)
        if is_multi:
            smiles_str = data[smiles_col].tolist()
            model_names = data['model_name'].tolist()
            mols = [dm.to_mol(x) for x in smiles_str]
            if has_mol_idx:
                mol_indices = data['mol_idx'].tolist()
                smiles = list(zip(smiles_str, model_names, mols, mol_indices))
            else:
                smiles = list(zip(smiles_str, model_names, mols))
        else:
            smiles = data[smiles_col].tolist()
            mols = [dm.to_mol(x) for x in smiles]
            if has_mol_idx:
                mol_indices = data['mol_idx'].tolist()
                smiles = list(zip(smiles, [None]*len(smiles), mols, mol_indices))

    elif input_type == "smi" or input_type == 'txt':
        with open(input_path, 'r') as file:
            lines = [line.rstrip('\n') for line in file]

        if subsample <= len(lines):
            lines = np.random.permutation(lines)[:subsample].tolist()

        smiles = []
        model_names = []
        for line in lines:
            parts = line.split(',')
            if len(parts) == 2:
                smi, model = parts
                smiles.append(smi)
                model_names.append(model)
            else:
                smiles.append(parts[0])
        mols = [dm.to_mol(x) for x in smiles]
        if len(model_names) == len(smiles):
            smiles = list(zip(smiles, model_names, mols))

    elif input_type == 'sdf':
        mols, smiles = sdf_to_mols(input_path, subsample)

    if isinstance(smiles[0], tuple):
        cleaned = []
        for item in smiles:
            if len(item) >= 3:
                smi, model, mol = item[0], item[1], item[2]
                if mol is not None:
                    if len(item) >= 4:
                        mol_idx = item[3]
                        cleaned.append((smi, model, mol, mol_idx))
                    else:
                        cleaned.append((smi, model, mol))
        smiles = cleaned
        mols = [it[2] for it in smiles]
    else:
        mols, smiles = dropna(mols, smiles)

    assert len(mols) == len(smiles), f"{len(mols)}, {len(smiles)}"
    if not isinstance(smiles[0], tuple):
        assert len(mols) <= subsample

    for mol, smi in zip(mols, smiles):
        smi_val = smi[0] if isinstance(smi, tuple) else smi
        assert mol is not None, f"{smi_val}"

    final_result = None

    if len(mols) > 0:
        if isinstance(smiles[0], tuple):
            final_result = apply_filter(config, mols, smiles)
        else:
            final_result = apply_filter(config, mols)

    return final_result


def add_model_name_col(final_result, smiles_with_model):
    if len(final_result) == len(smiles_with_model):
        if isinstance(smiles_with_model[0], tuple) and len(smiles_with_model[0]) >= 3:
            smiles_vals = []
            model_vals = []
            mol_idx_vals = []
            any_mol_idx = False
            for item in smiles_with_model:
                smi = item[0]
                model = item[1]
                smiles_vals.append(smi)
                model_vals.append(model if model is not None else 'single')
                if len(item) >= 4:
                    mol_idx_vals.append(item[3])
                    any_mol_idx = True
                else:
                    mol_idx_vals.append(None)
            final_result['smiles'] = smiles_vals
            final_result['model_name'] = model_vals
            if any_mol_idx:
                final_result['mol_idx'] = mol_idx_vals
        else:
            final_result['smiles'] = smiles_with_model
    else:
        logger.warning("Length mismatch between results and smiles_with_model! Attempting to align valid mols.")
        valid_records = []
        mols_in_result = list(final_result['mol'])
        for mol in mols_in_result:
            matched = None
            for item in smiles_with_model:
                if isinstance(item, tuple) and len(item) >= 3:
                    smi, model, orig_mol = item[0], (item[1] if item[1] is not None else 'single'), item[2]
                    mol_idx = item[3] if len(item) >= 4 else None
                else:
                    smi, model, orig_mol, mol_idx = item, 'single', None, None
                if orig_mol is not None and mol == orig_mol:
                    matched = (smi, model, mol_idx)
                    break
            if matched is not None:
                valid_records.append(matched)
        if len(valid_records) == len(final_result) and len(valid_records) > 0:
            s_list, m_list, idx_list = zip(*valid_records)
            final_result['smiles'] = s_list
            final_result['model_name'] = m_list
            if any(x is not None for x in idx_list):
                final_result['mol_idx'] = idx_list
        else:
            logger.error("Could not align all mols with smiles/model_name. Check your data integrity.")

    return final_result


def filter_function_applier(filter_name):
    if filter_name == 'common_alerts':
        return apply_structural_alerts
    elif filter_name == 'molgraph_stats':
        return apply_molgraph_stats
    elif filter_name == 'molcomplexity':
        return apply_molcomplexity_filters
    elif filter_name == 'NIBR':
        return apply_nibr_filter
    elif filter_name == 'bredt':
        return apply_bredt_filter
    elif filter_name == 'lilly':
        return apply_lilly_filter
    else:
        raise ValueError(f"Filter {filter_name} not found")


def apply_structural_alerts(config, mols, smiles_modelName_mols=None):
    def _apply_alerts(row):
        mol = row["mol"]
        row["smiles"] = dm.to_smiles(mol) if mol is not None else None

        config_structFilters = load_config(config['config_structFilters'])
        alert_data = filter_alerts(config_structFilters)
        df = alert_data.copy()
        df["matches"] = df.smarts.apply(lambda x, y: y.GetSubstructMatches(Chem.MolFromSmarts(x)), args=(mol,)) 
        grouped = df.groupby("rule_set_name").apply(
            lambda group: pd.Series({"matches": [match for matches in group["matches"] 
                                                       for match in matches] 
                                                       if any(matches 
                                                              for matches in group["matches"]) 
                                                       else (),
                                     "reasons": ";".join(group[group["matches"].apply(lambda x: len(x) > 0)]["description"].fillna('').tolist()) 
                                                                if any(matches 
                                                                       for matches in group["matches"]) 
                                                                else "",
                                    }),
            include_groups=False
        ).reset_index()
        grouped["pass_filter"] = grouped["matches"].apply(lambda x: True if not x else False)
        for _, g_row in grouped.iterrows():
            name = g_row["rule_set_name"]
            row[f"pass_filter_{name}"] = g_row["pass_filter"]
            row[f"reasons_{name}"] = g_row["reasons"]
        return row


    def _get_full_any_pass(row):
        full_pass = True
        any_pass  = False
        for col in row.index:
            if "pass_filter" in col:
                full_pass &= row[col]
                any_pass |= row[col]
        row["full_pass"] = full_pass
        row["any_pass"] = any_pass
        return row
    

    logger.info(f"Processing {len(mols)} filtered molecules")

    n_jobs = config['n_jobs']

    pandarallel.initialize(progress_bar=True, nb_workers=n_jobs)

    mols_df = pd.DataFrame({"mol" : mols})
    results = mols_df.parallel_apply(_apply_alerts, axis=1)
    results = results.parallel_apply(_get_full_any_pass, axis=1)

    if smiles_modelName_mols is not None:
        results = add_model_name_col(results, smiles_modelName_mols)
    return results


def apply_molgraph_stats(config, mols: list[Chem.rdchem.Mol], smiles_modelName_mols=None):
    logger.info(f"Processing {len(mols)} molecules to calculate Molecular Graph statistics")
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
        results = add_model_name_col(results, smiles_modelName_mols)
    return results


def apply_molcomplexity_filters(config, mols: list[Chem.Mol], smiles_modelName_mols=None):
    logger.info(f"Processing {len(mols)} molecules to calculate Complexity filters")
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
        final_result = add_model_name_col(final_result, smiles_modelName_mols)
    return final_result


def apply_bredt_filter(config, mols: list[Chem.Mol], smiles_modelName_mols=None):
    logger.info(f"Processing {len(mols)} molecules to calculate Bredt filter")
    out = mc.functional.bredt_filter(mols=mols,
                                     n_jobs=-1,
                                     progress=False,
                                     return_idx=False,
                                    )
    results = pd.DataFrame({'mol' : mols,
                            'full_pass' : out 
                            })    
    if smiles_modelName_mols is not None:
        results = add_model_name_col(results, smiles_modelName_mols) 
    return results


def apply_nibr_filter(config, mols: list[Chem.Mol], smiles_modelName_mols=None):
    logger.info(f"Processing {len(mols)} molecules to calculate NIBR filter")
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
        results = add_model_name_col(results, smiles_modelName_mols)      
    return results


def apply_lilly_filter(config, mols, smiles_modelName_mols=None):
    logger.info(f"Processing {len(mols)} molecules to calculate Lilly filter")
    n_jobs = config['n_jobs']

    config_strcuFilters = load_config(config['config_structFilters'])
    scheduler = config_strcuFilters['lilly_scheduler']

    dfilter = LillyDemeritsFilters()
    results = dfilter(mols=mols,
                      n_jobs=n_jobs,
                      scheduler=scheduler,
                      )
        
    if smiles_modelName_mols is not None:
        results = add_model_name_col(results, smiles_modelName_mols)
    return results 


def get_basic_stats(config, filter_results: pd.DataFrame, model_name:str, filter_name:str, stat=None, extend=None):
    is_multi = 'model_name' in filter_results.columns and filter_results['model_name'].nunique(dropna=True) > 1
    if is_multi:
        all_res = []
        all_extended = []
        for model, group in filter_results.groupby('model_name'):
            res_df, filter_extended = get_basic_stats(config, group.copy(), model, filter_name, stat, extend)
            all_res.append(res_df)
            all_extended.append(filter_extended)

        res_df = pd.concat(all_res, ignore_index=True)
        filter_extended = pd.concat(all_extended, ignore_index=True)
        return res_df, filter_extended

    num_mol = len(filter_results)

    filter_results.dropna(subset='mol', inplace=True)
    filter_results['model_name'] = model_name 

    if filter_name == 'common_alerts':
        any_banned_percent = (filter_results.full_pass == False).mean()
        all_banned_percent = (filter_results.any_pass  == False).mean()

        res_df = pd.DataFrame({"model_name" : [model_name],
                               "num_mol" : [num_mol],
                               "all_banned_ratio" : [all_banned_percent],
                               "any_banned_ratio" : [any_banned_percent]
                             })

        for name in config["include_rulesets"]:
            res_df[f"{name}_banned_ratio"] = 1 - filter_results[f"pass_filter_{name}"].mean()

        res_df, filter_extended = common_postprocessing_statistics(filter_results, res_df, stat, extend)
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


    elif filter_name == 'lilly':
        mean_noNA_demerit_score = filter_results.demerit_score.dropna().mean()
        
        res_df = pd.DataFrame({'model_name': [model_name], 
                               'num_mol': [num_mol],
                               'mean_noNA_demerit_score': mean_noNA_demerit_score
                             })
        res_df[f"banned_ratio"] = 1 - filter_results[f"pass_filter"].mean()
        
        res_df, filter_extended = common_postprocessing_statistics(filter_results, res_df, stat, extend)
        filter_extended.rename(columns={'pass_filter' : 'full_pass'}, inplace=True)
        return res_df, filter_extended
        
    else:
        raise ValueError(f"Filter {filter_name} not found")
    

def check_paths(config, paths):
    all_filters = {}
    for k, v in config.items():
        if 'calculate_' in k:
            k = k.replace('calculate_', '')
            all_filters[k] = v

    required_patterns = [''.join(k.split('_')) for k, v in all_filters.items() if v]
    missing_patterns = [pattern 
                        for pattern in required_patterns 
                        if not any(pattern.lower() in path.lower() 
                            for path in paths)]
    if len(missing_patterns) > 0:
        raise AssertionError(f"Invalid filter name(s) missing: {', '.join(missing_patterns)}")
    return True


def plot_calculated_stats(config, prefix):
    folder_to_save = process_path(config['folder_to_save'])

    config_structFilters = load_config(config['config_structFilters'])

    if prefix == 'beforeDescriptors':
        paths = glob.glob(folder_to_save + f'{prefix}_StructFilters/' + '*metrics.csv')
    else:
        paths = glob.glob(folder_to_save + f'StructFilters/' + '*metrics.csv')
    check_paths(config_structFilters, paths)

    model_name_set = get_model_name(config)
    datas=[]
    filter_names = []

    identity_map, src = build_identity_map_from_descriptors(config)

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
                for i, (model, passed) in enumerate(filter_results[known_filter].items()):
                    bar_center_x = data.loc[models, 'num_mol'].values[0] / 2
                    bar_center_y = x[i]
                    model_total = data.loc[model, 'num_mol'] if model in data.index else data['num_mol'].iloc[0]
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
        ax.set_xlim(left=0)

    if isinstance(model_name, str):
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
    folder_to_save = process_path(config['folder_to_save'])
    folder_name = config['folder_to_save'].split('/')[-1]

    config_structFilters = load_config(config['config_structFilters'])

    if prefix == 'beforeDescriptors':
        paths = glob.glob(folder_to_save + f'{prefix}_StructFilters/' + '*metrics.csv')
    else:
        paths = glob.glob(folder_to_save + f'StructFilters/' + '*metrics.csv')
    check_paths(config_structFilters, paths)
    
    if not paths:
        logger.error(f"No data files found in {folder_to_save}")
        return
    
    filter_data = {}
    model_names_filters = {}
    for path in paths:
        filter_name = path.split(f'{folder_name}/')[-1].split('_metrics.csv')[0]
        data = pd.read_csv(path)

        ratio_cols = [col for col in data.columns if 'banned_ratio' in col]
        model_names_filters[filter_name] = {}
        model_names_filters[filter_name] = dict(zip(data['model_name'].tolist(), data['num_mol'].tolist()))

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

    model_names_filters = pd.DataFrame(model_names_filters).reset_index()
    model_names_filters = model_names_filters.rename(columns={'index': 'model_name'})

    plt.style.use('default')
    sns.set_style("white")
    sns.set_context("talk")

    n_filters = len(filter_data)
    n_cols = min(2, n_filters) 
    n_rows = (n_filters + n_cols - 1) // n_cols

    fig = plt.figure(figsize=(16, 7*n_rows))
    fig.suptitle('Comparison of Restriction Ratios Across Different Filters', fontsize=16, y=0.98, fontweight='bold')

    for idx, (filter_name, data) in enumerate(filter_data.items()):
        number_of_mols = np.array(model_names_filters[filter_name].tolist())
        ax = plt.subplot(n_rows, n_cols, idx + 1)
        for col in data.columns:
            if col not in ['any', 'all', 'model_name', 'num_mol']:
                data[col] = number_of_mols * (1 - np.array(data[col].tolist()))
        if not data.empty and data.notna().any().any():
            if 'all' in data.columns:
                data.drop(columns=['all'], inplace=True) 
            if 'any' in data.columns:
                data.drop(columns=['any'], inplace=True)

            from matplotlib.colors import LinearSegmentedColormap
            custom_cmap = LinearSegmentedColormap.from_list('custom', ['white', '#b29eee'])
            if idx == 1:

                sns.heatmap(data.T, cmap=custom_cmap, 
                        cbar_kws={'label': 'Passed Molecules', 'format': '%d', }, ax=ax, vmin=0, vmax=max(data.max()),
                        fmt='.0f', 
                        annot=True, annot_kws={'size': 12, 'rotation': 0, 'color': 'black'}, cbar=True)
            else:
                sns.heatmap(data.T, cmap=custom_cmap, 
                        cbar_kws={'label': 'Passed Molecules', 'format': '%d', }, ax=ax, vmin=0, vmax=max(data.max()),
                        fmt='.0f', 
                        annot=True, annot_kws={'size': 12, 'rotation': 0, 'color': 'black'}, cbar=False)
            
            ax.set_title(f'{clean_name(filter_name)} Filter', fontsize=12, fontweight='bold')
            plt.setp(ax.get_yticklabels(), rotation=0, ha='right', fontsize=12)
            plt.setp(ax.get_xticklabels(), rotation=0, ha='right', fontsize=12)
            ax.set_xlabel('Model')

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
    config_structFilters = load_config(config['config_structFilters'])

    if prefix == 'beforeDescriptors':
        folder_to_save = process_path(config['folder_to_save'], key_word=f'{prefix}_StructFilters')
    else:
        folder_to_save = process_path(config['folder_to_save'], key_word='StructFilters')
    paths = glob.glob(folder_to_save + '*filteredMols.csv')

    columns_to_drop = ['full_pass', 'any_pass', 'name', 'pass_any']
    datas = []
    for path in paths:
        data = pd.read_csv(path)
        filter_name = path.split(f'{folder_to_save}')[-1].split('_filteredMols.csv')[0].strip('_')

        for col in columns_to_drop:
            if col in data.columns:
                data.drop(columns=[col], inplace=True)

        data = data.rename(columns={col: f"{filter_name}_{col}" for col in data.columns if col not in ['smiles', 'model_name']})

        if 'model_name' not in data.columns:
            data['model_name'] = 'single'
        datas.append(data)

    identity_map, _ = build_identity_map_from_descriptors(config)

    if len(datas) > 0:
        filtered_data = reduce(lambda x, y: pd.merge(x, y,
                                                    on=['smiles', 'model_name'],
                                                    how='inner'), datas)
    else:
        filtered_data = pd.DataFrame(columns=['smiles', 'model_name'])
    model_name_cols = [col for col in filtered_data.columns if 'model_name' in col.lower()]
    if model_name_cols:
        filtered_data['model_name'] = filtered_data[model_name_cols[0]]

    if identity_map and len(filtered_data) > 0:
        def _canon(s):
            return str(s)
        raw_vals = [identity_map.get((s, m)) for s, m in zip(filtered_data['smiles'], filtered_data['model_name'])]
        canon_vals = [identity_map.get((_canon(s), m)) for s, m in zip(filtered_data['smiles'], filtered_data['model_name'])]
        filled = []
        for r, c in zip(raw_vals, canon_vals):
            filled.append(r if pd.notna(r) and r != '' else c)
        filtered_data['mol_idx'] = filled
    elif identity_map and len(filtered_data) == 0:
        id_df = pd.DataFrame(list(identity_map.items()), columns=['key', 'mol_idx'])
        id_df[['smiles', 'model_name']] = pd.DataFrame(id_df['key'].tolist(), index=id_df.index)
        filtered_data = id_df[['smiles', 'model_name', 'mol_idx']]

    if 'model_name' not in filtered_data.columns:
        filtered_data['model_name'] = 'single'
    cols = ['smiles', 'model_name'] + (['mol_idx'] if 'mol_idx' in filtered_data.columns else [])
    out_df = filtered_data[cols].copy()
    out_df = out_df[['smiles', 'model_name'] + ([ 'mol_idx'] if 'mol_idx' in out_df.columns else [])]
    out_df.to_csv(folder_to_save + 'passStructFiltersSMILES.csv', index_label='smiles', index=False)

    extended_paths = glob.glob(folder_to_save + '*extended.csv')
    if extended_paths:
        failures_map = {}
        has_model_name = False
        for path in extended_paths:
            try:
                df_ext = pd.read_csv(path)
                filter_name = path.split('/')[-1].split('_extended.csv')[0].strip('_')
                if 'model_name' in df_ext.columns:
                    has_model_name = True
                if 'full_pass' in df_ext.columns:
                    failed_df = df_ext[df_ext['full_pass'] == False]
                    for _, row in failed_df.iterrows():
                        key = (row['smiles'], row['model_name']) if 'model_name' in failed_df.columns else (row['smiles'], None)
                        if key not in failures_map:
                            failures_map[key] = set()
                        failures_map[key].add(filter_name)
            except Exception:
                continue

        if failures_map:
            def _canon(s):
                return str(s)
            fail_rows = []
            for k in failures_map.keys():
                smi = k[0]
                model = (k[1] if has_model_name else 'single')
                entry = {'smiles': smi, 'model_name': model}
                mol_idx_val = None
                if identity_map:
                    mol_idx_val = identity_map.get((smi, model))
                    if mol_idx_val is None:
                        mol_idx_val = identity_map.get((_canon(smi), model))
                if mol_idx_val is not None:
                    entry['mol_idx'] = mol_idx_val
                fail_rows.append(entry)
            fail_df = pd.DataFrame(fail_rows)
            id_cols = ['smiles', 'model_name', 'mol_idx']
            ordered = [c for c in id_cols if c in fail_df.columns] + [c for c in fail_df.columns if c not in id_cols]
            fail_df = fail_df[ordered]
            fail_df.to_csv(folder_to_save + 'failStructFiltersSMILES.csv', index=False)

            records = []
            for (smi, model), filters in failures_map.items():
                rec = {'smiles': smi,
                       'model_name': (model if has_model_name else 'single'),
                       'failed_filters': ';'.join(sorted(filters))}
                if identity_map:
                    mol_idx_val = identity_map.get((smi, rec['model_name']))
                    if mol_idx_val is None:
                        mol_idx_val = identity_map.get((_canon(smi), rec['model_name']))
                    if mol_idx_val is not None:
                        rec['mol_idx'] = mol_idx_val
                records.append(rec)
            df_fail_summary = pd.DataFrame(records)
            cols = ['smiles', 'model_name'] + (['mol_idx'] if 'mol_idx' in df_fail_summary.columns else []) + ['failed_filters']
            df_fail_summary = df_fail_summary[cols]
            df_fail_summary.to_csv(folder_to_save + 'failStructFiltersMetrics.csv', index=False)

    return filtered_data


def inject_identity_columns_to_all_csvs(config, prefix):
    base_folder = process_path(config['folder_to_save'])
    id_path = base_folder + 'Descriptors/passDescriptorsSMILES.csv'
    if not os.path.exists(id_path):
        return

    try:
        id_df = pd.read_csv(id_path)
        if 'model_name' not in id_df.columns:
            id_df['model_name'] = 'single'
        keep_cols = ['smiles', 'model_name'] + (['mol_idx'] if 'mol_idx' in id_df.columns else [])
        id_df = id_df[keep_cols].copy()
    except Exception:
        return

    if prefix == 'beforeDescriptors':
        target_folder = process_path(config['folder_to_save'], key_word=f'{prefix}_StructFilters')
    else:
        target_folder = process_path(config['folder_to_save'], key_word='StructFilters')

    csv_paths = glob.glob(target_folder + '*.csv')
    for path in csv_paths:
        try:
            df = pd.read_csv(path)
            if 'smiles' not in df.columns:
                continue
            right = id_df.rename(columns={'mol_idx': 'mol_idx_id'}) if 'mol_idx' in id_df.columns else id_df.copy()
            if 'model_name' in df.columns:
                merged = df.merge(right, on=['smiles', 'model_name'], how='left')
            else:
                right_dedup = right.drop_duplicates('smiles')
                merged = df.merge(right_dedup, on='smiles', how='left')

            if 'mol_idx_id' in merged.columns:
                if 'mol_idx' in merged.columns:
                    merged['mol_idx'] = merged['mol_idx'].fillna(merged['mol_idx_id'])
                else:
                    merged['mol_idx'] = merged['mol_idx_id']
                merged.drop(columns=['mol_idx_id'], inplace=True)

            identity_order = ['smiles', 'model_name', 'mol_idx']
            ordered = [c for c in identity_order if c in merged.columns] + [c for c in merged.columns if c not in identity_order]
            merged = merged[ordered]

            merged.to_csv(path, index=False)
        except Exception:
            continue
