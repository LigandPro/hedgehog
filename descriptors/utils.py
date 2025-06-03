import ast
import json
import math 
import random

import numpy as np
import pandas as pd
import datamol as dm
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt

from rdkit import Chem, RDLogger, rdBase

from rdkit.Chem import DataStructs, Lipinski, rdMolDescriptors, QED, RDConfig, Draw
from rdkit.Chem import AllChem as Chem

import os 
import sys
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

from syba import syba
from modules.mce18 import MCE18

from modules.MolScore.moleval.metrics.metrics import GetMetrics

from medchem.rules._utils import n_fused_aromatic_rings

from logger_config import logger
from config_utils import load_config

RDLogger.DisableLog('rdApp.*')
rdBase.DisableLog('rdApp.*')
dm.disable_rdkit_log() 



config = load_config()
n_jobs = config['descriptors']['n_jobs']
device = config['descriptors']['device']
batch_size = config['descriptors']['batch_size']



def load_inhibitors(path_to_save):
    logger.info("Loading data...")
    df_inhibitors_seen = pd.read_csv(config['data']['molecules_to_compare_with_seen'], sep="\t", names=["smiles"], header=None)
    df_inhibitors_unseen = pd.read_csv(config['data']['molecules_to_compare_with_unseen'], sep="\t", names=["smiles"], header=None)

    def _preprocess_inhibitors(i, row):
        mol = dm.to_mol(row["smiles"], ordered=True)
        if not mol: return row
        
        with dm.without_rdkit_log():
            mol = dm.fix_mol(mol)
            if not mol: return row
            
        mol = dm.sanitize_mol(mol, sanifix=True, charge_neutral=False)
        if not mol: return row
        mol = dm.standardize_mol(
            mol, disconnect_metals=False, normalize=True, reionize=True,
            uncharge=False, stereo=True
        )
        if not mol: return row
        row["standard_smiles"] = dm.standardize_smiles(dm.to_smiles(mol))
        return row
    
    inhibitors_data = []
    for df in [df_inhibitors_seen, df_inhibitors_unseen]:
        inhib_data = dm.parallelized(_preprocess_inhibitors,
                                     df.iterrows(),
                                     arg_type="args",
                                     progress=True,
                                     total=len(df)
                                    )
        inhib_data = pd.DataFrame(inhib_data)
        inhib_data.dropna(subset=["standard_smiles"], inplace=True)
        inhib_data.drop_duplicates(subset=["standard_smiles"], inplace=True)
        inhibitors_data.append(inhib_data)

    df_inhibitors_clean_seen, df_inhibitors_clean_unseen = inhibitors_data

    df_inhibitors_total = pd.concat([df_inhibitors_clean_seen, df_inhibitors_clean_unseen])
    df_inhibitors_total.to_csv(f'{path_to_save}inhibs_total.smi', index=False)

    logger.info(f"Number of unique inhibitors seen: {len(df_inhibitors_clean_seen)}")
    logger.info(f"Number of unique inhibitors unseen: {len(df_inhibitors_clean_unseen)}")
    logger.info(f"Number of total inhibitors (seen + unseen): {len(df_inhibitors_total)}")

    return df_inhibitors_clean_seen, df_inhibitors_clean_unseen, df_inhibitors_total


def loading_generated_mols(generated_mols_path):
    logger.info("Loading generated SDF files...")

    def load_and_clean_sdf(sdf_file):
        import logging
        logging.getLogger("rdkit").setLevel(logging.ERROR)
       
        smiles = pd.read_csv(sdf_file, names=["smiles"])['smiles'].tolist()
        smiles = random.sample(smiles, 10000) if len(smiles) > 10000 else smiles
        df_tmp = pd.DataFrame({"smiles": smiles})
        df_tmp.dropna(subset=["smiles"], inplace=True)

        def _preprocess_sdf(i, row):
            mol = dm.to_mol(row["smiles"], ordered=True)
            if mol:
                row["standard_smiles"] = dm.standardize_smiles(dm.to_smiles(mol))
            return row


        df_proc = dm.parallelized(
            _preprocess_sdf,
            df_tmp.iterrows(),
            arg_type="args",
            progress=True,
            total=len(df_tmp)
        )
        df_proc = pd.DataFrame(df_proc)
        df_proc.dropna(subset=["standard_smiles"], inplace=True)
        return df_proc

    dict_generated = {}
    base = Path(generated_mols_path).name
    df_gen = load_and_clean_sdf(generated_mols_path)
    dict_generated[base] = df_gen
    logger.info(f"{base} unique molecules: {len(df_gen)}")

    return dict_generated


def check_intersection(dict_generated, df_inhibitors_clean_seen, df_inhibitors_clean_unseen, df_inhibitors_total):
    for name, df_gen in dict_generated.items():
        for i, df_inhibitors_clean in enumerate([df_inhibitors_clean_seen, df_inhibitors_clean_unseen, df_inhibitors_total]):
            merged = pd.merge(df_inhibitors_clean[["standard_smiles"]],
                              df_gen[["standard_smiles"]],
                              on="standard_smiles",
                              how="inner"
                             )
            if i == 0:
                logger.info(f"Intersection with known SEEN inhibitors. {name} intersection count: {len(merged)}")
            elif i == 1:
                logger.info(f"Intersection with known UNSEEN inhibitors. {name} intersection count: {len(merged)}")
            else:
                logger.info(f"Intersection with all inhibitors. {name} intersection count: {len(merged)}")


def tanimoto_similarity_claculation(dict_generated, df_inhibitors_total, path_to_save):
    logger.info("Computing fingerprints for inhibitors...")
    fps_inhibitors = []
    
    for smi in df_inhibitors_total["standard_smiles"]:
        m = dm.to_mol(smi)
        if m:
            fp = dm.to_fp(m, fp_type="ecfp", as_array=False)
            fps_inhibitors.append(fp)


    def tanimoto_to_inhibitors(df_gen):
        similarities = []
        for smi in df_gen["standard_smiles"]:
            mol = dm.to_mol(smi)
            if mol:
                fp = dm.to_fp(mol, fp_type="ecfp", as_array=False)
                tanimoto_values = [DataStructs.TanimotoSimilarity(fp, fi) for fi in fps_inhibitors]
                max_sim = max(tanimoto_values) if tanimoto_values else 0
                similarities.append(max_sim)
            else:
                similarities.append(0)
        return similarities


    logger.info("Calculating Tanimoto similarity distribution...")
    for name, df_gen in dict_generated.items():
        sims = tanimoto_to_inhibitors(df_gen)
        dict_generated[name]["tanimoto_to_inhibitors"] = sims
        logger.info(f"{name:<20} Tanimoto average: {np.mean(sims):.2f}" if sims else f"{name:<20} Tanimoto average: 0.00")

    plt.figure(figsize=(10, 6))
    for name, df_gen in dict_generated.items():
        plt.hist(df_gen["tanimoto_to_inhibitors"],
                 bins=100,
                 alpha=0.4,
                 label=name
                )
    plt.title("Tanimoto Similarity to Inhibitors")
    plt.xlabel("Tanimoto")
    plt.ylabel("Frequency")
    plt.legend()
    plt.savefig(f"{path_to_save}tanimotoSim.png")

    return


def compute_descriptors(dict_generated):
    logger.info("Computing mean descriptors for molecules...")
    desc_dict = {}
    for name, df_gen in dict_generated.items():
        mols = [mol for mol in (dm.to_mol(smi) 
                                for smi in df_gen["standard_smiles"]) 
                if mol is not None]
        df_desc = dm.descriptors.batch_compute_many_descriptors(mols, progress=True)
        df_desc = df_desc.dropna()
        desc_dict[name] = df_desc

    return desc_dict


def compute_mce18_score(dict_generated):
    if dict_generated:
        max_name_len = max(len(name)  for name in dict_generated.keys())
        results = []
        for name, df_gen in dict_generated.items():
            mols = [dm.to_mol(smi) for smi in df_gen["standard_smiles"]]
            mce_scores = dm.parallelized(_compute_mce18, mols, n_jobs=n_jobs, progress=True)
            mce_scores = [s for s in mce_scores if s is not None]
            mean_mce = np.mean(mce_scores) if mce_scores else 0
            dict_generated[name]["mce18_scores"] = mce_scores
            results.append((name, mean_mce))
        
        for name, mean_mce in results:
            logger.info(f"Mean MCE-18 values: {name:<{max_name_len}} | {mean_mce:.2f}")

    return dict_generated


def _compute_mce18(mol):
    try:
        return MCE18(mol).CalculateMCE18() if mol else 0
    except Exception as e:
        print(f"Error calculating MCE18 for molecule: {e}")
        return 0
        

def compute_syba_score(dict_generated, syba_model):
    if syba_model is None:
        for name, df_gen in dict_generated.items():
            dict_generated[name]["syba_scores"] = [0] * len(df_gen)
        logger.warning("SYBA model not available, setting all scores to 0")
        return dict_generated

    for name, df_gen in dict_generated.items():
        mols = [dm.to_mol(smi) for smi in df_gen["standard_smiles"]]
        smiles_list = [Chem.MolToSmiles(m) if m is not None
                        else None
                        for m in mols]

        syba_scores = dm.parallelized(_compute_syba, 
                                      [(smiles, syba_model) for smiles in smiles_list],
                                      arg_type="args",
                                      n_jobs=n_jobs, 
                                      progress=True,
                                      )
        syba_scores = [v if not np.isnan(v) 
                       else 0
                       for v in syba_scores 
                       ]
        mean_syba = np.mean(syba_scores) if len(syba_scores) > 0 else 0
        dict_generated[name]["syba_scores"] = syba_scores

    max_name_len = max(len(name) for name in dict_generated.keys())
    for name, df in dict_generated.items():
        syba_scores = df["syba_scores"]
        mean_syba = np.mean(syba_scores) if len(syba_scores) > 0 else 0
        logger.info(f"Mean SYBA values: {name:<{max_name_len}} | {mean_syba:.2f}")

    return dict_generated


def _compute_syba(smiles, syba_model):
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return syba_model.predict(mol=mol)
    return np.nan
    

def computing_target_descriptors(dataset_to_compare_with, syba_model):
    logger.info("Computing descriptors for target dataset...")
    target_desc = dm.descriptors.batch_compute_many_descriptors(
        [dm.to_mol(smi) for smi in dataset_to_compare_with["standard_smiles"]]
    ).dropna()

    target_mols = [dm.to_mol(smi) for smi in dataset_to_compare_with["standard_smiles"]]
    target_mce_values = dm.parallelized(_compute_mce18, target_mols, n_jobs=n_jobs, progress=True)
    target_mce_scores = [s for s in target_mce_values if s is not None]
    mean_target_mce = np.mean(target_mce_scores) if target_mce_scores else 0

    if syba_model is not None:
        target_smiles = [Chem.MolToSmiles(m) for m in target_mols if m is not None]
        target_syba_values = dm.parallelized(_compute_syba, 
                                             [(smiles, syba_model) for smiles in target_smiles],
                                             arg_type="args",
                                             n_jobs=n_jobs, 
                                             progress=True)
        target_syba_scores = [v for v in target_syba_values if not np.isnan(v)]
        mean_target_syba = np.mean(target_syba_scores) if target_syba_scores else 0
        logger.info(f"Target dataset: Mean SYBA: {mean_target_syba:.2f}")
    else:
        target_syba_scores = [0] * len(target_mols)
        logger.warning("SYBA model not available, setting all target scores to 0")

    logger.info(f"Target dataset: Mean MCE-18: {mean_target_mce:.2f}")

    df_target = pd.DataFrame({"standard_smiles": dataset_to_compare_with["standard_smiles"],
                              "mce18_scores": target_mce_scores,
                              "syba_scores": target_syba_scores
                            })

    return df_target



def collect_metrics_dict(dict_generated, desc_dict, target_desc):
    metrics_dict = {}
    for name, df_gen in dict_generated.items():
        metrics_dict[name] = {}

        if name in desc_dict:
            numeric_desc = desc_dict[name].select_dtypes(include=[float, int]).to_dict('list')
            metrics_dict[name].update(numeric_desc)

        if isinstance(df_gen, pd.DataFrame):
            for col in df_gen.columns:
                if pd.api.types.is_numeric_dtype(df_gen[col]) or pd.api.types.is_bool_dtype(df_gen[col]):
                    series_unique = df_gen[col].dropna().unique()
                    if len(series_unique) == 1:
                        val = series_unique[0]
                        if col not in metrics_dict[name]:
                            metrics_dict[name][col] = [val]
                    else:
                        if col not in metrics_dict[name]:
                            metrics_dict[name][col] = df_gen[col].dropna().tolist()

        for key, value in dict_generated[name].items():
            if key not in metrics_dict[name] and isinstance(value, (int, float)):
                metrics_dict[name][key] = [value]

    target_key = "Target"
    metrics_dict[target_key] = {}

    if not target_desc.empty:
        numeric_desc_target = target_desc.select_dtypes(include=[float, int]).to_dict('list')
        metrics_dict[target_key].update(numeric_desc_target)

    if target_key in dict_generated:
        df_target = dict_generated[target_key]
        if isinstance(df_target, pd.DataFrame):
            for col in df_target.columns:
                if pd.api.types.is_numeric_dtype(df_target[col]) or pd.api.types.is_bool_dtype(df_target[col]):
                    series_unique = df_target[col].dropna().unique()
                    if len(series_unique) == 1:
                        val = series_unique[0]
                        if col not in metrics_dict[target_key]:
                            metrics_dict[target_key][col] = [val]
                    else:
                        if col not in metrics_dict[target_key]:
                            metrics_dict[target_key][col] = df_target[col].dropna().tolist()
        for k, v in dict_generated[target_key].items():
            if k not in metrics_dict[target_key] and isinstance(v, (int, float)):
                metrics_dict[target_key][k] = [v]

    return metrics_dict


def save_plots(metrics_dict, dict_generated, path_to_save, df_to_compare_with):
    fig_name = f'{path_to_save}molevalDescriptors.png'

    all_metrics = set()
    for name_metrics in metrics_dict.values():
        all_metrics.update(name_metrics.keys())

    ncols = 3
    nrows = math.ceil(len(all_metrics) / ncols)
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(25, 50))
    axs = axs.flatten()

    handles, labels = [], []

    for i, metric in enumerate(all_metrics):
        for name, data_dict in metrics_dict.items():
            if metric in data_dict:
                if name == "Target":
                    hist = axs[i].hist(data_dict[metric], bins=30, alpha=1, label=name, 
                                    edgecolor='black', linewidth=2, facecolor='gray')
                else:
                    hist = axs[i].hist(data_dict[metric], bins=30, alpha=0.4, label=name)
                if i == 0:
                    handles.append(hist[2][0])  
                    labels.append(name)
        axs[i].set_title(metric, fontsize=40)

    for j in range(i + 1, len(axs)):
        fig.delaxes(axs[j])

    plt.tight_layout()
    plt.subplots_adjust(right=0.85)
    fig.legend(handles, labels, loc='center right', bbox_to_anchor=(0.98, 0.5), fontsize=12)
    plt.savefig(fig_name)


    all_data = []
    for dataset_name, metrics_map in metrics_dict.items():
        row = {"dataset": dataset_name}
        for metric_name, values in metrics_map.items():

            if isinstance(values, list) and len(values) > 0:
                row[metric_name + "_mean"] = np.mean(values)
            else:
                row[metric_name + "_mean"] = values
        all_data.append(row)

    df_means = pd.DataFrame(all_data)
    df_means.to_csv(f"{path_to_save}allChptsDescriptors.csv", index=False)

    TARGET_SMILES = df_to_compare_with["standard_smiles"].tolist()
    MetricEngine = GetMetrics(
        n_jobs=n_jobs,
        device=device,
        batch_size=batch_size,
        target=TARGET_SMILES
    )

    moleval_metrics = {}
    for name, df_gen in dict_generated.items():
        gen_smiles = df_gen["standard_smiles"].tolist()  
        metrics_list = MetricEngine.calculate(
            gen_smiles,
            calc_valid=True,
            calc_unique=True,
            se_k=1,
            sp_k=1,
            properties=True,
            return_stats=True
        )
        moleval_metrics[name] = {
            item["metric"]: item["value"] 
            for item in metrics_list 
            if item.get("metric") is not None
        }

    rows = []
    for dataset_name, metrics in moleval_metrics.items():
        rows.append({'Dataset': dataset_name, **metrics})
        

    df_metrics = pd.DataFrame(rows).set_index('Dataset')
    df_metrics.to_csv(f"{path_to_save}molevalMetrics.csv", index=False)
    df_means = pd.read_csv(f"{path_to_save}allChptsDescriptors.csv")

    df_metrics_reset = df_metrics.reset_index().rename(columns={'Dataset': 'dataset'})
    df_merged = pd.merge(df_means, df_metrics_reset, on='dataset', how='outer')

    target_row = df_merged[df_merged['dataset'] == 'Target']
    other_rows = df_merged[df_merged['dataset'] != 'Target']
    df_merged = pd.concat([target_row, other_rows]).reset_index(drop=True)
    df_merged.to_csv(f"{path_to_save}allMetrics.csv", index=False)

    n_datasets = len(df_merged['dataset'].unique())
    colors = ['#FF9999', '#66B2FF', '#99FF99', '#FFCC99', '#FF99CC'][:n_datasets]
    numeric_cols = df_merged.select_dtypes(include=[np.number]).columns
    ncols = 6
    nrows = math.ceil(len(numeric_cols) / ncols)
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(25, 120))
    axs = axs.flatten()

    for i, col in enumerate(numeric_cols):
        sns.barplot(data=df_merged, x='dataset', y=col, hue='dataset', ax=axs[i], palette=colors, legend=False)
        axs[i].set_title(col, fontsize=30, pad=40)
        axs[i].tick_params(axis='x', rotation=90, labelsize=14)
        axs[i].tick_params(axis='y', labelsize=12)
        
        
    for j in range(i + 1, len(axs)):
        fig.delaxes(axs[j])
        
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig(f"{path_to_save}molevalMetrics.png", bbox_inches='tight')

    import json
    with open(f"{path_to_save}metricsDict.json", "w") as f:
        json.dump(metrics_dict, f)
    with open(f"{path_to_save}molevalMetrics.json", "w") as f:
        json.dump(moleval_metrics, f)


def canonicalize_smarts(smarts):
    if smarts == 'clean molecule': return smarts
    mol = Chem.MolFromSmarts(smarts)
    if mol: return Chem.MolToSmarts(mol)
    else: raise ValueError(f"Invalid SMARTS: {smarts}")
   

def clean_path_name(path_name, patterns_to_remove=None):
    if patterns_to_remove is None:
        patterns_to_remove = ['.csv', '.sdf', '.txt', 
                              '_smiles_metrics.csv', '_smiles_metrics_seen.csv', '_smiles_metrics_filtered.csv',
                              '_10000_mols', '_10000', '_sampling_seen', 'sampling_10000_', 
                              'chkpt', '_sampling', '_smiles_additional_metrics.csv', '_smiles', '_merged', 
                              ]
    
    name = os.path.basename(path_name)
    for pattern in patterns_to_remove:
        if pattern == 'chkpt': name = name.replace(pattern, 'ckpt')
        elif pattern == 'reinvent4_epochs': name = name.replace(pattern, 'epochs')
        else: name = name.replace(pattern, '')

    return name


def get_model_colors(model_names, cmap=None):
    return dict(zip(model_names, plt.cm.YlOrRd(np.linspace(1, 0, len(model_names) + 1)) if cmap is None else plt.cm.get_cmap(cmap)(np.linspace(1, 0, len(model_names) + 1))))


def get_smarts_to_regid_mapping(pains_file_path):
    pains_df = pd.read_csv(pains_file_path, names=['SMARTS', 'regID'])
    pains_col = list(canonicalize_smarts(row.strip()) for row in pains_df['SMARTS'])
    regid_col = list(row.strip('"<>') for row in pains_df['regID'])
    smarts_to_regid = dict(zip(pains_col, regid_col))
    return pains_col, smarts_to_regid


def get_structure_filters_data(is_pains):
    if is_pains:
        pains_file_path = config['data']['pains_file_path']
        pains_list, smarts_to_regid = get_smarts_to_regid_mapping(pains_file_path)
        mcf_list = None 
    else:
        smarts_to_regid = {}
        pains_list = None
        mcf_file_path = config['data']['mcf_file_path']
        mcf_list_not_canonized_smarts = pd.read_csv(mcf_file_path, names=['names', 'smarts'])['smarts'][1:].tolist()
        mcf_list = list(canonicalize_smarts(row.strip()) for row in mcf_list_not_canonized_smarts)
    return smarts_to_regid, pains_list, mcf_list


def get_counts_from_csv(df, target_column, passed_column_name, value_column_name, plot_true_vals):
    if target_column == passed_column_name:
        if plot_true_vals:
            true_count = df[passed_column_name].value_counts()
            if True in true_count: return true_count[True]
            else: return 0

        else:
            false_count = df[passed_column_name].value_counts()
            if False in false_count: return false_count[False]
            else: return 0
    
    elif target_column == value_column_name:
        values_raw = df[value_column_name].dropna()
        values = []
        for value in values_raw:
            if type(value) in (int, float, np.int64, np.float64): values.append(value)
            elif type(ast.literal_eval(value)) in (int, float): values.append(value)
            elif isinstance(value, str):
                if type(ast.literal_eval(value)) == list:
                    values.extend(ast.literal_eval(value))

        if type(values) == list:
            values = pd.Series(list(map(float, values)))
            return values


def sort_data(data_dict, reverse=True):
    sorted_items = sorted(data_dict.items(), key=lambda x: x[1], reverse=reverse)
    model_names = [item[0] for item in sorted_items]
    counts = [item[1] for item in sorted_items]
    return model_names, counts


def refilter_data(paths, borders=None, save=True, path_to_save=None):
    if borders is None:
        borders = {'molWt_max': 1000, 
           'logP_min': -2, 
           'logP_max': 10, 
           'hbd': 6, 
           'hba': 15, 
           'tpsa': 250, 
           'num_rot_bonds_max': 20}
    logger.info(f'Borders: {borders}')
    datas = []
    for path in paths: 
        if path_to_save is None: 
            out_path = path.replace('.csv', '_filtered.csv')
        else: 
            name = path.split('/')[-1].split('.csv')[0]
            out_path = path_to_save + name + '_filtered.csv'
    
        data = pd.read_csv(path, index_col=0)
        
        if 'filters_logP_value' in data.columns:
            data['filters_logP_passed'] = ((data['filters_logP_value'] >= borders['logP_min']) & (data['filters_logP_value'] <= borders['logP_max']))
        if 'filters_molWt_value' in data.columns:
            data['filters_molWt_passed'] = data['filters_molWt_value'] <= borders['molWt_max']
        if 'filters_hbd_value' in data.columns: 
            data['filters_hbd_passed'] = data['filters_hbd_value'] <= borders['hbd']
        if 'filters_hba_value' in data.columns:
            data['filters_hba_passed'] = data['filters_hba_value'] <= borders['hba']
        if 'filters_tpsa_value' in data.columns:
            data['filters_tpsa_passed'] = data['filters_tpsa_value'] <= borders['tpsa']
        if 'filters_num_rot_bonds_value' in data.columns:
            data['filters_num_rot_bonds_passed'] = data['filters_num_rot_bonds_value'] <= borders['num_rot_bonds_max']
        
        if save: 
            data.to_csv(out_path)
            print(f'{path} saved')
        
        datas.append(data)
    
    assert len(datas) == len(paths), 'Number of filtered data is not equal to number of paths.'

    if save: return
    else: return datas


def draw_smarts(smarts):
    mol = Chem.MolFromSmarts(smarts)

    if mol is not None:
        img = Draw.MolToImage(mol)
        img.show()
    else:
        logger.error("Invalid SMARTS pattern!")
    return 


def sdf_to_txt(sdf_file_path, txt_file_path):
    supplier = Chem.SDMolSupplier(sdf_file_path)
    with open(txt_file_path, "w") as out:
        for mol in supplier:
            if mol is not None:
                smiles = Chem.MolToSmiles(mol) 
                out.write(smiles + "\n")
    logger.info(f'{txt_file_path} saved')
    return


def calculate_metrics(df, syba_model, save_path=None):
    if type(df) == dict:
        df = pd.DataFrame(df)
    df_name = df.columns[0]
    metrics = {}
    skipped_molecules = []
    for smiles in df[df.columns[0]]:
        mol = Chem.MolFromSmiles(smiles)
        mol_metrics = {}
        if mol:
            mol_metrics['hbd'] = Lipinski.NumHDonors(mol)
            mol_metrics['hba'] = Lipinski.NumHAcceptors(mol)
            mol_metrics['tpsa'] = rdMolDescriptors.CalcTPSA(mol)
            mol_metrics['qed'] = QED.qed(mol)
            mol_metrics['sp3'] = rdMolDescriptors.CalcFractionCSP3(mol)
            mol_metrics['num_heavy_atoms'] = rdMolDescriptors.CalcNumHeavyAtoms(mol)
            mol_metrics['sas_score'] = _calculate_sas(smiles)
            mol_metrics['syba_score'] = _calculate_syba(smiles, syba_model)
            metrics[smiles] = mol_metrics
        else:
            skipped_molecules.append(smiles)
        
    if save_path:
        if save_path.endswith('/'):
            path_to_save = save_path + f'{df_name}additionalMetrics.csv'
        else:
            path_to_save = save_path + f'/{df_name}additionalMetrics.csv'
        metrics_df = pd.DataFrame.from_dict(metrics, orient='index')
        metrics_df.to_csv(path_to_save, index_label=df_name)

    if skipped_molecules:
        logger.warning(f'Skipped {len(skipped_molecules)} molecules: {skipped_molecules}')

    return metrics


def calculate_cycle_metrics(df, save_path=None):
    assert type(df) == pd.DataFrame or type(df) == pd.Series, 'df must be a pandas DataFrame or Series'

    df_name = df.columns[0]
    metrics = {}
    skipped_molecules = []
    for smiles in df[df.columns[0]]:
        mol = Chem.MolFromSmiles(smiles)
        mol_metrics = {}
        if mol:
            mol_metrics['#rings'] = mol.GetRingInfo().NumRings()
            mol_metrics['#aromatic_rings'] = rdMolDescriptors.CalcNumAromaticRings(mol)
            mol_metrics['#heteroatoms'] = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in (1, 6))
            mol_metrics['#cycles_together'] = n_fused_aromatic_rings(mol)

            total_bonds = mol.GetNumBonds()
            rotatable_bonds = Lipinski.NumRotatableBonds(mol)
            rigid_bonds = total_bonds - rotatable_bonds
            mol_metrics['#rigid_bonds'] = rigid_bonds
            metrics[smiles] = mol_metrics
        else:
            skipped_molecules.append(smiles)
    if save_path:
        if save_path.endswith('/'):
            path_to_save = save_path + f'{df_name}cycleMetrics.csv'
        else:
            path_to_save = save_path + f'/{df_name}cycleMetrics.csv'
        metrics_df = pd.DataFrame.from_dict(metrics, orient='index')
        metrics_df.to_csv(path_to_save, index_label=df_name)

    if skipped_molecules:
        logger.warning(f'Skipped {len(skipped_molecules)} molecules: {skipped_molecules}')
    return     

          
def _calculate_sas(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return sascorer.calculateScore(mol)


def _calculate_syba(smiles, syba_model):
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return syba_model.predict(mol=mol)
    return np.nan

