import ast
import json
import glob

import math 
import random

import numpy as np
import pandas as pd
import datamol as dm
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt

from rdkit import Chem, RDLogger, rdBase

from rdkit.Chem import DataStructs, Lipinski, rdMolDescriptors, QED, RDConfig, Draw, Descriptors
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
from configs.config_utils import load_config

RDLogger.DisableLog('rdApp.*')
rdBase.DisableLog('rdApp.*')
dm.disable_rdkit_log() 


def read_syba():
    syba_model = syba.SybaClassifier()
    try:
        syba_model.fitDefaultScore()
    except Exception as e:
        logger.error(f"Failed to load SYBA model: {str(e)}")
        syba_model = None
    return syba_model


def load_inhibitors(target_molecules_path):
    logger.info("Loading target molecules data...")

    df_target_mols = pd.read_csv(target_molecules_path, sep="\t", names=["smiles"], header=None)


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
    
    target_data = dm.parallelized(_preprocess_inhibitors,
                                  df_target_mols.iterrows(),
                                  arg_type="args",
                                  progress=True,
                                  total=len(df_target_mols)
                                 )
    
    target_data = pd.DataFrame(target_data)
    target_data.dropna(subset=["standard_smiles"], inplace=True)
    target_data.drop_duplicates(subset=["standard_smiles"], inplace=True)

    logger.info(f"Number of unique target molecules: {len(target_data)}")
    return target_data


def loading_generated_mols(data, model_name):
    logger.info("Loading generated molecules...")

    def load_and_clean(data):
        import logging
        logging.getLogger("rdkit").setLevel(logging.ERROR)
       
        smiles = data["smiles"].tolist()
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
    df_gen = load_and_clean(data)
    dict_generated[model_name] = df_gen
    logger.info(f"{model_name} unique molecules: {len(df_gen)}")

    return dict_generated


def check_intersection(dict_generated, target_data):
    for name, df_gen in dict_generated.items():
        
        merged = pd.merge(target_data[["standard_smiles"]],
                            df_gen[["standard_smiles"]],
                            on="standard_smiles",
                            how="inner"
                            )
        
        logger.info(f"Intersection with target molecules: {len(merged)}")


def tanimoto_similarity_claculation(dict_generated, target_data, folder_to_save):
    logger.info("Computing fingerprints for inhibitors...")
    fps_inhibitors = []
    
    for smi in target_data["standard_smiles"]:
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
    plt.savefig(f"{folder_to_save}tanimotoSim.png")

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


def compute_mce18_score(dict_generated, n_jobs):
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
        

def compute_syba_score(dict_generated, syba_model, n_jobs):
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
    

def computing_target_descriptors(dataset_to_compare_with, syba_model, n_jobs):

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


def save_plots(metrics_dict, dict_generated, folder_to_save, target_data, n_jobs, batch_size):
    fig_name = f'{folder_to_save}meanDescriptors.png'

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
    df_means.to_csv(f"{folder_to_save}meanMetrics.csv", index=False)

    TARGET_SMILES = target_data["standard_smiles"].tolist()
    MetricEngine = GetMetrics(
        n_jobs=n_jobs,
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
    df_metrics.to_csv(f"{folder_to_save}relativeToTargetMetrics.csv", index=False)
    df_means = pd.read_csv(f"{folder_to_save}meanMetrics.csv")

    df_metrics_reset = df_metrics.reset_index().rename(columns={'Dataset': 'dataset'})
    df_merged = pd.merge(df_means, df_metrics_reset, on='dataset', how='outer')

    target_row = df_merged[df_merged['dataset'] == 'Target']
    other_rows = df_merged[df_merged['dataset'] != 'Target']
    df_merged = pd.concat([target_row, other_rows]).reset_index(drop=True)
    df_merged.to_csv(f"{folder_to_save}allMetrics.csv", index=False)

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
    plt.savefig(f"{folder_to_save}relativeToTargetMetrics.png", bbox_inches='tight')

    with open(f"{folder_to_save}metricsDict.json", "w") as f:
        json.dump(metrics_dict, f)
    with open(f"{folder_to_save}relativeToTargetMetrics.json", "w") as f:
        json.dump(moleval_metrics, f)


def canonicalize_smarts(smarts):
    if smarts == 'clean molecule': return smarts
    mol = Chem.MolFromSmarts(smarts)
    if mol: return Chem.MolToSmarts(mol)
    else: raise ValueError(f"Invalid SMARTS: {smarts}")


def get_model_colors(model_names, cmap=None):
    return dict(zip(model_names, plt.cm.YlOrRd(np.linspace(1, 0, len(model_names) + 1)) if cmap is None else plt.cm.get_cmap(cmap)(np.linspace(1, 0, len(model_names) + 1))))


def get_smarts_to_regid_mapping(pains_file_path):
    pains_df = pd.read_csv(pains_file_path, names=['SMARTS', 'regID'])
    pains_col = list(canonicalize_smarts(row.strip()) for row in pains_df['SMARTS'])
    regid_col = list(row.strip('"<>') for row in pains_df['regID'])
    smarts_to_regid = dict(zip(pains_col, regid_col))
    return pains_col, smarts_to_regid


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


def compute_metrics(df, syba_model, model_name, save_path=None):
    if type(df) == dict:
        df = pd.DataFrame(df)
    metrics = {}
    skipped_molecules = []
    for smiles in df[df.columns[0]]:
        mol = Chem.MolFromSmiles(smiles)
        mol_metrics = {}
        if mol:
            symbols = list(set(atom.GetSymbol() 
                         for atom in mol.GetAtoms() 
                         if atom.GetSymbol()))
            
            
            charged_mol = False if any(atom.GetFormalCharge() != 0 for atom in mol.GetAtoms()) else True

            ring_info = mol.GetRingInfo()
            rings = [len(x) for x in ring_info.AtomRings()]

            total_bonds = mol.GetNumBonds()
            rotatable_bonds = Lipinski.NumRotatableBonds(mol)
            rigid_bonds = total_bonds - rotatable_bonds

            mol_metrics['chars'] = symbols
            mol_metrics['n_atoms'] = Chem.AddHs(mol).GetNumAtoms()
            mol_metrics['n_heavy_atoms'] = rdMolDescriptors.CalcNumHeavyAtoms(mol)
            mol_metrics['charged_mol'] = charged_mol
            mol_metrics['molWt'] = Descriptors.ExactMolWt(mol)
            mol_metrics['logP'] = Descriptors.MolLogP(mol)
            mol_metrics['ring_size'] = rings
            mol_metrics['n_rings'] = mol.GetRingInfo().NumRings()
            mol_metrics['n_aroma_rings'] = rdMolDescriptors.CalcNumAromaticRings(mol)
            mol_metrics['n_rings_together'] = n_fused_aromatic_rings(mol)
            mol_metrics['n_het_atoms'] = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in (1, 6))
            mol_metrics['n_rigid_bonds'] = rigid_bonds
            mol_metrics['n_rot_bonds'] = rdMolDescriptors.CalcNumRotatableBonds(mol)
            mol_metrics['hbd'] = Lipinski.NumHDonors(mol)
            mol_metrics['hba'] = Lipinski.NumHAcceptors(mol)
            mol_metrics['fsp3'] = rdMolDescriptors.CalcFractionCSP3(mol)
            mol_metrics['tpsa'] = rdMolDescriptors.CalcTPSA(mol)
            mol_metrics['qed'] = QED.qed(mol)
            mol_metrics['sa_score'] = _calculate_sas(smiles)
            mol_metrics['syba_score'] = _calculate_syba(smiles, syba_model)
            metrics[smiles] = mol_metrics
        else:
            skipped_molecules.append(smiles)
        
    if save_path:
        if save_path.endswith('/'):
            folder_to_save = save_path + f'perMolMetrics.csv'
        else:
            folder_to_save = save_path + f'/perMolMetrics.csv'
        metrics_df = pd.DataFrame.from_dict(metrics, orient='index').reset_index().rename(columns={'index': 'smiles'})
        metrics_df.to_csv(folder_to_save, index_label=model_name)

    if skipped_molecules:
        logger.warning(f'Skipped {len(skipped_molecules)} molecules: {skipped_molecules}')

    return metrics_df  


def filter_molecules(df, borders, folder_to_save):
    logger.info(f'Borders: {borders}')

    allowed_chars = borders['allowed_chars']
    n_atoms_min = borders['n_atoms_min']
    n_atoms_max = borders['n_atoms_max']
    n_heavy_atoms_min = borders['n_heavy_atoms_min']
    n_heavy_atoms_max = borders['n_heavy_atoms_max']
    charged_mol_allowed = borders['charged_mol_allowed']
    molWt_min = borders['molWt_min']
    molWt_max = borders['molWt_max']
    logP_min = borders['logP_min']
    logP_max = borders['logP_max']
    ring_size_min = borders['ring_size_min']
    ring_size_max = borders['ring_size_max']
    n_rings_min = borders['n_rings_min']
    n_rings_max = borders['n_rings_max']
    n_aroma_rings_min = borders['n_aroma_rings_min']
    n_aroma_rings_max = borders['n_aroma_rings_max']
    n_rings_together_min = borders['n_rings_together_min']
    n_rings_together_max = borders['n_rings_together_max']  
    n_het_atoms_min = borders['n_het_atoms_min']
    n_het_atoms_max = borders['n_het_atoms_max']  
    n_rigid_bonds_min = borders['n_rigid_bonds_min']
    n_rigid_bonds_max = borders['n_rigid_bonds_max']
    n_rot_bonds_min = borders['n_rot_bonds_min']
    n_rot_bonds_max = borders['n_rot_bonds_max']
    hbd_min = borders['hbd_min']
    hbd_max = borders['hbd_max']
    hba_min = borders['hba_min']
    hba_max = borders['hba_max']
    fsp3_min = borders['fsp3_min']
    fsp3_max = borders['fsp3_max']
    tpsa_min = borders['tpsa_min']
    tpsa_max = borders['tpsa_max']
    qed_min = borders['qed_min']
    qed_max = borders['qed_max']
    sa_score_min = borders['sa_score_min']
    sa_score_max = borders['sa_score_max']
    syba_score_min = borders['syba_score_min']
    syba_score_max = borders['syba_score_max']


    filtered_data = {}
    filtered_data['smiles'] = df['smiles']
    
    filtered_data['chars'] = df['chars']
    filtered_data['chars_passed'] = df['chars'].apply(lambda x: 
                                                             all(str(char).strip() in allowed_chars 
                                                                 for char in (x if isinstance(x, list) else eval(x))))

    filtered_data['n_atoms'] = df['n_atoms']
    filtered_data['n_atoms_passed'] = (df['n_atoms'] >= n_atoms_min) & (df['n_atoms'] <= n_atoms_max)
    
    filtered_data['n_heavy_atoms'] = df['n_heavy_atoms']
    filtered_data['n_heavy_atoms_passed'] = (df['n_heavy_atoms'] >= n_heavy_atoms_min) & (df['n_heavy_atoms'] <= n_heavy_atoms_max)
    
    filtered_data['charged_mol'] = df['charged_mol']
    filtered_data['charged_mol_passed'] = df['charged_mol'] == charged_mol_allowed
    
    filtered_data['molWt'] = df['molWt']
    filtered_data['molWt_passed'] = (df['molWt'] >= molWt_min) & (df['molWt'] <= molWt_max)
    
    filtered_data['logP'] = df['logP']
    filtered_data['logP_passed'] = (df['logP'] >= logP_min) & (df['logP'] <= logP_max)
    
    filtered_data['ring_size'] = df['ring_size']
    filtered_data['ring_size_passed'] = df['ring_size'].apply(lambda x: 
                                                                all(float(ring_size) >= ring_size_min and float(ring_size) <= ring_size_max 
                                                                   for ring_size in (x if isinstance(x, list) else eval(x))))
    
    filtered_data['n_rings'] = df['n_rings']
    filtered_data['n_rings_passed'] = (df['n_rings'] >= n_rings_min) & (df['n_rings'] <= n_rings_max)

    filtered_data['n_aroma_rings'] = df['n_aroma_rings']
    filtered_data['n_aroma_rings_passed'] = (df['n_aroma_rings'] >= n_aroma_rings_min) & (df['n_aroma_rings'] <= n_aroma_rings_max)

    filtered_data['n_rings_together'] = df['n_rings_together']
    filtered_data['n_rings_together_passed'] = (df['n_rings_together'] >= n_rings_together_min) & (df['n_rings_together'] <= n_rings_together_max)

    filtered_data['n_het_atoms'] = df['n_het_atoms']
    filtered_data['n_het_atoms_passed'] = (df['n_het_atoms'] >= n_het_atoms_min) & (df['n_het_atoms'] <= n_het_atoms_max)

    filtered_data['n_rigid_bonds'] = df['n_rigid_bonds']
    filtered_data['n_rigid_bonds_passed'] = (df['n_rigid_bonds'] >= n_rigid_bonds_min) & (df['n_rigid_bonds'] <= n_rigid_bonds_max)
    
    filtered_data['n_rot_bonds'] = df['n_rot_bonds']
    filtered_data['n_rot_bonds_passed'] = (df['n_rot_bonds'] >= n_rot_bonds_min) & (df['n_rot_bonds'] <= n_rot_bonds_max)

    filtered_data['hbd'] = df['hbd']
    filtered_data['hbd_passed'] = (df['hbd'] <= hbd_max) & (df['hbd'] >= hbd_min)

    filtered_data['hba'] = df['hba']
    filtered_data['hba_passed'] = (df['hba'] <= hba_max) & (df['hba'] >= hba_min)

    filtered_data['fsp3'] = df['fsp3']
    filtered_data['fsp3_passed'] = (df['fsp3'] <= fsp3_max) & (df['fsp3'] >= fsp3_min)

    
    filtered_data['tpsa'] = df['tpsa']
    filtered_data['tpsa_passed'] = (df['tpsa'] <= tpsa_max) & (df['tpsa'] >= tpsa_min)

    filtered_data['qed'] = df['qed']
    filtered_data['qed_passed'] = (df['qed'] <= qed_max) & (df['qed'] >= qed_min)

    filtered_data['sa_score'] = df['sa_score']
    filtered_data['sa_score_passed'] = (df['sa_score'] <= sa_score_max) & (df['sa_score'] >= sa_score_min)

    filtered_data['syba_score'] = df['syba_score']
    if syba_score_max == 'inf':
        filtered_data['syba_score_passed'] = (df['syba_score'] >= syba_score_min)
    else:
        filtered_data['syba_score_passed'] = (df['syba_score'] <= syba_score_max) & (df['syba_score'] >= syba_score_min)
    
    filtered_data = pd.DataFrame(filtered_data)
    filtered_data.to_csv(folder_to_save + 'filtered_mols.csv', index_label='SMILES')
    return filtered_data


def _parse_chars_in_mol_column(series):
    parsed = []
    for val in series.dropna():
        try:
            if isinstance(val, list):
                chars = val
            else:
                chars = eval(val)
            parsed.extend(chars)
        except Exception as e:
            print(f"Error processing value: {val}, error: {e}")
    return parsed

def _parse_ring_size_column(series):
    parsed = []
    for val in series.dropna():
        try:
            if isinstance(val, list):
                sizes = val
            else:
                sizes = eval(val)
            parsed.extend([float(size) for size in sizes])
        except Exception as e:
            print(f"Error processing value: {val}, error: {e}")
    return parsed


def draw_filtered_mols(df, folder_to_save, config): 
    model_name = config['model_name']
    
    descriptors_config = load_config(config['config_descriptors'])
    borders = descriptors_config['borders']
    borders['charged_mol_allowed'] = int(borders['charged_mol_allowed'])
    cols_to_plot = descriptors_config['filtered_cols_to_plot']
    discrete_feats = descriptors_config['discrete_features_to_plot']
    not_to_smooth_by_sides_cols = descriptors_config['not_to_smooth_plot_by_sides']
    renamer = descriptors_config['renamer']

    colors = get_model_colors([model_name])

    nrows = len(cols_to_plot) // 5  + len(cols_to_plot) % 5 if len(cols_to_plot) > 5 else 1
    ncols = 5 if len(cols_to_plot) > 5 else len(cols_to_plot)
    fig, axes = plt.subplots(nrows=nrows, ncols=5, figsize=(nrows * ncols, nrows * ncols))
    axes = axes.flatten()

    for i, col in enumerate(cols_to_plot):
        relevant_keys = [k for k in borders.keys() if col.lower() in k.lower()]
        if not relevant_keys:
            continue
        
        min_val = next((borders[k] for k in relevant_keys if 'min' in k), None)
        max_val = next((borders[k] for k in relevant_keys if 'max' in k), None)
        minmax_str = f"min: {min_val}, max: {max_val}"

        data_to_plot_after = {}
        counts_after, total_mols = {}, {}

        if col == 'chars':
            values_raw = df[col].dropna()
            values_before = _parse_chars_in_mol_column(values_raw)
        elif col == 'ring_size':
            values_raw = df[col].dropna()
            values_before = _parse_ring_size_column(values_raw)
        else:
            values_before = df[col].dropna().tolist()

        if col == 'syba_score' and max_val == 'inf':
            values_after = [v for v in values_before if 
                            (min_val is None or v >= min_val)]
        else:
            values_after = [v for v in values_before if 
                            (min_val is None or v >= min_val) and 
                            (max_val is None or v <= max_val)]

        data_to_plot_after[model_name] = values_before
        counts_after[model_name] = len(values_after)
        total_mols[model_name] = len(values_before)

        ax = axes[i]
        for model_name, values in data_to_plot_after.items():
            total = total_mols.get(model_name, 0)
            count = counts_after.get(model_name, 0)
            mols_passed = count / total * 100 if total > 0 else 0

            label = f'{model_name}, passed: {mols_passed:.1f}%'
            color = colors[model_name]

            if len(values) > 1:
                if col in discrete_feats:
                    if col == 'charged_mol':
                        value_counts = pd.Series(values).value_counts().sort_index()
                        complete_counts = pd.Series([0, 0], index=[False, True])
                        complete_counts.update(value_counts)
                        value_counts = complete_counts.sort_index()

                        ax.bar(x=['Not charged', 'Charged'], height=value_counts.values, 
                            alpha=0.5, color=color, edgecolor='black', linewidth=0.3, 
                            label=label)
                        total = len(values)
                        passed = sum(v == borders['charged_mol_allowed'] for v in values)
                        mols_passed = (passed / total * 100) if total > 0 else 0
                        ax.legend([label], loc='upper right', fontsize=8)
                    else:
                        value_counts = pd.Series(values).value_counts().sort_index()
                        ax.bar(value_counts.index, value_counts.values, alpha=0.5, color=color, edgecolor='black', linewidth=0.3, align='edge', width=0.8, label=label)
                elif col in not_to_smooth_by_sides_cols:
                    sns.kdeplot(values, label=label, fill=True, alpha=0.3, ax=ax, color=color, clip=(0, None))
                else:
                    sns.kdeplot(values, label=label, fill=True, alpha=0.3, ax=ax, color=color)
            else:
                ax.scatter(values, [0.01]*len(values), label=label, alpha=0.7, color=color)

        title = f"{renamer[col] if col in renamer else col} ({minmax_str})"
        ax.set_title(title, fontsize=12)
        ax.set_xlabel(renamer[col], fontsize=10)
        
        if col in discrete_feats:
            ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
        ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))


        if min_val is not None:
            ax.axvline(min_val, color='red', linestyle='--', linewidth=1.5, label=f'min: {min_val}')
        if max_val is not None and max_val != 'inf':
            ax.axvline(max_val, color='blue', linestyle='--', linewidth=1.5, label=f'max: {max_val}')
        
        x_min, x_max = ax.get_xlim()
        if min_val is not None:
            ax.axvspan(x_min, min_val, color='grey', alpha=0.2, zorder=0)
        if max_val is not None and max_val != 'inf':
            ax.axvspan(max_val, x_max, color='grey', alpha=0.2, zorder=0)

        ax.tick_params(axis='both', labelsize=10)  
        ax.legend(fontsize=8, loc='upper right')

    for j in range(len(cols_to_plot), len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.subplots_adjust(top=0.93, bottom=0.05, left=0.05, right=0.98)

    plt.savefig(f"{folder_to_save}filteredMetricsDistribution.png", dpi=300, bbox_inches='tight', format='png')
