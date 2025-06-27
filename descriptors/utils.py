import glob

import numpy as np
import pandas as pd
import datamol as dm
import seaborn as sns
import matplotlib.pyplot as plt

from rdkit import Chem, RDLogger, rdBase

from rdkit.Chem import Lipinski, rdMolDescriptors, QED, Descriptors, Crippen
from rdkit.Chem import AllChem as Chem

from medchem.rules._utils import n_fused_aromatic_rings

from logger_config import logger
from configs.config_utils import load_config

RDLogger.DisableLog('rdApp.*')
rdBase.DisableLog('rdApp.*')
dm.disable_rdkit_log() 



def get_model_name(path=None, config=None, df=None, mode='single_comparison'):
    assert mode in ['single_comparison', 'multi_comparison'], "Mode must be either 'single_comparison' or 'multi_comparison'"
    assert config or path or (df is not None), "Config OR path OR df must be provided. Choose one of them."
    assert (config and not path and not (df is not None)) or (not config and path and not (df is not None)) or (not config and not path and (df is not None)), "Config OR path OR df must be provided. Provide only one."

    if config:  
        if mode == 'single_comparison':
            return config['generated_mols_path'].split('/')[-1].split('.')[0]
        
        else:
            paths = glob.glob(config['generated_mols_path'])
            model_names = [path.split('/')[-1].split('.')[0] for path in paths]
            return model_names
        

    if path:
        model_name = path.split('/')[-1].split('.')[0]
        return model_name
    

    if df is not None:
        if mode == 'single_comparison':
            model_name = config['generated_mols_path'].split('/')[-1].split('.')[0]
            return model_name
        
        elif mode == 'multi_comparison':
            model_names = df['model_name'].unique().tolist()
            return model_names
        else:
            raise ValueError(f"Invalid mode: {mode}")
        
    

def get_model_colors(model_names, cmap=None):
    return dict(zip(model_names, plt.cm.YlOrRd(np.linspace(1, 0, len(model_names) + 1)) if cmap is None else plt.colormaps.get_cmap(cmap)(np.linspace(1, 0, len(model_names) + 1))))


def drop_false_rows(df, borders):
    passed_cols = []
    filter_charged_mol = borders['filter_charged_mol']
    charged_mol_col = None

    for col in df.columns:
        if 'pass' in col:
            if 'charged_mol' in col:
                if filter_charged_mol:
                    passed_cols.append(col)
                else:
                    charged_mol_col = col
            else:
                passed_cols.append(col)
    
    mask = df[passed_cols].all(axis=1)
    df_masked = df[mask].copy()

    if not filter_charged_mol and charged_mol_col is not None and charged_mol_col not in df_masked.columns:
        df_masked[charged_mol_col] = df.loc[mask, charged_mol_col]
    return df_masked

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
            logger.error(f"Error processing value: {val}, error: {e}")
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
            logger.error(f"Error processing value: {val}, error: {e}")
    return parsed


def compute_metrics(df, save_path, config):
    mode = config['mode']
    if mode == 'single_comparison':
        model_name = get_model_name(config=config, mode=mode)
    else:
        model_name = get_model_name(df=df, mode=mode)
    metrics = {}
    
    skipped_molecules = []

    for smiles in df[df.columns[0]]:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol_metrics = {}

            symbols = list(set(atom.GetSymbol() 
                                for atom in mol.GetAtoms() 
                                if atom.GetSymbol()))

            charged_mol = False if any(atom.GetFormalCharge() != 0 for atom in mol.GetAtoms()) else True

            ring_info = mol.GetRingInfo()
            rings = [len(x) for x in ring_info.AtomRings()]

            total_bonds = mol.GetNumBonds()
            n_rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
            n_rigid_bonds = total_bonds - n_rot_bonds

            n_heavy_atoms = rdMolDescriptors.CalcNumHeavyAtoms(mol)
            n_aromatic_atoms = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic() and a.GetAtomicNum()>1)
            molWt = Descriptors.ExactMolWt(mol)
            clogp = Crippen.MolLogP(mol)


            if mode == 'multi_comparison': 
                model_name = df.loc[df['smiles'] == smiles, 'model_name'].iloc[0]
                mol_metrics['model_name'] = model_name

            mol_metrics['chars'] = symbols
            mol_metrics['n_atoms'] = Chem.AddHs(mol).GetNumAtoms()
            mol_metrics['n_heavy_atoms'] = n_heavy_atoms
            mol_metrics['n_N_atoms'] = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
            mol_metrics['fN_atoms'] = mol_metrics['n_N_atoms'] / mol_metrics['n_heavy_atoms']
            mol_metrics['charged_mol'] = charged_mol
            mol_metrics['molWt'] = molWt
            mol_metrics['logP'] = Descriptors.MolLogP(mol)
            mol_metrics['clogP'] = clogp
            mol_metrics['sw'] = 0.16 - 0.63*clogp - 0.0062*molWt + 0.066*n_rot_bonds - 0.74*n_aromatic_atoms
            mol_metrics['ring_size'] = rings
            mol_metrics['n_rings'] = mol.GetRingInfo().NumRings()
            mol_metrics['n_aroma_rings'] = rdMolDescriptors.CalcNumAromaticRings(mol)
            mol_metrics['n_fused_aromatic_rings'] = n_fused_aromatic_rings(mol)
            mol_metrics['n_aromatic_atoms'] = n_aromatic_atoms
            mol_metrics['n_het_atoms'] = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in (1, 6))
            mol_metrics['n_rigid_bonds'] = n_rigid_bonds
            mol_metrics['n_rot_bonds'] = n_rot_bonds
            mol_metrics['hbd'] = Lipinski.NumHDonors(mol)
            mol_metrics['hba'] = Lipinski.NumHAcceptors(mol)
            mol_metrics['fsp3'] = rdMolDescriptors.CalcFractionCSP3(mol)
            mol_metrics['tpsa'] = rdMolDescriptors.CalcTPSA(mol)
            mol_metrics['qed'] = QED.qed(mol)
            metrics[smiles] = mol_metrics
            
        else:
            if mode == 'single_comparison':
                skipped_molecules.append(smiles)
            else:
                skipped_molecules.append((smiles, model_name))

    if not save_path.endswith('/'):
        save_path = save_path + f'/'

    save_path = save_path + f'Descriptors/'

    if skipped_molecules:
        logger.warning(f'Skipped {len(skipped_molecules)} molecules: {skipped_molecules}')
        if mode == 'single_comparison':
            skipped_df = pd.DataFrame({'smiles': skipped_molecules})
            skipped_df.to_csv(save_path + f'skippedMolsDescriptors.csv', index=False)
        else:
            skipped_df = pd.DataFrame({'smiles': [smiles for smiles, _ in skipped_molecules], 
                                       'model_name': [model_name for _, model_name in skipped_molecules]})
            skipped_df.to_csv(save_path + f'skippedMolsDescriptors.csv', index=False)
        
    folder_to_save = save_path + f'perMolMetrics.csv'
    metrics_df = pd.DataFrame.from_dict(metrics, orient='index').reset_index().rename(columns={'index': 'smiles'})
    metrics_df.to_csv(folder_to_save, index_label=model_name, index=False)
    return metrics_df  


def filter_molecules(df, borders, folder_to_save, mode):
    folder_to_save = folder_to_save + f'Descriptors/'
    logger.info(f'Borders: {borders}')
    filtered_data = {}
    for col in df.columns.tolist():
        col_in_borders = any(col.lower() in k.lower() for k in borders.keys())

        if col == 'smiles' or col == 'model_name':
            filtered_data[col] = df[col]

        elif col_in_borders:
            filtered_data[col] = df[col]

            relevant_keys = [k for k in borders.keys() if col.lower() in k.lower()]
            min_border = next((borders[k] for k in relevant_keys if 'min' in k), None)
            max_border = next((borders[k] for k in relevant_keys if 'max' in k), None)

            if col == 'chars':
                filtered_data[f'{col}_passed'] = df[col].apply(lambda x: 
                                                                all(str(char).strip() in borders['allowed_chars'] 
                                                                    for char in (x if isinstance(x, list) else eval(x))))
            elif col == 'ring_size':
                filtered_data[f'{col}_passed'] = df[col].apply(lambda x: 
                                                                all(float(ring_size) >= min_border and float(ring_size) <= max_border 
                                                                    for ring_size in (x if isinstance(x, list) else eval(x))))
            elif col == 'syba_score':
                if max_border == 'inf':
                    filtered_data[f'{col}_passed'] = df[col] >= min_border
                else:
                    filtered_data[f'{col}_passed'] = (df[col] >= min_border) & (df[col] <= max_border)
            
            else:
                filtered_data[f'{col}_passed'] = (df[col] >= min_border) & (df[col] <= max_border)

    filtered_data_withFalse = pd.DataFrame(filtered_data)
    filtered_data_withFalse.to_csv(folder_to_save + 'filterResultsAfterDescriptors.csv', index_label='SMILES', index=False)
    left_mol_after_filters = drop_false_rows(filtered_data_withFalse, borders)

    if len(left_mol_after_filters) > 0:
        left_mol_after_filters.to_csv(folder_to_save + 'leftMolsAfterDescriptorsMetircs.csv', index_label='SMILES', index=False)
        if mode == 'single_comparison':
            left_mol_after_filters['smiles'].to_csv(folder_to_save + 'leftMolsAfterDescriptorsSMILES.csv', index_label='SMILES', index=False)
        else:
            left_mol_after_filters[['smiles', 'model_name']].to_csv(folder_to_save + 'leftMolsAfterDescriptorsSMILES.csv', index_label='SMILES', index=False)
    else:
        logger.warning(f'No molecules passed Descriptors Filters')
    return



def draw_filtered_mols(df, folder_to_save, config): 
    mode = config['mode']
    if mode == 'single_comparison':
        model_name = get_model_name(config=config, mode=mode)
    else:   
        model_name = get_model_name(df=df, mode=mode)
    folder_to_save = folder_to_save + '/Descriptors/'
    
    descriptors_config = load_config(config['config_descriptors'])
    borders = descriptors_config['borders']
    borders['charged_mol_allowed'] = int(borders['charged_mol_allowed'])
    cols_to_plot = descriptors_config['filtered_cols_to_plot']
    discrete_feats = descriptors_config['discrete_features_to_plot']
    not_to_smooth_by_sides_cols = descriptors_config['not_to_smooth_plot_by_sides']
    renamer = descriptors_config['renamer']

    if mode == 'multi_comparison':
        model_names = df['model_name'].unique().tolist()
        colors = get_model_colors(model_names, cmap='tab20')
    else:
        model_names = [model_name]
        colors = get_model_colors(model_names, cmap='tab20')

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

        ax = axes[i]
        for model in model_names:
            if mode == 'multi_comparison':
                model_df = df[df['model_name'] == model]
            else:
                model_df = df
            if col == 'chars':
                values_raw = model_df[col].dropna()
                values_before = _parse_chars_in_mol_column(values_raw)
            elif col == 'ring_size':
                values_raw = model_df[col].dropna()
                values_before = _parse_ring_size_column(values_raw)
            else:
                values_before = model_df[col].dropna().tolist()

            if max_val == 'inf':
                values_after = [v for v in values_before if (min_val is None or v >= min_val)]
            else:
                values_after = [v for v in values_before if (min_val is None or v >= min_val) and (max_val is None or v <= max_val)]

            if min_val is None:
                values_after = [v for v in values_before if (max_val is None or v <= max_val)]
            else:
                values_after = [v for v in values_before if (min_val is None or v >= min_val) and (max_val is None or v <= max_val)]

            total = len(values_before)
            count = len(values_after)
            mols_passed = count / total * 100 if total > 0 else 0

            label = f'{model}, passed: {mols_passed:.1f}%' 
            color = colors[model]

            if len(values_before) > 1:
                if col in discrete_feats:
                    if col == 'charged_mol':
                        value_counts = pd.Series(values_before).value_counts().sort_index()
                        complete_counts = pd.Series([0, 0], index=[False, True])
                        complete_counts.update(value_counts)
                        value_counts = complete_counts.sort_index()
                        ax.bar(x=['Not charged', 'Charged'], height=value_counts.values, 
                            alpha=0.5, color=color, edgecolor='black', linewidth=0.3, 
                            label=label)
                    else:
                        value_counts = pd.Series(values_before).value_counts().sort_index()
                        ax.bar(value_counts.index, value_counts.values, alpha=0.5, color=color, edgecolor='black', linewidth=0.3, align='edge', width=0.8, label=label)
                elif col in not_to_smooth_by_sides_cols:
                    sns.kdeplot(values_before, label=label, fill=True, alpha=0.3, ax=ax, color=color, clip=(0, None))
                else:
                    sns.kdeplot(values_before, label=label, fill=True, alpha=0.3, ax=ax, color=color)
            else:
                ax.scatter(values_before, [0.01]*len(values_before), label=label, alpha=0.7, color=color)

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

        ax.legend(fontsize=8, loc='upper right')
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_linewidth(0.5)
        ax.spines['left'].set_linewidth(0.5)

        ax.tick_params(axis='both', which='both', 
                      bottom=True, top=False, left=True, right=False,
                      labelbottom=True, labeltop=False, 
                      labelleft=True, labelright=False,
                      length=4, width=0.5, colors='black', labelsize=10)

    for j in range(len(cols_to_plot), len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.subplots_adjust(top=0.93, bottom=0.05, left=0.05, right=0.98)
    plt.savefig(f"{folder_to_save}filteredMetricsDistribution.png", dpi=300, bbox_inches='tight', format='png')
    
    return
