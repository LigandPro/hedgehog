import os 
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


def get_model_name(path=None, config=None, df=None):
    assert config or path or (df is not None), "Config OR path OR df must be provided. Choose one of them."
    assert (config and not path and not (df is not None)) or (not config and path and not (df is not None)) or (not config and not path and (df is not None)), "Config OR path OR df must be provided. Provide only one."

    if config:
        paths = glob.glob(config['generated_mols_path'])
        if len(paths) > 1:
            return [p.split('/')[-1].split('.')[0] for p in paths]
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

    if path:
        return path.split('/')[-1].split('.')[0]
    
    if df is not None:
        if 'model_name' in df.columns and df['model_name'].nunique(dropna=True) > 1:
            return sorted(df['model_name'].dropna().unique().tolist())
        return 'single'
    

def process_path(folder_to_save, key_word=None):
    if not folder_to_save.endswith('/'):
        folder_to_save = folder_to_save + '/'

    if key_word:
        folder_to_save = folder_to_save + f'{key_word}/'
    os.makedirs(folder_to_save, exist_ok=True)
    return folder_to_save


def backfill_identity(df_to_fill: pd.DataFrame, src_df: pd.DataFrame) -> pd.DataFrame:
    """Preserve existing mol_idx and backfill missing values from src_df identities.
    src_df should contain at least 'smiles' and optionally 'model_name', 'mol_idx'.
    """
    if 'smiles' not in df_to_fill.columns:
        return df_to_fill
    if 'smiles' not in src_df.columns:
        smiles_col = next((c for c in src_df.columns if c.lower() == 'smiles'), None)
        if smiles_col:
            src_df = src_df.rename(columns={smiles_col: 'smiles'})
        else:
            return df_to_fill
    right = src_df.copy()
    if 'model_name' not in right.columns:
        right['model_name'] = 'single'
    left = df_to_fill.copy()
    if 'model_name' not in left.columns:
        left['model_name'] = 'single'
    merged = left.merge(right[['smiles','model_name'] + ([ 'mol_idx'] if 'mol_idx' in right.columns else [])],
                       on=['smiles','model_name'], how='left', suffixes=('', '_id'))
    if 'mol_idx_id' in merged.columns:
        if 'mol_idx' in merged.columns:
            merged['mol_idx'] = merged['mol_idx'].fillna(merged['mol_idx_id'])
        else:
            merged['mol_idx'] = merged['mol_idx_id']
        merged.drop(columns=['mol_idx_id'], inplace=True)
    return merged

def drop_false_rows(df, borders):
    passed_cols = []
    filter_charged_mol = borders['filter_charged_mol'] if 'filter_charged_mol' in borders else False
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


def compute_metrics(df, save_path, config=None):
    is_multi = 'model_name' in df.columns and df['model_name'].nunique(dropna=True) > 1
    model_name = get_model_name(df=df) if is_multi else (get_model_name(config=config) if config else get_model_name(df=df))
    metrics = {}
    skipped_molecules = []

    has_mol_idx = 'mol_idx' in df.columns
    has_model_name_col = 'model_name' in df.columns
    for smiles in df[df.columns[0]]:
        mol_n = Chem.MolFromSmiles(smiles)
        if mol_n:
            mol_metrics = {}
            mol = Chem.AddHs(mol_n)
            symbols = list(set(atom.GetSymbol() 
                                for atom in mol.GetAtoms() 
                                if atom.GetSymbol()))

            charged_mol = False if any(atom.GetFormalCharge() != 0 for atom in mol.GetAtoms()) else True

            ring_info = mol.GetRingInfo()
            rings = [len(x) for x in ring_info.AtomRings()]

            total_bonds = mol_n.GetNumBonds()
            n_rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol_n)
            n_rigid_bonds = total_bonds - n_rot_bonds

            n_heavy_atoms = rdMolDescriptors.CalcNumHeavyAtoms(mol_n)
            n_aromatic_atoms = sum(1 for a in mol_n.GetAtoms() if a.GetIsAromatic() and a.GetAtomicNum()>1)
            molWt = Descriptors.ExactMolWt(mol_n)
            clogp = Crippen.MolLogP(mol_n)


            if has_model_name_col:
                model_for_smiles = df.loc[df['smiles'] == smiles, 'model_name'].iloc[0]
            else:
                model_for_smiles = model_name
            mol_metrics['model_name'] = model_for_smiles
            if has_mol_idx:
                if has_model_name_col:
                    sel = df[(df['smiles'] == smiles) & (df['model_name'] == model_for_smiles)]
                    if not sel.empty and 'mol_idx' in sel.columns:
                        mol_metrics['mol_idx'] = sel['mol_idx'].iloc[0]
                else:
                    mol_idx_val = df.loc[df['smiles'] == smiles, 'mol_idx'].iloc[0]
                    mol_metrics['mol_idx'] = mol_idx_val

            mol_metrics['chars'] = symbols
            mol_metrics['n_atoms'] = Chem.AddHs(mol).GetNumAtoms()
            mol_metrics['n_heavy_atoms'] = n_heavy_atoms
            mol_metrics['n_het_atoms'] = sum(1 for atom in mol_n.GetAtoms() if atom.GetAtomicNum() not in (1, 6))
            mol_metrics['n_N_atoms'] = sum(1 for atom in mol_n.GetAtoms() if atom.GetAtomicNum() == 7)
            mol_metrics['fN_atoms'] = mol_metrics['n_N_atoms'] / mol_metrics['n_heavy_atoms']
            mol_metrics['charged_mol'] = charged_mol
            mol_metrics['molWt'] = molWt
            mol_metrics['logP'] = Descriptors.MolLogP(mol_n)
            mol_metrics['clogP'] = clogp
            mol_metrics['sw'] = 0.16 - 0.63*clogp - 0.0062*molWt + 0.066*n_rot_bonds - 0.74*n_aromatic_atoms
            mol_metrics['ring_size'] = rings
            mol_metrics['n_rings'] = mol_n.GetRingInfo().NumRings()
            mol_metrics['n_aroma_rings'] = rdMolDescriptors.CalcNumAromaticRings(mol_n)
            mol_metrics['n_fused_aromatic_rings'] = n_fused_aromatic_rings(mol_n)
            mol_metrics['n_aromatic_atoms'] = n_aromatic_atoms
            mol_metrics['n_rigid_bonds'] = n_rigid_bonds
            mol_metrics['n_rot_bonds'] = n_rot_bonds
            mol_metrics['hbd'] = Lipinski.NumHDonors(mol_n)
            mol_metrics['hba'] = Lipinski.NumHAcceptors(mol_n)
            mol_metrics['fsp3'] = rdMolDescriptors.CalcFractionCSP3(mol_n)
            mol_metrics['tpsa'] = rdMolDescriptors.CalcTPSA(mol_n)
            mol_metrics['qed'] = QED.qed(mol_n)
            metrics[smiles] = mol_metrics
            
        else:
            if is_multi:
                skipped_molecules.append((smiles, df.loc[df['smiles'] == smiles, 'model_name'].iloc[0]))
            else:
                skipped_molecules.append(smiles)

    save_path = process_path(save_path, key_word='Descriptors')

    if skipped_molecules:
        logger.warning(f'Skipped {len(skipped_molecules)} molecules: {skipped_molecules}')
        if is_multi:
            skipped_df = pd.DataFrame({'smiles': [s for s, _ in skipped_molecules], 
                                       'model_name': [m for _, m in skipped_molecules]})
        else:
            skipped_df = pd.DataFrame({'smiles': skipped_molecules})
            skipped_df['model_name'] = model_name
        if 'mol_idx' in df.columns:
            skipped_df = skipped_df.merge(df[['smiles', 'mol_idx']].drop_duplicates('smiles'), on='smiles', how='left')
        skipped_df.to_csv(save_path + f'skippedMolsDescriptors.csv', index=False)
        
    folder_to_save = save_path + f'perMoleculeDescriptors.csv'
    metrics_df = pd.DataFrame.from_dict(metrics, orient='index').reset_index().rename(columns={'index': 'smiles'})

    if 'model_name' not in metrics_df.columns:
        metrics_df['model_name'] = model_name
    metrics_df = backfill_identity(metrics_df, df)

    id_cols = ['smiles', 'model_name', 'mol_idx']
    for col in ['model_name', 'mol_idx']:
        if col not in metrics_df.columns:
            if col == 'model_name':
                metrics_df[col] = model_name
    ordered = id_cols + [c for c in metrics_df.columns if c not in id_cols]
    metrics_df = metrics_df[ordered]
    metrics_df.to_csv(folder_to_save, index=False)
    return metrics_df  


def filter_molecules(df, borders, folder_to_save):
    folder_to_save = process_path(folder_to_save, key_word='Descriptors')
    is_multi = 'model_name' in df.columns and df['model_name'].nunique(dropna=True) > 1
    has_mol_idx = 'mol_idx' in df.columns
    logger.info(f'Borders: {borders}')
    filtered_data = {}
    for col in df.columns.tolist():
        col_in_borders = any(col.lower() in k.lower() for k in borders.keys())

        if col == 'smiles' or col == 'model_name' or col == 'mol_idx':
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

    if 'model_name' not in filtered_data_withFalse.columns:
        filtered_data_withFalse['model_name'] = model_name

    filtered_data_withFalse = backfill_identity(filtered_data_withFalse, df)

    id_cols = ['smiles', 'model_name', 'mol_idx']
    for col in ['model_name', 'mol_idx']:
        if col not in filtered_data_withFalse.columns:
            if col == 'model_name':
                filtered_data_withFalse[col] = model_name
    ordered = id_cols + [c for c in filtered_data_withFalse.columns if c not in id_cols]
    filtered_data_withFalse = filtered_data_withFalse[ordered]
    filtered_data_withFalse.to_csv(folder_to_save + 'descriptorsPassFlags.csv', index_label='SMILES', index=False)
    pass_filters = drop_false_rows(filtered_data_withFalse, borders)

    if len(pass_filters) > 0:

        if 'model_name' not in pass_filters.columns:
            pass_filters['model_name'] = model_name

        id_cols = ['smiles', 'model_name', 'mol_idx']
        ordered_full = id_cols + [c for c in pass_filters.columns if c not in id_cols]
        pass_filters = pass_filters[ordered_full]
        pass_filters.to_csv(folder_to_save + 'passDescriptorsMetircs.csv', index_label='SMILES', index=False)
        cols = ['smiles', 'model_name'] + (['mol_idx'] if has_mol_idx else [])
        pass_filters[cols].to_csv(folder_to_save + 'passDescriptorsSMILES.csv', index_label='SMILES', index=False)
    else:
        logger.warning(f'No molecules passed Descriptors Filters')

    passed_flag_cols = [col for col in filtered_data_withFalse.columns if col.endswith('_passed')]
    if passed_flag_cols:
        all_pass_mask = filtered_data_withFalse[passed_flag_cols].all(axis=1)
        fail_filters = filtered_data_withFalse.loc[~all_pass_mask].copy()
        if len(fail_filters) > 0:

            if 'model_name' not in fail_filters.columns:
                fail_filters['model_name'] = model_name

            id_cols = ['smiles', 'model_name', 'mol_idx']
            ordered_full = id_cols + [c for c in fail_filters.columns if c not in id_cols]
            fail_filters = fail_filters[ordered_full]
            fail_filters.to_csv(folder_to_save + 'failDescriptorsMetircs.csv', index_label='SMILES', index=False)
            cols = ['smiles', 'model_name'] + (['mol_idx'] if has_mol_idx else [])
            fail_filters[cols].to_csv(folder_to_save + 'failDescriptorsSMILES.csv', index_label='SMILES', index=False)
    try:

        id_path = folder_to_save + 'passDescriptorsSMILES.csv'
        if os.path.exists(id_path):
            id_df = pd.read_csv(id_path)
            if 'model_name' not in id_df.columns:
                id_df['model_name'] = 'single'
            keep_cols = ['smiles', 'model_name'] + (['mol_idx'] if 'mol_idx' in id_df.columns else [])
            id_df = id_df[keep_cols]

            if 'model_name' in id_df.columns:
                id_df = id_df.drop_duplicates(['smiles','model_name'])
            else:
                id_df = id_df.drop_duplicates('smiles')

            per_path = folder_to_save + 'perMoleculeDescriptors.csv'
            if os.path.exists(per_path):
                per = pd.read_csv(per_path)
                merged = backfill_identity(per, id_df)
                identity_order = ['smiles', 'model_name', 'mol_idx']
                ordered = [c for c in identity_order if c in merged.columns] + [c for c in merged.columns if c not in identity_order]
                merged = merged[ordered]
                merged.to_csv(per_path, index=False)

            flags_path = folder_to_save + 'descriptorsPassFlags.csv'
            if os.path.exists(flags_path):
                flags = pd.read_csv(flags_path)
                merged = backfill_identity(flags, id_df)
                identity_order = ['smiles', 'model_name', 'mol_idx']
                ordered = [c for c in identity_order if c in merged.columns] + [c for c in merged.columns if c not in identity_order]
                merged = merged[ordered]
                merged.to_csv(flags_path, index=False)
    except Exception:
        pass
    return



def draw_filtered_mols(df, folder_to_save, config): 
    is_multi = 'model_name' in df.columns and df['model_name'].nunique(dropna=True) > 1
    model_name = get_model_name(df=df) if is_multi else get_model_name(config=config)
    folder_to_save = process_path(folder_to_save, key_word='Descriptors')
    
    descriptors_config = load_config(config['config_descriptors'])
    borders = descriptors_config['borders']
    borders['charged_mol_allowed'] = int(borders['charged_mol_allowed']) if 'charged_mol_allowed' in borders else False
    cols_to_plot = descriptors_config['filtered_cols_to_plot']
    discrete_feats = descriptors_config['discrete_features_to_plot']
    not_to_smooth_by_sides_cols = descriptors_config['not_to_smooth_plot_by_sides']
    renamer = descriptors_config['renamer']

    if is_multi:
        model_names = sorted(df['model_name'].dropna().unique().tolist())
        model_names = [model.lower() for model in model_names]
    else:
        model_names = [model_name]
    distinct_colors = ['brown', 'green', 'blue', 'cyan', 'yellow', 'pink', 'orange', 
                       '#dd37fa', '#ad5691', '#f46fa1', '#89cff0', '#93c83e', ]
    colors = {model:color for model, color in zip(sorted(model_names), distinct_colors)}

    nrows = len(cols_to_plot) // 5  + len(cols_to_plot) % 5 if len(cols_to_plot) > 5 else 1
    ncols = 5 if len(cols_to_plot) > 5 else len(cols_to_plot)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(nrows * ncols, nrows * ncols))
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
            if is_multi:
                model_df = df[df['model_name'].str.lower() == model]
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

            values_after = values_before.copy()
            if min_val is not None:
                values_after = [v for v in values_after if v >= min_val]
            if max_val is not None and max_val != 'inf':
                values_after = [v for v in values_after if v <= max_val]

            total = len(values_before)
            count = len(values_after)
            mols_passed = count / total * 100 if total > 0 else 0

            if 'model_name' in df.columns:
                name_map = {str(m): str(m).upper() for m in sorted(df['model_name'].dropna().unique())}
            else:
                name_map = {}
            name_map_lc = {str(k).lower(): v for k, v in name_map.items()}
            label_name = name_map_lc.get(str(model).lower(), str(model).upper())
            label = f'{label_name}, pass: {mols_passed:.1f}%'
            color = colors[model] 

            if len(values_before) > 1:
                if col in discrete_feats:
                    bar_width = 0.1 
                    model_index = model_names.index(model)
                    offset = model_index * bar_width

                    if col == 'chars':
                        value_counts = pd.Series(values_before).value_counts()
                        desired_order = ['C', 'N', 'S', 'O', 'F', 'Cl', 'Br', 'H']
                        all_chars = borders.get('allowed_chars', desired_order)
                        sorted_chars = [char for char in desired_order if char in all_chars] + [char for char in all_chars if char not in desired_order]
                        complete_counts = pd.Series(0, index=sorted_chars)
                        complete_counts.update(value_counts)
                        x_positions = [i + offset for i in range(len(complete_counts.index))]
                        ax.bar(x_positions, complete_counts.values, width=bar_width, alpha=0.4, color=color, 
                                edgecolor='black', linewidth=0.3, label=label)
                        if model_index == 0:
                            tick_positions = list(range(len(complete_counts.index)))
                            ax.set_xticks(tick_positions)
                            ax._discrete_tick_values = complete_counts.index
                            ax.set_xticklabels(complete_counts.index)
                    else:
                        value_counts = pd.Series(values_before).value_counts().sort_index()
                        
                        if max_val is not None and max_val != 'inf':
                            extended_max = int(max_val) + 5
                            full_range = list(range(0, extended_max + 1))
                            complete_counts = pd.Series(0, index=full_range)
                            for val in value_counts.index:
                                if val in complete_counts.index:
                                    complete_counts[val] = value_counts[val]
                        else:
                            complete_counts = value_counts
                            
                        x_positions = [i + offset for i in range(len(complete_counts.index))]
                        ax.bar(x_positions, complete_counts.values, width=bar_width, alpha=0.4, color=color, 
                                edgecolor='black', linewidth=0.3, label=label)
                        if model_index == 0:
                            if col == 'n_rigid_bonds':
                                all_values = list(complete_counts.index)
                                tick_values = [x for x in all_values if x % 5 == 0]
                                if max(all_values) not in tick_values:
                                    tick_values.append(max(all_values))
                                tick_positions = [all_values.index(val) for val in tick_values if val in all_values]
                                ax.set_xticks(tick_positions)
                                ax.set_xticklabels([str(int(val)) for val in tick_values])
                                ax._discrete_tick_values = complete_counts.index
                            else:
                                tick_positions = list(range(len(complete_counts.index)))
                                ax.set_xticks(tick_positions)
                                ax._discrete_tick_values = complete_counts.index
                                ax.set_xticklabels([str(int(x)) for x in complete_counts.index])
             
                elif col in not_to_smooth_by_sides_cols:
                    if col in ['fsp3', 'qed']:
                        sns.kdeplot(values_before, label=label, fill=True, alpha=0.4, ax=ax, color=color, clip=(0, 1.0))
                    else:
                        sns.kdeplot(values_before, label=label, fill=True, alpha=0.4, ax=ax, color=color, clip=(0, None))
                else:
                    if col in ['fsp3', 'qed']:
                        sns.kdeplot(values_before, label=label, fill=True, alpha=0.4, ax=ax, color=color, clip=(0, 1.0))
                    else:
                        sns.kdeplot(values_before, label=label, fill=True, alpha=0.4, ax=ax, color=color)
            else:
                ax.scatter(values_before, [0.01]*len(values_before), label=label, alpha=0.4, color=color)

        if col == 'chars':
            title = f"{renamer[col] if col in renamer else col}"
        else:   
            title = f"{renamer[col] if col in renamer else col} ({minmax_str})"
        ax.set_title(title, fontsize=12)
        ax.set_xlabel(renamer[col], fontsize=10)
        
        if col not in discrete_feats:
            if col == 'fsp3': ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
            elif col == 'qed': pass  
            else: ax.xaxis.set_major_locator(plt.MaxNLocator(integer=False)) 
        ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

        if min_val is not None:
            if col in discrete_feats:
                if hasattr(ax, '_discrete_tick_values'):
                    try:
                        min_pos = list(ax._discrete_tick_values).index(min_val)
                        ax.axvline(min_pos - 0.5, color='red', linestyle='--', linewidth=1.5, label=f'min: {min_val}')
                    except ValueError:
                        ax.axvline(0, color='red', linestyle='--', linewidth=1.5, label=f'min: {min_val}')
                else:
                    ax.axvline(min_val, color='red', linestyle='--', linewidth=1.5, label=f'min: {min_val}')
            else:
                ax.axvline(min_val, color='red', linestyle='--', linewidth=1.5, label=f'min: {min_val}')
        if max_val is not None and max_val != 'inf':
            if col in discrete_feats:
                if hasattr(ax, '_discrete_tick_values'):
                    try:
                        max_pos = list(ax._discrete_tick_values).index(max_val)
                        ax.axvline(max_pos + 0.5, color='blue', linestyle='--', linewidth=1.5, label=f'max: {max_val}')
                    except ValueError:
                        ax.axvline(len(ax._discrete_tick_values) - 1, color='blue', linestyle='--', linewidth=1.5, label=f'max: {max_val}')
                else:
                    if col == 'fsp3':
                        ax.axvline(max_val + 0.01, color='blue', linestyle='--', linewidth=1.5, label=f'max: {max_val}')
                    elif col == 'n_rigid_bonds':
                        ax.axvline(max_val + 0.0000001, color='blue', linestyle='--', linewidth=1.5, label=f'max: {max_val}')
                    else:
                        ax.axvline(max_val + 0.5, color='blue', linestyle='--', linewidth=1.5, label=f'max: {max_val}')
            else:
                ax.axvline(max_val, color='blue', linestyle='--', linewidth=1.5, label=f'max: {max_val}')
                
        x_min, x_max = ax.get_xlim()
        if min_val is not None:
            if col in discrete_feats:
                if hasattr(ax, '_discrete_tick_values'):
                    try:
                        min_pos = list(ax._discrete_tick_values).index(min_val)
                        ax.axvspan(min_pos - 0.5, x_min, color='grey', alpha=0.2, zorder=0)
                    except ValueError:
                        ax.axvspan(x_min, 0, color='grey', alpha=0.2, zorder=0)
                else:
                    if col == 'fsp3':
                        ax.axvspan(min_val - 0.01, x_max, color='grey', alpha=0.2, zorder=0)
                    elif col == 'n_rigid_bonds':
                        ax.axvspan(min_val - 0.0000001, x_max, color='grey', alpha=0.2, zorder=0)
                    else:
                        ax.axvspan(min_val - 0.5, x_max, color='grey', alpha=0.2, zorder=0)
            else:
                ax.axvspan(x_min, min_val, color='grey', alpha=0.2, zorder=0)
        if max_val is not None and max_val != 'inf':
            if col in discrete_feats:
                if hasattr(ax, '_discrete_tick_values'):
                    try:
                        max_pos = list(ax._discrete_tick_values).index(max_val)
                        ax.axvspan(max_pos + 0.5, x_max, color='grey', alpha=0.2, zorder=0)
                    except ValueError:
                        ax.axvspan(len(ax._discrete_tick_values) - 0.5, x_max, color='grey', alpha=0.2, zorder=0)
                else:
                    if col == 'fsp3':
                        ax.axvspan(max_val + 0.01, x_max, color='grey', alpha=0.2, zorder=0)
                    elif col == 'n_rigid_bonds':
                        ax.axvspan(max_val + 0.0000001, x_max, color='grey', alpha=0.2, zorder=0)
                    else:
                        ax.axvspan(max_val + 0.5, x_max, color='grey', alpha=0.2, zorder=0)
            else:
                ax.axvspan(max_val, x_max, color='grey', alpha=0.2, zorder=0)

        handles, labels = ax.get_legend_handles_labels()
        sorted_labels, sorted_handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
        ax.legend(sorted_handles, sorted_labels, fontsize=8, loc='upper right')
        
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
    