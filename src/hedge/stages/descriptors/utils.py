import os
import json
import ast

import numpy as np
import polars as pl
import pandas as pd  # Keep for plotting compatibility
import seaborn as sns
import matplotlib.pyplot as plt

from rdkit import Chem, RDLogger, rdBase

from rdkit.Chem import Lipinski, rdMolDescriptors, QED, Descriptors, Crippen
from rdkit.Chem import AllChem as Chem

from medchem.rules._utils import n_fused_aromatic_rings

from hedge.configs.logger import logger, load_config

# Disable RDKit warnings
RDLogger.DisableLog('rdApp.*')
rdBase.DisableLog('rdApp.*')


def _stringify_nested_value(value):
    """Convert nested polars values (Series, list, struct) to a JSON string."""
    if value is None:
        return None
    if hasattr(value, 'to_list'):
        value = value.to_list()
    return json.dumps(value)


def write_csv_safe(df: pl.DataFrame, path: str) -> None:
    """
    Write a DataFrame to CSV, stringifying nested columns first so Polars can serialize them.
    """
    nested_cols = [
        name for name, dtype in zip(df.columns, df.dtypes)
        if isinstance(dtype, (pl.List, pl.Struct, pl.Array))
    ]
    if nested_cols:
        df = df.with_columns(
            [
                pl.col(name).map_elements(
                    _stringify_nested_value, return_dtype=pl.Utf8
                ).alias(name)
                for name in nested_cols
            ]
        )
    df.write_csv(path)


def _chars_pass(value, allowed_chars) -> bool:
    if value is None:
        return False
    if isinstance(value, list):
        items = value
    else:
        try:
            parsed = ast.literal_eval(value)
            items = parsed if isinstance(parsed, (list, tuple, set)) else [parsed]
        except Exception:
            items = [value]
    try:
        return all(str(char).strip() in allowed_chars for char in items)
    except Exception:
        return False


def _ring_size_pass(value, min_border, max_border) -> bool:
    if value is None:
        return False
    if isinstance(value, list):
        items = value
    else:
        try:
            parsed = ast.literal_eval(value)
            items = parsed if isinstance(parsed, (list, tuple, set)) else [parsed]
        except Exception:
            items = [value]

    def _within_bounds(size):
        try:
            size_val = float(size)
        except Exception:
            return False
        if min_border is not None and size_val < float(min_border):
            return False
        if max_border is not None and max_border != 'inf' and size_val > float(max_border):
            return False
        return True

    try:
        return all(_within_bounds(size) for size in items)
    except Exception:
        return False


def process_path(folder_to_save, key_word=None):
    """
    Args:
        folder_to_save: Base folder path
        key_word: Optional subfolder name
        
    Returns:
        str: Processed path with trailing slash
    """
    if not folder_to_save.endswith('/'):
        folder_to_save = folder_to_save + '/'

    if key_word:
        folder_to_save = folder_to_save + f'{key_word}/'
    os.makedirs(folder_to_save, exist_ok=True)
    return folder_to_save


def order_identity_columns(df):
    """Reorder dataframe columns with identity columns first."""
    id_cols = ['smiles', 'model_name', 'mol_idx']
    ordered = id_cols + [c for c in df.columns if c not in id_cols]
    return df.select(ordered)


def drop_false_rows(df, borders):
    """
    Filter rows that passed all descriptor filters.

    Args:
        df: DataFrame with '_pass' columns
        borders: Configuration dict with filter settings

    Returns:
        pl.DataFrame: Filtered dataframe with only passed molecules
    """
    passed_cols = []
    filter_charged_mol = borders.get('filter_charged_mol', False)
    charged_mol_col = None

    for col in df.columns:
        if col.endswith('_pass') or col == 'pass':
            if 'charged_mol' in col:
                if filter_charged_mol:
                    passed_cols.append(col)
                else:
                    charged_mol_col = col
            else:
                passed_cols.append(col)

    # Create mask by combining all pass columns with AND
    if passed_cols:
        mask = pl.lit(True)
        for col in passed_cols:
            mask = mask & pl.col(col)
        df_masked = df.filter(mask)
    else:
        df_masked = df.clone()

    if not filter_charged_mol and charged_mol_col is not None and charged_mol_col not in df_masked.columns:
        df_masked = df_masked.with_columns(df.filter(mask).select(charged_mol_col))
    return df_masked


def _parse_chars_in_mol_column(series):
    """Parse character lists from series for plotting."""
    parsed = []
    for val in series.dropna():
        try:
            chars = val if isinstance(val, list) else eval(val)
            parsed.extend(chars)
        except Exception as e:
            logger.error(f"Error parsing chars: {val}, {e}")
    return parsed


def _parse_ring_size_column(series):
    """Parse ring size lists from series for plotting."""
    parsed = []
    for val in series.dropna():
        try:
            sizes = val if isinstance(val, list) else eval(val)
            parsed.extend([float(size) for size in sizes])
        except Exception as e:
            logger.error(f"Error parsing ring sizes: {val}, {e}")
    return parsed


def compute_metrics(df, save_path, config=None):
    """
    Compute 22 physicochemical descriptors for each molecule.
    model_name and mol_idx are already in df from sampledMols.csv.

    Args:
        df: DataFrame with molecules (must have 'smiles', 'model_name', 'mol_idx')
        save_path: Output folder path
        config: Configuration dictionary

    Returns:
        pl.DataFrame: Dataframe with computed descriptors per molecule
    """
    if df is None or len(df) == 0:
        logger.warning('Empty DataFrame provided to compute_metrics. Returning empty DataFrame.')
        return pl.DataFrame()

    metrics = {}
    skipped_molecules = []

    for row in df.iter_rows(named=True):
        smiles = row['smiles']
        model_name = row['model_name']
        mol_idx = row['mol_idx']
        
        mol_n = Chem.MolFromSmiles(smiles)
        if mol_n:
            mol_metrics = {}
            mol = Chem.AddHs(mol_n)
            symbols = list(set(atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetSymbol()))
            charged_mol = False if any(atom.GetFormalCharge() != 0 for atom in mol.GetAtoms()) else True

            ring_info = mol.GetRingInfo()
            rings = [len(x) for x in ring_info.AtomRings()]

            n_rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol_n)
            n_rigid_bonds = mol_n.GetNumBonds() - n_rot_bonds
            n_heavy_atoms = rdMolDescriptors.CalcNumHeavyAtoms(mol_n)
            n_aromatic_atoms = sum(1 for a in mol_n.GetAtoms() if a.GetIsAromatic() and a.GetAtomicNum() > 1)
            molWt = Descriptors.ExactMolWt(mol_n)
            clogp = Crippen.MolLogP(mol_n)

            mol_metrics['model_name'] = model_name
            mol_metrics['mol_idx'] = mol_idx
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
            mol_metrics['n_rigid_bonds'] = n_rigid_bonds
            mol_metrics['n_rot_bonds'] = n_rot_bonds
            mol_metrics['hbd'] = Lipinski.NumHDonors(mol_n)
            mol_metrics['hba'] = Lipinski.NumHAcceptors(mol_n)
            mol_metrics['fsp3'] = rdMolDescriptors.CalcFractionCSP3(mol_n)
            mol_metrics['tpsa'] = rdMolDescriptors.CalcTPSA(mol_n)
            mol_metrics['qed'] = QED.qed(mol_n)
            metrics[smiles] = mol_metrics
        else:
            skipped_molecules.append((smiles, model_name, mol_idx))

    save_path = process_path(save_path)
    if skipped_molecules:
        logger.warning(f'Skipped {len(skipped_molecules)} molecules that failed to parse')
        skipped_data = {
            'smiles': [s for s, _, _ in skipped_molecules],
            'model_name': [m for _, m, _ in skipped_molecules]
        }
        if any(idx is not None for _, _, idx in skipped_molecules):
            skipped_data['mol_idx'] = [idx for _, _, idx in skipped_molecules]
        skipped_df = pl.DataFrame(skipped_data)
        write_csv_safe(skipped_df, save_path + 'skipped_molecules.csv')

    # Convert dict to list of records
    metrics_list = [{'smiles': smiles, **metrics_dict} for smiles, metrics_dict in metrics.items()]
    metrics_df = pl.DataFrame(metrics_list)
    metrics_df = order_identity_columns(metrics_df)
    write_csv_safe(metrics_df, save_path + 'descriptors_all.csv')
    return metrics_df  


def filter_molecules(df, borders, folder_to_save):
    """
    Filter molecules based on descriptor thresholds.
    
    Args:
        df: DataFrame with computed descriptors
        borders: Dictionary with min/max thresholds for each descriptor
        folder_to_save: Output folder path (should already include 'Descriptors' subfolder)
    """
    folder_to_save = process_path(folder_to_save)
    logger.info(f'Borders: {borders}')
    identity_cols = ['smiles', 'model_name', 'mol_idx']
    allowed_chars = borders.get('allowed_chars', [])

    select_exprs = []
    descriptor_exprs = []

    for col in df.columns:
        col_lower = col.lower()
        col_in_borders = any(col_lower in key.lower() for key in borders.keys())

        if col in identity_cols:
            select_exprs.append(pl.col(col))
            continue

        if not col_in_borders:
            continue

        select_exprs.append(pl.col(col))
        relevant_keys = [k for k in borders.keys() if col_lower in k.lower()]
        min_border = next((borders[k] for k in relevant_keys if 'min' in k), None)
        max_border = next((borders[k] for k in relevant_keys if 'max' in k), None)
        pass_col = f'{col}_pass'

        if col == 'chars' and allowed_chars:
            descriptor_exprs.append(
                pl.col(col)
                .map_elements(lambda x, allowed=allowed_chars: _chars_pass(x, allowed), return_dtype=pl.Boolean)
                .fill_null(False)
                .alias(pass_col)
            )
        elif col == 'ring_size':
            descriptor_exprs.append(
                pl.col(col)
                .map_elements(
                    lambda x, min_b=min_border, max_b=max_border: _ring_size_pass(x, min_b, max_b),
                    return_dtype=pl.Boolean,
                )
                .fill_null(False)
                .alias(pass_col)
            )
        elif col == 'syba_score':
            expr = pl.lit(True)
            if min_border is not None:
                expr = expr & (pl.col(col) >= min_border)
            if max_border is not None and max_border != 'inf':
                expr = expr & (pl.col(col) <= max_border)
            descriptor_exprs.append(expr.cast(pl.Boolean).fill_null(False).alias(pass_col))
        else:
            expr = pl.lit(True)
            if min_border is not None:
                expr = expr & (pl.col(col) >= min_border)
            if max_border is not None and max_border != 'inf':
                expr = expr & (pl.col(col) <= max_border)
            descriptor_exprs.append(expr.cast(pl.Boolean).fill_null(False).alias(pass_col))

    if select_exprs:
        filtered_data_withFalse = df.select(select_exprs)
    else:
        filtered_data_withFalse = df.clone()

    if descriptor_exprs:
        filtered_data_withFalse = filtered_data_withFalse.with_columns(descriptor_exprs)

    filtered_data_withFalse = order_identity_columns(filtered_data_withFalse)
    write_csv_safe(filtered_data_withFalse, folder_to_save + 'pass_flags.csv')
    pass_filters = drop_false_rows(filtered_data_withFalse, borders)

    if len(pass_filters) > 0:
        pass_filters = order_identity_columns(pass_filters)

        id_cols = ['smiles', 'model_name', 'mol_idx']
        descriptor_cols = [col for col in pass_filters.columns if col not in id_cols and not col.endswith('_pass')]

        ordered_cols = id_cols + sorted(descriptor_cols)
        ordered_cols = [col for col in ordered_cols if col in pass_filters.columns]
        write_csv_safe(pass_filters.select(ordered_cols), folder_to_save + 'descriptors_passed.csv')

        cols = [col for col in ['smiles', 'model_name', 'mol_idx'] if col in pass_filters.columns]
        if cols:
            write_csv_safe(pass_filters.select(cols), folder_to_save + 'filtered_molecules.csv')
    else:
        logger.warning('No molecules pass Descriptors Filters')

    if len(pass_filters) > 0:
        all_computed_path = folder_to_save + 'descriptors_all.csv'
        if os.path.exists(all_computed_path):
            all_computed = pl.read_csv(all_computed_path)
            merge_cols = ['smiles', 'model_name']
            merged = all_computed.join(pass_filters.select(merge_cols), on=merge_cols, how='left', suffix='_pass')
            fail_filters = merged.filter(pl.col('smiles_pass').is_null()).drop([c for c in merged.columns if c.endswith('_pass')])

            if len(fail_filters) > 0:
                flags_path = folder_to_save + 'pass_flags.csv'
                if os.path.exists(flags_path):
                    flags_df = pl.read_csv(flags_path)
                    merge_cols = ['smiles', 'model_name']
                    pass_cols = [col for col in flags_df.columns if col.endswith('_pass') or col == 'pass']
                    if pass_cols:
                        fail_filters = fail_filters.join(flags_df.select(merge_cols + pass_cols), on=merge_cols, how='left', suffix='_flags')
                        for col in pass_cols:
                            if f'{col}_flags' in fail_filters.columns:
                                fail_filters = fail_filters.with_columns(
                                    pl.col(f'{col}_flags').fill_null(pl.col(col) if col in fail_filters.columns else pl.lit(False)).alias(col)
                                ).drop(f'{col}_flags')
                
                fail_filters = order_identity_columns(fail_filters)
                
                id_cols = ['smiles', 'model_name', 'mol_idx']
                descriptor_cols = [col for col in fail_filters.columns if col not in id_cols and not col.endswith('_pass')]
                pass_cols = [col for col in fail_filters.columns if col.endswith('_pass')]
                
                ordered_cols = id_cols.copy()
                for desc_col in sorted(descriptor_cols):
                    if desc_col in fail_filters.columns:
                        ordered_cols.append(desc_col)
                        pass_col = f'{desc_col}_pass'
                        if pass_col in pass_cols:
                            ordered_cols.append(pass_col)
                
                for pass_col in sorted(pass_cols):
                    if pass_col not in ordered_cols:
                        ordered_cols.append(pass_col)
                
                ordered_cols = [col for col in ordered_cols if col in fail_filters.columns]
                write_csv_safe(fail_filters.select(ordered_cols), folder_to_save + 'descriptors_failed.csv')

                id_cols_smiles = ['smiles', 'model_name', 'mol_idx']
                write_csv_safe(fail_filters.select(id_cols_smiles), folder_to_save + 'failed_molecules.csv')
    else:
        all_computed_path = folder_to_save + 'descriptors_all.csv'
        if os.path.exists(all_computed_path):
            all_computed = pl.read_csv(all_computed_path)
            fail_filters = all_computed.clone()
            flags_path = folder_to_save + 'pass_flags.csv'
            if os.path.exists(flags_path):
                flags_df = pl.read_csv(flags_path)
                merge_cols = ['smiles', 'model_name']
                pass_cols = [col for col in flags_df.columns if col.endswith('_pass') or col == 'pass']
                if pass_cols:
                    fail_filters = fail_filters.join(flags_df.select(merge_cols + pass_cols), on=merge_cols, how='left', suffix='_flags')
                    for col in pass_cols:
                        if f'{col}_flags' in fail_filters.columns:
                            fail_filters = fail_filters.with_columns(
                                pl.col(f'{col}_flags').fill_null(pl.col(col) if col in fail_filters.columns else pl.lit(False)).alias(col)
                            ).drop(f'{col}_flags')
            
            fail_filters = order_identity_columns(fail_filters)
            
            id_cols = ['smiles', 'model_name', 'mol_idx']
            descriptor_cols = [col for col in fail_filters.columns if col not in id_cols and not col.endswith('_pass')]
            pass_cols = [col for col in fail_filters.columns if col.endswith('_pass')]
            
            ordered_cols = id_cols.copy()
            for desc_col in sorted(descriptor_cols):
                if desc_col in fail_filters.columns:
                    ordered_cols.append(desc_col)
                    pass_col = f'{desc_col}_pass'
                    if pass_col in pass_cols:
                        ordered_cols.append(pass_col)
            
            for pass_col in sorted(pass_cols):
                if pass_col not in ordered_cols:
                    ordered_cols.append(pass_col)
            
            ordered_cols = [col for col in ordered_cols if col in fail_filters.columns]
            write_csv_safe(fail_filters.select(ordered_cols), folder_to_save + 'descriptors_failed.csv')

            id_cols_smiles = ['smiles', 'model_name', 'mol_idx']
            write_csv_safe(fail_filters.select(id_cols_smiles), folder_to_save + 'failed_molecules.csv')

    id_path = folder_to_save + 'filtered_molecules.csv'
    if os.path.exists(id_path):
        id_df = pl.read_csv(id_path)
        id_df = id_df.unique(subset=['smiles', 'model_name'])

        per_path = folder_to_save + 'descriptors_all.csv'
        if os.path.exists(per_path):
            per = pl.read_csv(per_path)
            per = order_identity_columns(per)
            write_csv_safe(per, per_path)

        flags_path = folder_to_save + 'pass_flags.csv'
        if os.path.exists(flags_path):
            flags = pl.read_csv(flags_path)
            flags = order_identity_columns(flags)
            write_csv_safe(flags, flags_path)


def draw_filtered_mols(df, folder_to_save, config):
    """
    Generate distribution plots for descriptor filters.

    Args:
        df: DataFrame with computed descriptors
        folder_to_save: Output folder path (should already include 'Descriptors' subfolder)
        config: Configuration dictionary
    """
    # Convert polars to pandas for plotting compatibility
    if isinstance(df, pl.DataFrame):
        df = pd.DataFrame(df.to_dict(as_series=False))

    is_multi = df['model_name'].nunique(dropna=True) > 1
    folder_to_save = process_path(folder_to_save)

    descriptors_config = load_config(config['config_descriptors'])
    borders = descriptors_config['borders']
    borders['charged_mol_allowed'] = int(borders['charged_mol_allowed']) if 'charged_mol_allowed' in borders else False
    cols_to_plot = descriptors_config['filtered_cols_to_plot']
    discrete_feats = descriptors_config['discrete_features_to_plot']
    not_to_smooth_by_sides_cols = descriptors_config['not_to_smooth_plot_by_sides']
    renamer = descriptors_config['renamer']

    model_names = sorted(df['model_name'].dropna().unique().tolist())
    model_names = [m.lower() for m in model_names]
    is_multi = len(model_names) > 1
    distinct_colors = ['brown', 'green', 'blue', 'cyan', 'yellow', 'pink', 'orange', '#dd37fa', '#ad5691', '#f46fa1', '#89cff0', '#93c83e'] \
                      if len(model_names) <= 12 else plt.cm.YlOrRd(np.linspace(1, 0, len(model_names) + 1))
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

            name_map = {str(m): str(m).upper() for m in sorted(df['model_name'].dropna().unique())}
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
                        ax.bar(x_positions, complete_counts.values, width=bar_width, alpha=0.4, color=color, edgecolor='black', linewidth=0.3, label=label)
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
    plt.savefig(f"{folder_to_save}descriptors_distribution.png", dpi=300, bbox_inches='tight', format='png')

    return
    
