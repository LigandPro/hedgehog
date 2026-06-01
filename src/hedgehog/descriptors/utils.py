import ast
import json
import math
from collections import deque
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from medchem.rules._utils import (
    fraction_atom_in_scaff,
    has_spider_chains,
    n_fused_aromatic_rings,
)
from rdkit import Chem, RDConfig, RDLogger, rdBase
from rdkit.Chem import (
    QED,
    ChemicalFeatures,
    Descriptors,
    Lipinski,
    rdMolDescriptors,
)

from hedgehog._constants import CFG_DESCRIPTORS
from hedgehog.configs.logger import load_config, logger
from hedgehog.descriptors.constants import (
    _DESCRIPTOR_KEY_MAP,
    STRUCTURAL_ELEMENT_LIMIT_MAP,
    TYPE_ALIAS_COLUMNS,
)
from hedgehog.utils.mce18 import compute_mce18
from hedgehog.utils.parallel import parallel_map, resolve_n_jobs
from hedgehog.utils.paths import process_path

# Disable RDKit warnings
RDLogger.DisableLog("rdApp.*")
rdBase.DisableLog("rdApp.*")

try:
    _FEATURE_FACTORY = ChemicalFeatures.BuildFeatureFactory(
        str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef")
    )
except Exception:
    _FEATURE_FACTORY = None


def _collect_donor_acceptor_atom_sets(mol):
    """Collect donor/acceptor atom ids from RDKit feature definitions."""
    if _FEATURE_FACTORY is None:
        return set(), set()

    donors = set()
    acceptors = set()
    for feature in _FEATURE_FACTORY.GetFeaturesForMol(mol):
        family = feature.GetFamily()
        if family == "Donor":
            donors.update(feature.GetAtomIds())
        elif family == "Acceptor":
            acceptors.update(feature.GetAtomIds())
    return donors, acceptors


def _count_small_rings_3_4(rings):
    """Count 3/4-membered rings."""
    return sum(1 for ring_size in rings if ring_size in (3, 4))


def _compute_max_acyclic_chain_length(mol):
    """Compute maximum atom count along a non-ring chain."""
    adjacency = {}
    for bond in mol.GetBonds():
        if bond.IsInRing():
            continue
        begin = bond.GetBeginAtom()
        end = bond.GetEndAtom()
        if begin.GetAtomicNum() == 1 or end.GetAtomicNum() == 1:
            continue
        a_idx = begin.GetIdx()
        b_idx = end.GetIdx()
        adjacency.setdefault(a_idx, set()).add(b_idx)
        adjacency.setdefault(b_idx, set()).add(a_idx)

    if not adjacency:
        return 0

    def _farthest(start, allowed):
        visited = {start}
        queue = deque([(start, 0)])
        farthest_node = start
        farthest_dist = 0
        while queue:
            node, dist = queue.popleft()
            if dist > farthest_dist:
                farthest_node = node
                farthest_dist = dist
            for neighbor in adjacency.get(node, ()):
                if neighbor not in allowed or neighbor in visited:
                    continue
                visited.add(neighbor)
                queue.append((neighbor, dist + 1))
        return farthest_node, farthest_dist

    remaining = set(adjacency)
    best_diameter_edges = 0

    while remaining:
        start = next(iter(remaining))
        component = {start}
        queue = deque([start])
        while queue:
            node = queue.popleft()
            for neighbor in adjacency.get(node, ()):
                if neighbor in component:
                    continue
                component.add(neighbor)
                queue.append(neighbor)
        remaining.difference_update(component)

        endpoint, _ = _farthest(start, component)
        _, diameter_edges = _farthest(endpoint, component)
        best_diameter_edges = max(best_diameter_edges, diameter_edges)

    return best_diameter_edges + 1


def _compute_type_alias_counts(mol):
    """Compute configured atom type alias counts."""
    donors, acceptors = _collect_donor_acceptor_atom_sets(mol)
    counts = {alias: 0 for alias in TYPE_ALIAS_COLUMNS}

    so2_sulfur_ids = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16:
            continue
        n_double_bonded_oxygens = 0
        for bond in atom.GetBonds():
            neighbor = bond.GetOtherAtom(atom)
            if (
                neighbor.GetAtomicNum() == 8
                and bond.GetBondType() == Chem.BondType.DOUBLE
            ):
                n_double_bonded_oxygens += 1
        if n_double_bonded_oxygens >= 2:
            so2_sulfur_ids.add(atom.GetIdx())

    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        atom_idx = atom.GetIdx()
        in_ring = atom.IsInRing()
        is_aromatic = atom.GetIsAromatic()
        hyb = atom.GetHybridization()
        formal_charge = atom.GetFormalCharge()
        total_h = atom.GetTotalNumHs()

        if atomic_num == 8:
            if hyb == Chem.HybridizationType.SP2 and atom_idx in acceptors:
                counts[".=O"] += 1
            if (
                hyb == Chem.HybridizationType.SP3
                and atom_idx in acceptors
                and total_h == 0
            ):
                counts["O_a"] += 1
            if hyb == Chem.HybridizationType.SP3 and atom_idx in donors and total_h > 0:
                counts["O_d"] += 1

        elif atomic_num == 6:
            if is_aromatic:
                counts["Car"] += 1
            if in_ring and not is_aromatic and hyb == Chem.HybridizationType.SP2:
                counts["C2r"] += 1
            if in_ring and hyb == Chem.HybridizationType.SP3:
                counts["C3r"] += 1
            if not in_ring and not is_aromatic and hyb == Chem.HybridizationType.SP2:
                counts["Cs2"] += 1
            if not in_ring and hyb == Chem.HybridizationType.SP3:
                counts["Cs3"] += 1
            if hyb == Chem.HybridizationType.SP:
                counts["Csp"] += 1

        elif atomic_num == 7:
            if formal_charge == 0 and atom_idx in acceptors:
                counts["Nac"] += 1
            if formal_charge > 0 and atom_idx in donors and total_h > 0:
                counts["Nd+"] += 1
            if formal_charge == 0 and atom_idx in donors and total_h > 0:
                counts["Nd0"] += 1

        elif atomic_num == 16:
            if atom_idx in so2_sulfur_ids:
                counts["SO2"] += 1
            elif atom.GetTotalValence() == 2:
                counts["Sul"] += 1

        elif atomic_num in (9, 17, 35, 53):
            counts["Hal"] += 1

    return counts


def order_identity_columns(df):
    """Reorder dataframe columns with identity columns first."""
    id_cols = ["smiles", "model_name", "mol_idx"]
    existing_id_cols = [c for c in id_cols if c in df.columns]
    ordered = existing_id_cols + [c for c in df.columns if c not in id_cols]
    return df[ordered]


def drop_false_rows(df, borders):
    """
    Filter rows that passed all descriptor filters.

    Args:
        df: DataFrame with '_pass' columns
        borders: Configuration dict with filter settings

    Returns:
        pd.DataFrame: Filtered dataframe with only passed molecules
    """
    passed_cols = [col for col in df.columns if col.endswith("_pass") or col == "pass"]
    mask = df[passed_cols].all(axis=1)
    return df[mask].copy()


def _parse_literal_list(value, label):
    """Parse list-like values from strings safely using ast.literal_eval.

    This function deserializes lists that were saved as strings in CSV files
    (e.g., ['C', 'N', 'O'] for atom symbols or [5, 6] for ring sizes).

    Security: ast.literal_eval is safe for untrusted input. Unlike Python's
    built-in eval function, literal_eval only parses literal data structures
    (strings, numbers, tuples, lists, dicts, sets, booleans, None) and
    rejects any other Python code, preventing code injection attacks.

    Args:
        value: Value to parse (string representation of a list or list itself)
        label: Label for error messages (e.g., "ring sizes", "type limits")

    Returns:
        Parsed list, or empty list if parsing fails
    """
    if isinstance(value, list):
        return value
    if isinstance(value, str):
        try:
            parsed = ast.literal_eval(value)
            if isinstance(parsed, list):
                return parsed
        except Exception as e:
            logger.error("Error parsing %s: %s, %s", label, value, e)
            return []
    logger.error("Error parsing %s: %s, unsupported type", label, value)
    return []


def _parse_ring_size_column(series):
    """Parse ring size lists from series for plotting."""
    parsed = []
    for val in series.dropna():
        sizes = _parse_literal_list(val, "ring sizes")
        try:
            parsed.extend([float(size) for size in sizes])
        except Exception as e:
            logger.error("Error parsing ring sizes: %s, %s", val, e)
    return parsed


def _compute_single_molecule_descriptors(mol_n, model_name, mol_idx):
    """Compute all descriptors for a single molecule.

    Args:
        mol_n: RDKit molecule object (without hydrogens)
        model_name: Model identifier string
        mol_idx: Molecule index

    Returns:
        dict: Dictionary of computed descriptors
    """
    mol = Chem.AddHs(mol_n)

    ring_info = mol.GetRingInfo()
    rings = [len(x) for x in ring_info.AtomRings()]

    n_rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol_n)
    n_rigid_bonds = mol_n.GetNumBonds() - n_rot_bonds
    n_heavy_atoms = rdMolDescriptors.CalcNumHeavyAtoms(mol_n)
    n_aromatic_atoms = sum(
        1 for a in mol_n.GetAtoms() if a.GetIsAromatic() and a.GetAtomicNum() > 1
    )
    mol_wt = Descriptors.ExactMolWt(mol_n)
    log_p = Descriptors.MolLogP(mol_n)
    n_N_atoms = sum(1 for atom in mol_n.GetAtoms() if atom.GetAtomicNum() == 7)
    n_O_atoms = sum(1 for atom in mol_n.GetAtoms() if atom.GetAtomicNum() == 8)
    n_S_atoms = sum(1 for atom in mol_n.GetAtoms() if atom.GetAtomicNum() == 16)
    n_NO_atoms = n_N_atoms + n_O_atoms
    n_small_rings_3_4 = _count_small_rings_3_4(rings)
    max_acyclic_chain_length = _compute_max_acyclic_chain_length(mol_n)
    type_alias_counts = _compute_type_alias_counts(mol_n)

    return {
        "model_name": model_name,
        "mol_idx": mol_idx,
        "n_atoms": mol.GetNumAtoms(),
        "n_heavy_atoms": n_heavy_atoms,
        "n_het_atoms": sum(
            1 for atom in mol_n.GetAtoms() if atom.GetAtomicNum() not in (1, 6)
        ),
        "n_N_atoms": n_N_atoms,
        "n_O_atoms": n_O_atoms,
        "n_S_atoms": n_S_atoms,
        "n_NO_atoms": n_NO_atoms,
        "fN_atoms": n_N_atoms / n_heavy_atoms if n_heavy_atoms > 0 else 0,
        "fNS_atoms": (n_N_atoms + n_S_atoms) / n_heavy_atoms
        if n_heavy_atoms > 0
        else 0,
        "molWt": mol_wt,
        "logP": log_p,
        "sw": 0.16
        - 0.63 * log_p
        - 0.0062 * mol_wt
        + 0.066 * n_rot_bonds
        - 0.74 * n_aromatic_atoms,
        "ring_size": rings,
        "n_rings": mol_n.GetRingInfo().NumRings(),
        "n_small_rings_3_4": n_small_rings_3_4,
        "max_acyclic_chain_length": max_acyclic_chain_length,
        "has_spider_side_chains": int(has_spider_chains(mol_n)),
        "fraction_ring_system": fraction_atom_in_scaff(mol_n),
        "n_aroma_rings": rdMolDescriptors.CalcNumAromaticRings(mol_n),
        "n_fused_aromatic_rings": n_fused_aromatic_rings(mol_n),
        "n_rigid_bonds": n_rigid_bonds,
        "n_rot_bonds": n_rot_bonds,
        "hbd": Lipinski.NumHDonors(mol_n),
        "hba": Lipinski.NumHAcceptors(mol_n),
        "fsp3": rdMolDescriptors.CalcFractionCSP3(mol_n),
        "mce18": compute_mce18(mol_n),
        "tpsa": rdMolDescriptors.CalcTPSA(mol_n),
        "qed": QED.qed(mol_n),
        **type_alias_counts,
    }


def _compute_descriptors_for_row(args):
    smiles, model_name, mol_idx, remove_radicals, remove_stereo = args
    if not isinstance(smiles, str):
        if pd.isna(smiles):
            return None, ("", model_name, mol_idx)
        smiles = str(smiles)

    smiles = smiles.strip()
    if smiles.lower() in {"", "nan", "none"}:
        return None, (smiles, model_name, mol_idx)

    mol_n = Chem.MolFromSmiles(smiles)
    if not mol_n:
        return None, (smiles, model_name, mol_idx)

    # Descriptor-level preprocessing is kept as an optional fallback.
    # Prefer MolPrep stage for standardization to avoid double-processing.
    if remove_radicals:
        has_radical = any(
            atom.GetNumRadicalElectrons() > 0 for atom in mol_n.GetAtoms()
        )
        if has_radical:
            return None, (smiles, model_name, mol_idx)

    if remove_stereo:
        try:
            Chem.RemoveStereochemistry(mol_n)
        except Exception:
            pass

    if remove_radicals or remove_stereo:
        try:
            smiles = Chem.MolToSmiles(mol_n, canonical=True, isomericSmiles=False)
        except Exception:
            return None, (smiles, model_name, mol_idx)

    row_metrics = _compute_single_molecule_descriptors(mol_n, model_name, mol_idx)
    row_metrics["smiles"] = smiles
    return row_metrics, None


def compute_metrics(
    df,
    save_path,
    config=None,
    config_descriptors=None,
    reporter=None,
    progress_stage_total: int | None = None,
    progress_completed_base: int = 0,
    progress_completed_span: int | None = None,
):
    """Compute 28 base physicochemical descriptors for each molecule.

    model_name and mol_idx are already in df from sampled_molecules.csv.

    Args:
        df: DataFrame with molecules (must have 'smiles', 'model_name', 'mol_idx')
        save_path: Output folder path
        config: Configuration dictionary

    Returns:
        pd.DataFrame: Dataframe with computed descriptors per molecule
    """
    if df is None or len(df) == 0:
        logger.warning(
            "Empty DataFrame provided to compute_metrics. Returning empty DataFrame."
        )
        return pd.DataFrame()

    metrics = []
    skipped_molecules = []

    preprocess_cfg = {}
    if isinstance(config_descriptors, dict):
        preprocess_cfg = config_descriptors.get("preprocess", {}) or {}

    remove_radicals = bool(preprocess_cfg.get("remove_radicals", False))
    remove_stereo = bool(preprocess_cfg.get("remove_stereochemistry", False))
    do_preprocess = remove_radicals or remove_stereo

    n_jobs = resolve_n_jobs(config_descriptors, config)
    logger.info("Descriptors workers: %d", n_jobs)
    items = list(
        zip(
            df["smiles"].tolist(),
            df["model_name"].tolist(),
            df["mol_idx"].tolist(),
            strict=False,
        )
    )
    items = [(s, m, idx, remove_radicals, remove_stereo) for (s, m, idx) in items]

    progress_cb = None
    if reporter is not None:

        def _progress_cb(done: int, total: int) -> None:
            mapped_done = done
            mapped_total = total
            if (
                progress_stage_total is not None
                and progress_completed_span is not None
                and total > 0
            ):
                ratio = min(1.0, max(0.0, done / total))
                mapped_done = progress_completed_base + int(
                    round(ratio * progress_completed_span)
                )
                mapped_total = progress_stage_total
            reporter.progress(
                mapped_done, mapped_total, message="Computing descriptors"
            )

        progress_cb = _progress_cb

    results = parallel_map(
        _compute_descriptors_for_row, items, n_jobs, progress=progress_cb
    )
    for row_metrics, skipped in results:
        if row_metrics is not None:
            metrics.append(row_metrics)
        elif skipped is not None:
            skipped_molecules.append(skipped)

    save_path = Path(process_path(save_path))
    if skipped_molecules:
        if do_preprocess:
            logger.warning(
                "Skipped %d molecules during preprocessing or SMILES parsing",
                len(skipped_molecules),
            )
        else:
            logger.warning(
                "Skipped %d molecules that failed to parse", len(skipped_molecules)
            )
        skipped_df = pd.DataFrame(
            {
                "smiles": [s for s, _, _ in skipped_molecules],
                "model_name": [m for _, m, _ in skipped_molecules],
            }
        )
        if any(idx is not None for _, _, idx in skipped_molecules):
            skipped_df["mol_idx"] = [idx for _, _, idx in skipped_molecules]
        skipped_df.to_csv(save_path / "skipped_molecules.csv", index=False)

    metrics_df = pd.DataFrame(metrics)
    if metrics_df.shape[1] == 0:
        metrics_df = pd.DataFrame(columns=["smiles", "model_name", "mol_idx"])
    metrics_df = order_identity_columns(metrics_df)
    metrics_df.to_csv(save_path / "descriptors_all.csv", index=False)
    if (
        reporter is not None
        and progress_stage_total is not None
        and progress_completed_span is not None
    ):
        reporter.progress(
            progress_completed_base + progress_completed_span,
            progress_stage_total,
            message="Descriptors computed",
        )
    return metrics_df


def _get_border_values(col, borders):
    """Get min and max border values for a column from borders config.

    Uses explicit key mapping to prevent case-insensitive confusion between
    similar descriptor names with different casing.

    Args:
        col: Column name
        borders: Dictionary with border configurations

    Returns:
        tuple: (min_border, max_border)
    """
    # Try exact match first (preserves case distinction)
    min_key = f"{col}_min"
    max_key = f"{col}_max"
    if min_key in borders or max_key in borders:
        return borders.get(min_key), borders.get(max_key)

    # Try canonical mapping for case-insensitive lookup
    col_lower = col.lower()
    canonical = _DESCRIPTOR_KEY_MAP.get(col_lower, col)
    min_key = f"{canonical}_min"
    max_key = f"{canonical}_max"
    if min_key in borders or max_key in borders:
        return borders.get(min_key), borders.get(max_key)

    # Fallback: search all keys (legacy behavior)
    min_border = None
    max_border = None
    for key, value in borders.items():
        if key.endswith("_min"):
            base = key[: -len("_min")]
            if base.lower() == col_lower:
                min_border = value
        elif key.endswith("_max"):
            base = key[: -len("_max")]
            if base.lower() == col_lower:
                max_border = value
    return min_border, max_border


def _apply_column_filter(df, col, borders):
    """Apply filter to a single column based on borders configuration.

    Args:
        df: DataFrame with data
        col: Column name to filter
        borders: Dictionary with border configurations

    Returns:
        pd.Series: Boolean series indicating pass/fail for each row
    """
    min_border, max_border = _get_border_values(col, borders)

    # Normalize common YAML/string sentinels (e.g., "inf") for numeric comparisons.
    def _normalize_border_value(val):
        if val is None:
            return None
        if isinstance(val, str) and val.lower() == "inf":
            return math.inf
        return val

    min_border = _normalize_border_value(min_border)
    max_border = _normalize_border_value(max_border)

    if col == "ring_size":
        # Missing bounds mean "no bound" on that side.
        if min_border is None:
            min_border = -math.inf
        if max_border is None:
            max_border = math.inf
        return df[col].apply(
            lambda x: all(
                min_border <= float(ring_size) <= max_border
                for ring_size in _parse_literal_list(x, "ring sizes")
            )
        )

    # Generic numeric range filter.
    #
    # Allow one-sided bounds:
    # - only min -> value >= min
    # - only max -> value <= max
    if min_border is None and max_border is None:
        return pd.Series(True, index=df.index)
    if min_border is None:
        return df[col] <= max_border
    if max_border is None or max_border == math.inf:
        return df[col] >= min_border
    return (df[col] >= min_border) & (df[col] <= max_border)


def _order_descriptor_columns(df):
    """Order DataFrame columns with identity cols, then descriptors with their pass flags.

    Args:
        df: DataFrame to reorder

    Returns:
        list: Ordered column names
    """
    id_cols = ["smiles", "model_name", "mol_idx"]
    descriptor_cols = [
        col for col in df.columns if col not in id_cols and not col.endswith("_pass")
    ]
    pass_cols = [col for col in df.columns if col.endswith("_pass")]

    ordered_cols = id_cols.copy()
    for desc_col in sorted(descriptor_cols):
        if desc_col in df.columns:
            ordered_cols.append(desc_col)
            pass_col = f"{desc_col}_pass"
            if pass_col in pass_cols:
                ordered_cols.append(pass_col)

    for pass_col in sorted(pass_cols):
        if pass_col not in ordered_cols:
            ordered_cols.append(pass_col)

    return [col for col in ordered_cols if col in df.columns]


def _merge_pass_flags(df, flags_path):
    """Merge pass flags from flags file into dataframe.

    Args:
        df: DataFrame to merge flags into
        flags_path: Path to flags CSV file

    Returns:
        pd.DataFrame: DataFrame with merged pass flags
    """
    if not Path(flags_path).exists():
        return df

    flags_df = pd.read_csv(flags_path)
    merge_cols = ["smiles", "model_name"]
    pass_cols = [
        col for col in flags_df.columns if col.endswith("_pass") or col == "pass"
    ]

    if not pass_cols:
        return df

    df = df.merge(
        flags_df[merge_cols + pass_cols],
        on=merge_cols,
        how="left",
        suffixes=("", "_flags"),
    )
    for col in pass_cols:
        if f"{col}_flags" in df.columns:
            df[col] = df[f"{col}_flags"].fillna(df.get(col, False))
            df = df.drop(columns=[f"{col}_flags"])

    return df


def _save_failed_molecules(fail_filters, folder_to_save, flags_path):
    """Save failed molecules to CSV files.

    Args:
        fail_filters: DataFrame with failed molecules
        folder_to_save: Output folder path (Path object)
        flags_path: Path to pass flags CSV
    """
    fail_filters = _merge_pass_flags(fail_filters, flags_path)
    fail_filters = order_identity_columns(fail_filters)

    ordered_cols = _order_descriptor_columns(fail_filters)
    fail_filters[ordered_cols].to_csv(
        folder_to_save / "descriptors_failed.csv", index=False
    )

    id_cols = ["smiles", "model_name", "mol_idx"]
    fail_filters[id_cols].to_csv(folder_to_save / "failed_molecules.csv", index=False)


def _extract_borders_and_constraints(borders, structural_constraints):
    """Normalize descriptors config into plain borders + structural constraints."""
    if not isinstance(borders, dict):
        if isinstance(structural_constraints, dict):
            return {}, structural_constraints
        return {}, {}

    normalized_borders = borders
    constraints = structural_constraints

    if "borders" in borders and isinstance(borders.get("borders"), dict):
        normalized_borders = borders.get("borders", {}) or {}
        if constraints is None:
            constraints = borders.get("structural_constraints")
    elif constraints is None:
        constraints = borders.get("structural_constraints")

    normalized_borders = {
        key: value
        for key, value in normalized_borders.items()
        if key != "structural_constraints"
    }
    if not isinstance(constraints, dict):
        constraints = {}
    return normalized_borders, constraints


def _build_structural_constraint_borders(structural_constraints):
    """Convert structural constraints block to regular *_max border entries."""
    if not structural_constraints.get("enabled", False):
        return {}

    translated = {}

    type_limits = structural_constraints.get("type_limits") or {}
    if isinstance(type_limits, dict):
        for alias, max_value in type_limits.items():
            if max_value is None:
                continue
            alias_key = str(alias)
            canonical_alias = _DESCRIPTOR_KEY_MAP.get(alias_key.lower(), alias_key)
            translated[f"{canonical_alias}_max"] = max_value

    element_limits = structural_constraints.get("element_limits") or {}
    if isinstance(element_limits, dict):
        for key, max_value in element_limits.items():
            if max_value is None:
                continue
            key_str = str(key)
            element_col = STRUCTURAL_ELEMENT_LIMIT_MAP.get(
                key_str, STRUCTURAL_ELEMENT_LIMIT_MAP.get(key_str.lower())
            )
            if element_col is not None:
                translated[f"{element_col}_max"] = max_value

    direct_limits = {
        "max_n_or_o_atoms": "n_NO_atoms",
        "max_small_rings_3_4": "n_small_rings_3_4",
        "max_acyclic_chain_length": "max_acyclic_chain_length",
    }
    for config_key, column in direct_limits.items():
        max_value = structural_constraints.get(config_key)
        if max_value is not None:
            translated[f"{column}_max"] = max_value

    return translated


def filter_molecules(df, borders, folder_to_save, structural_constraints=None):
    """Filter molecules based on descriptor thresholds.

    Args:
        df: DataFrame with computed descriptors
        borders: Dictionary with min/max thresholds for each descriptor
        folder_to_save: Output folder path (should already include 'Descriptors' subfolder)
    """
    folder_to_save = Path(process_path(folder_to_save))
    id_cols = ["smiles", "model_name", "mol_idx"]
    normalized_borders, constraints = _extract_borders_and_constraints(
        borders, structural_constraints
    )
    effective_borders = dict(normalized_borders)
    effective_borders.update(_build_structural_constraint_borders(constraints))

    logger.info("[#B29EEE]Applied Descriptor Filters:[/#B29EEE]")
    for line in json.dumps(effective_borders, indent=2, ensure_ascii=False).split("\n"):
        logger.info("  %s", line)

    # Build filtered data with pass flags
    filtered_data = {}
    for col in df.columns.tolist():
        if col in id_cols:
            filtered_data[col] = df[col]
            continue

        col_in_borders = any(
            v is not None for v in _get_border_values(col, effective_borders)
        )
        if col_in_borders:
            filtered_data[col] = df[col]
            filtered_data[f"{col}_pass"] = _apply_column_filter(
                df, col, effective_borders
            )

    filtered_data_df = pd.DataFrame(filtered_data)
    filtered_data_df = order_identity_columns(filtered_data_df)
    filtered_data_df.to_csv(
        folder_to_save / "pass_flags.csv", index_label="SMILES", index=False
    )

    pass_filters = drop_false_rows(filtered_data_df, effective_borders)

    # Save passed molecules
    if len(pass_filters) > 0:
        pass_filters = order_identity_columns(pass_filters)
        descriptor_cols = [
            col
            for col in pass_filters.columns
            if col not in id_cols and not col.endswith("_pass")
        ]
        ordered_cols = [
            col
            for col in id_cols + sorted(descriptor_cols)
            if col in pass_filters.columns
        ]
        pass_filters[ordered_cols].to_csv(
            folder_to_save / "descriptors_passed.csv", index=False
        )
        pass_filters[id_cols].to_csv(
            folder_to_save / "filtered_molecules.csv", index=False
        )
    else:
        logger.warning("No molecules pass Descriptors Filters")

    # Save failed molecules
    all_computed_path = folder_to_save / "descriptors_all.csv"
    if not all_computed_path.exists():
        # descriptors_all.csv is written into the sibling "metrics/" folder
        # by compute_metrics(). Fall back to that location for provenance.
        sibling_metrics = folder_to_save.parent / "metrics" / "descriptors_all.csv"
        if sibling_metrics.exists():
            all_computed_path = sibling_metrics
    flags_path = folder_to_save / "pass_flags.csv"

    if all_computed_path.exists():
        all_computed = pd.read_csv(all_computed_path)

        if len(pass_filters) > 0:
            merge_cols = ["smiles", "model_name"]
            merged = all_computed.merge(
                pass_filters[merge_cols], on=merge_cols, how="left", indicator=True
            )
            fail_filters = (
                merged[merged["_merge"] == "left_only"].drop(columns=["_merge"]).copy()
            )
        else:
            fail_filters = all_computed.copy()

        if len(fail_filters) > 0:
            _save_failed_molecules(fail_filters, folder_to_save, flags_path)

    # Re-order existing CSV files
    if (folder_to_save / "filtered_molecules.csv").exists():
        if all_computed_path.exists():
            per = pd.read_csv(all_computed_path)
            order_identity_columns(per).to_csv(all_computed_path, index=False)

        if flags_path.exists():
            flags = pd.read_csv(flags_path)
            order_identity_columns(flags).to_csv(flags_path, index=False)


def _get_model_colors(model_names):
    """Get color mapping for models.

    Args:
        model_names: List of model name strings

    Returns:
        dict: Mapping of model names to colors
    """
    distinct_colors = [
        "brown",
        "green",
        "blue",
        "cyan",
        "yellow",
        "pink",
        "orange",
        "#dd37fa",
        "#ad5691",
        "#f46fa1",
        "#89cff0",
        "#93c83e",
    ]
    if len(model_names) <= 12:
        colors = distinct_colors
    else:
        colors = plt.cm.YlOrRd(np.linspace(1, 0, len(model_names) + 1))
    return dict(zip(sorted(model_names), colors, strict=False))


def _get_column_values(model_df, col):
    """Extract values from a column, handling special columns.

    Args:
        model_df: DataFrame for a specific model
        col: Column name

    Returns:
        list: Extracted values
    """
    if col == "ring_size":
        return _parse_ring_size_column(model_df[col].dropna())
    return model_df[col].dropna().tolist()


def _filter_values_by_bounds(values, min_val, max_val):
    """Filter values within min/max bounds.

    Args:
        values: List of values
        min_val: Minimum value (or None)
        max_val: Maximum value (or None, or 'inf')

    Returns:
        list: Filtered values
    """
    result = values.copy()
    if min_val is not None:
        result = [v for v in result if v >= min_val]
    if max_val is not None and max_val != "inf":
        result = [v for v in result if v <= max_val]
    return result


def _plot_discrete_numeric(
    ax, values, col, offset, bar_width, color, label, max_val, model_index
):
    """Plot discrete numeric values as bar chart."""
    value_counts = pd.Series(values).value_counts().sort_index()

    if max_val is not None and max_val != "inf":
        extended_max = int(max_val) + 5
        full_range = list(range(0, extended_max + 1))
        complete_counts = pd.Series(0, index=full_range)
        for val in value_counts.index:
            if val in complete_counts.index:
                complete_counts[val] = value_counts[val]
    else:
        complete_counts = value_counts

    x_positions = [i + offset for i in range(len(complete_counts.index))]
    ax.bar(
        x_positions,
        complete_counts.values,
        width=bar_width,
        alpha=0.4,
        color=color,
        edgecolor="black",
        linewidth=0.3,
        label=label,
    )

    if model_index == 0:
        if col == "n_rigid_bonds":
            all_values = list(complete_counts.index)
            tick_values = [x for x in all_values if x % 5 == 0]
            if max(all_values) not in tick_values:
                tick_values.append(max(all_values))
            tick_positions = [
                all_values.index(val) for val in tick_values if val in all_values
            ]
            ax.set_xticks(tick_positions)
            ax.set_xticklabels([str(int(val)) for val in tick_values])
        else:
            ax.set_xticks(list(range(len(complete_counts.index))))
            ax.set_xticklabels([str(int(x)) for x in complete_counts.index])
        ax._discrete_tick_values = complete_counts.index


def _plot_continuous(ax, values, col, color, label):
    """Plot continuous distribution using KDE."""
    clip = (0, 1.0) if col in ["fsp3", "qed"] else (0, None)
    sns.kdeplot(
        values, label=label, fill=True, alpha=0.4, ax=ax, color=color, clip=clip
    )


def _get_boundary_position(ax, val, col, discrete_feats, fallback_offset=0):
    """Calculate boundary position for vertical lines and spans.

    Args:
        ax: Matplotlib axis
        val: Boundary value
        col: Column name
        discrete_feats: List of discrete feature names
        fallback_offset: Offset to add for continuous features

    Returns:
        float or None: Position for the boundary
    """
    if col not in discrete_feats:
        return val

    if not hasattr(ax, "_discrete_tick_values"):
        return val

    try:
        return list(ax._discrete_tick_values).index(val)
    except ValueError:
        return None


def _draw_boundary_lines(ax, col, min_val, max_val, discrete_feats):
    """Draw vertical lines at min/max boundaries."""
    # Small nudge offsets for continuous max boundaries that need
    # a tiny shift to avoid sitting exactly on the data edge.
    _MAX_NUDGE = {"fsp3": 0.01, "n_rigid_bonds": 0.0000001}

    boundaries = [
        (min_val, "red", "min", -0.5, 0),
        (max_val, "blue", "max", +0.5, None),
    ]
    for val, color, label_prefix, discrete_offset, _fallback_pos in boundaries:
        if val is None or val == "inf":
            continue

        pos = _get_boundary_position(ax, val, col, discrete_feats)
        if pos is not None:
            if col in discrete_feats:
                x = pos + discrete_offset
            elif label_prefix == "max" and col in _MAX_NUDGE:
                x = val + _MAX_NUDGE[col]
            else:
                x = val
        else:
            # Fallback position when value is not found among discrete ticks
            if label_prefix == "min":
                x = 0
            else:
                x = (
                    len(ax._discrete_tick_values) - 1
                    if hasattr(ax, "_discrete_tick_values")
                    else val
                )

        ax.axvline(
            x,
            color=color,
            linestyle="--",
            linewidth=1.5,
            label=f"{label_prefix}: {val}",
        )


def _draw_boundary_spans(ax, col, min_val, max_val, discrete_feats):
    """Draw shaded spans for excluded regions."""
    x_min, x_max = ax.get_xlim()

    # Draw min boundary span (shade area below min)
    if min_val is not None:
        pos = _get_boundary_position(ax, min_val, col, discrete_feats)
        if col in discrete_feats:
            if pos is not None:
                ax.axvspan(x_min, pos - 0.5, color="grey", alpha=0.2, zorder=0)
            else:
                ax.axvspan(x_min, 0, color="grey", alpha=0.2, zorder=0)
        else:
            ax.axvspan(x_min, min_val, color="grey", alpha=0.2, zorder=0)

    # Draw max boundary span (shade area above max)
    if max_val is not None and max_val != "inf":
        pos = _get_boundary_position(ax, max_val, col, discrete_feats)
        if col in discrete_feats:
            if pos is not None:
                ax.axvspan(pos + 0.5, x_max, color="grey", alpha=0.2, zorder=0)
            else:
                tick_len = (
                    len(ax._discrete_tick_values) - 0.5
                    if hasattr(ax, "_discrete_tick_values")
                    else x_max
                )
                ax.axvspan(tick_len, x_max, color="grey", alpha=0.2, zorder=0)
        else:
            ax.axvspan(max_val, x_max, color="grey", alpha=0.2, zorder=0)


def _style_axis(ax):
    """Apply consistent styling to axis."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.spines["left"].set_linewidth(0.5)
    ax.tick_params(
        axis="both",
        which="both",
        bottom=True,
        top=False,
        left=True,
        right=False,
        labelbottom=True,
        labeltop=False,
        labelleft=True,
        labelright=False,
        length=4,
        width=0.5,
        colors="black",
        labelsize=10,
    )


def _plot_single_column(
    ax, df, col, model_names, colors, borders, discrete_feats, renamer, is_multi
):
    """Plot distribution for a single column across all models.

    Args:
        ax: Matplotlib axis to plot on
        df: Full DataFrame
        col: Column name to plot
        model_names: List of model names
        colors: Color mapping for models
        borders: Border configurations
        discrete_feats: List of discrete feature names
        renamer: Column name mapping for display
        is_multi: Whether multiple models exist
    """
    min_val, max_val = _get_border_values(col, borders)
    minmax_str = f"min: {min_val}, max: {max_val}"
    bar_width = 0.1
    name_map = {
        str(m): str(m).upper() for m in sorted(df["model_name"].dropna().unique())
    }
    name_map_lc = {str(k).lower(): v for k, v in name_map.items()}

    for model_index, model in enumerate(model_names):
        model_df = df[df["model_name"].str.lower() == model] if is_multi else df
        values_before = _get_column_values(model_df, col)

        if not values_before:
            continue

        values_after = _filter_values_by_bounds(values_before, min_val, max_val)
        total = len(values_before)
        mols_passed = len(values_after) / total * 100 if total > 0 else 0

        label_name = name_map_lc.get(str(model).lower(), str(model).upper())
        label = f"{label_name}, pass: {mols_passed:.1f}%"
        color = colors[model]
        offset = model_index * bar_width

        if len(values_before) <= 1:
            ax.scatter(
                values_before,
                [0.01] * len(values_before),
                label=label,
                alpha=0.4,
                color=color,
            )
        elif col in discrete_feats:
            _plot_discrete_numeric(
                ax,
                values_before,
                col,
                offset,
                bar_width,
                color,
                label,
                max_val,
                model_index,
            )
        else:
            _plot_continuous(ax, values_before, col, color, label)

    # Set title and labels
    display_name = renamer.get(col, col)
    title = f"{display_name} ({minmax_str})"
    ax.set_title(title, fontsize=12)
    ax.set_xlabel(display_name, fontsize=10)

    # Set tick locators
    if col not in discrete_feats:
        if col == "fsp3":
            ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        elif col != "qed":
            ax.xaxis.set_major_locator(plt.MaxNLocator(integer=False))
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

    # Draw boundaries and styling
    _draw_boundary_lines(ax, col, min_val, max_val, discrete_feats)
    _draw_boundary_spans(ax, col, min_val, max_val, discrete_feats)

    # Add sorted legend
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        sorted_pairs = sorted(zip(labels, handles, strict=False), key=lambda t: t[0])
        sorted_labels, sorted_handles = zip(*sorted_pairs, strict=False)
        ax.legend(sorted_handles, sorted_labels, fontsize=8, loc="upper right")

    _style_axis(ax)


def draw_filtered_mols(df, folder_to_save, config, progress_cb=None):
    """Generate distribution plots for descriptor filters.

    Args:
        df: DataFrame with computed descriptors
        folder_to_save: Output folder path (should already include 'Descriptors' subfolder)
        config: Configuration dictionary
    """
    folder_to_save = Path(process_path(folder_to_save))

    descriptors_config = load_config(config[CFG_DESCRIPTORS])
    normalized_borders, constraints = _extract_borders_and_constraints(
        descriptors_config.get("borders", {}),
        descriptors_config.get("structural_constraints"),
    )
    borders = dict(normalized_borders)
    borders.update(_build_structural_constraint_borders(constraints))
    cols_to_plot = descriptors_config["filtered_cols_to_plot"]
    discrete_feats = descriptors_config["discrete_features_to_plot"]
    renamer = descriptors_config["renamer"]

    model_names = sorted(
        [m.lower() for m in df["model_name"].dropna().unique().tolist()]
    )
    is_multi = len(model_names) > 1
    colors = _get_model_colors(model_names)

    # Calculate grid dimensions
    n_cols = min(5, len(cols_to_plot))
    n_rows = (len(cols_to_plot) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(
        nrows=n_rows, ncols=n_cols, figsize=(n_rows * n_cols, n_rows * n_cols)
    )
    axes = axes.flatten()

    # Plot each column
    plotted_count = 0
    total_to_plot = max(len(cols_to_plot), 1)
    for i, col in enumerate(cols_to_plot):
        col_lower = col.lower()
        canonical = _DESCRIPTOR_KEY_MAP.get(col_lower, col).lower()
        relevant_keys = [
            k
            for k in borders
            if k.endswith(("_min", "_max"))
            and k.rsplit("_", 1)[0].lower() in (col_lower, canonical)
        ]
        if not relevant_keys:
            if progress_cb is not None:
                progress_cb(i + 1, total_to_plot)
            continue
        _plot_single_column(
            axes[i],
            df,
            col,
            model_names,
            colors,
            borders,
            discrete_feats,
            renamer,
            is_multi,
        )
        plotted_count += 1
        if progress_cb is not None:
            progress_cb(i + 1, total_to_plot)

    # Remove unused axes
    for j in range(len(cols_to_plot), len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.subplots_adjust(top=0.93, bottom=0.05, left=0.05, right=0.98)
    plt.savefig(
        folder_to_save / "descriptors_distribution.png",
        dpi=300,
        bbox_inches="tight",
        format="png",
    )
    if progress_cb is not None and plotted_count == 0:
        progress_cb(total_to_plot, total_to_plot)
