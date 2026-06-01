"""Descriptor computation functions."""

from collections import deque
from pathlib import Path

import pandas as pd
from medchem.rules._utils import (
    fraction_atom_in_scaff,
    has_spider_chains,
    n_fused_aromatic_rings,
)
from rdkit import Chem, RDConfig
from rdkit.Chem import (
    QED,
    ChemicalFeatures,
    Descriptors,
    Lipinski,
    rdMolDescriptors,
)

from hedgehog.configs.logger import logger
from hedgehog.descriptors.constants import TYPE_ALIAS_COLUMNS
from hedgehog.descriptors.io import order_identity_columns, process_path
from hedgehog.utils.mce18 import compute_mce18
from hedgehog.utils.parallel import parallel_map, resolve_n_jobs

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
    """Compute 28 physicochemical descriptors for each molecule.

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
