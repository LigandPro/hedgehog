"""Descriptor computation functions."""

from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import QED, Crippen, Descriptors, Lipinski, rdMolDescriptors
from medchem.rules._utils import n_fused_aromatic_rings

from hedgehog.configs.logger import logger
from hedgehog.descriptors.io import order_identity_columns, process_path
from hedgehog.utils.mce18 import compute_mce18
from hedgehog.utils.parallel import parallel_map, resolve_n_jobs


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

    symbols = list({atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetSymbol()})
    has_formal_charge = any(atom.GetFormalCharge() != 0 for atom in mol.GetAtoms())
    is_neutral = not has_formal_charge
    # DEPRECATED: charged_mol has inverted semantics (True = neutral molecule).
    # Kept for backwards compatibility. Use is_neutral or has_formal_charge instead.
    charged_mol = is_neutral

    ring_info = mol.GetRingInfo()
    rings = [len(x) for x in ring_info.AtomRings()]

    n_rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol_n)
    n_rigid_bonds = mol_n.GetNumBonds() - n_rot_bonds
    n_heavy_atoms = rdMolDescriptors.CalcNumHeavyAtoms(mol_n)
    n_aromatic_atoms = sum(
        1 for a in mol_n.GetAtoms() if a.GetIsAromatic() and a.GetAtomicNum() > 1
    )
    mol_wt = Descriptors.ExactMolWt(mol_n)
    clogp = Crippen.MolLogP(mol_n)
    n_N_atoms = sum(1 for atom in mol_n.GetAtoms() if atom.GetAtomicNum() == 7)
    n_S_atoms = sum(1 for atom in mol_n.GetAtoms() if atom.GetAtomicNum() == 16)

    return {
        "model_name": model_name,
        "mol_idx": mol_idx,
        "chars": symbols,
        "n_atoms": mol.GetNumAtoms(),
        "n_heavy_atoms": n_heavy_atoms,
        "n_het_atoms": sum(
            1 for atom in mol_n.GetAtoms() if atom.GetAtomicNum() not in (1, 6)
        ),
        "n_N_atoms": n_N_atoms,
        "fN_atoms": n_N_atoms / n_heavy_atoms if n_heavy_atoms > 0 else 0,
        "fNS_atoms": (n_N_atoms + n_S_atoms) / n_heavy_atoms
        if n_heavy_atoms > 0
        else 0,
        "is_neutral": is_neutral,
        "has_formal_charge": has_formal_charge,
        "charged_mol": charged_mol,  # DEPRECATED: use is_neutral instead
        "molWt": mol_wt,
        "logP": Descriptors.MolLogP(mol_n),
        "clogP": clogp,
        "sw": 0.16
        - 0.63 * clogp
        - 0.0062 * mol_wt
        + 0.066 * n_rot_bonds
        - 0.74 * n_aromatic_atoms,
        "ring_size": rings,
        "n_rings": mol_n.GetRingInfo().NumRings(),
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
    }


def _compute_descriptors_for_row(args):
    smiles, model_name, mol_idx, remove_charges, remove_radicals, remove_stereo = args
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

    if remove_charges:
        try:
            from rdkit.Chem.MolStandardize import rdMolStandardize

            uncharger = rdMolStandardize.Uncharger()
            mol_n = uncharger.uncharge(mol_n)
        except Exception:
            # Best-effort: if uncharger isn't available, fall back to charge check only.
            pass

        still_charged = any(atom.GetFormalCharge() != 0 for atom in mol_n.GetAtoms())
        if still_charged:
            return None, (smiles, model_name, mol_idx)

    if remove_stereo:
        try:
            Chem.RemoveStereochemistry(mol_n)
        except Exception:
            pass

    if remove_charges or remove_radicals or remove_stereo:
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
    """Compute 22 physicochemical descriptors for each molecule.

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

    remove_charges = bool(preprocess_cfg.get("remove_charges", False))
    remove_radicals = bool(preprocess_cfg.get("remove_radicals", False))
    remove_stereo = bool(preprocess_cfg.get("remove_stereochemistry", False))
    do_preprocess = remove_charges or remove_radicals or remove_stereo

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
    items = [
        (s, m, idx, remove_charges, remove_radicals, remove_stereo)
        for (s, m, idx) in items
    ]

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
