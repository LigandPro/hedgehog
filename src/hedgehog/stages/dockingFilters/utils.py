"""Utility functions for docking pose filtering."""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from hedgehog.configs.logger import logger

# Suppress common warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)


def _project_root() -> Path:
    # src/hedgehog/stages/dockingFilters/utils.py -> project root
    return Path(__file__).resolve().parent.parent.parent.parent.parent


def _resolve_existing_path(base_folder: Path, path: str | Path) -> Path:
    """Resolve a path (possibly relative) to an existing absolute path."""
    p = Path(path)
    if p.is_absolute():
        return p

    candidates = [
        (base_folder / p).resolve(),
        (_project_root() / p).resolve(),
        (Path.cwd() / p).resolve(),
    ]
    for c in candidates:
        if c.exists():
            return c
    return (base_folder / p).resolve()


def _autobox_from_reference_ligand(
    base_folder: Path, docking_cfg: dict[str, Any]
) -> tuple[tuple[float, float, float], tuple[float, float, float]] | None:
    """Build a Vina/GNINA-style search box from autobox_ligand + autobox_add.

    Returns:
        (center_xyz, size_xyz) or None if insufficient config.
    """
    tool_cfg = docking_cfg.get("gnina_config") or docking_cfg.get("smina_config") or {}
    autobox_ligand = tool_cfg.get("autobox_ligand") or docking_cfg.get("autobox_ligand")
    if not autobox_ligand:
        return None

    ref_path = _resolve_existing_path(base_folder, autobox_ligand)
    if not ref_path.exists():
        return None

    suppl = Chem.SDMolSupplier(str(ref_path))
    ref_mol = next((m for m in suppl if m is not None), None)
    if ref_mol is None or ref_mol.GetNumConformers() == 0:
        return None

    conf = ref_mol.GetConformer()
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    for i in range(ref_mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        xs.append(float(p.x))
        ys.append(float(p.y))
        zs.append(float(p.z))

    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    min_z, max_z = min(zs), max(zs)

    add = float(tool_cfg.get("autobox_add", 0.0) or 0.0)
    min_x -= add
    max_x += add
    min_y -= add
    max_y += add
    min_z -= add
    max_z += add

    center = ((min_x + max_x) / 2, (min_y + max_y) / 2, (min_z + max_z) / 2)
    size = (max_x - min_x, max_y - min_y, max_z - min_z)
    return center, size


def apply_search_box_filter(
    mols: list[Chem.Mol],
    base_folder: Path,
    docking_cfg: dict[str, Any],
    config: dict[str, Any],
) -> pd.DataFrame:
    """Fast filter: ensure docked poses are inside the configured search box.

    Supports explicit center/size, or autobox_ligand-based boxes.

    Returns:
        DataFrame with columns: mol_idx, frac_atoms_outside_box, pass_search_box
    """
    enabled = bool(config.get("enabled", True))
    if not enabled:
        return pd.DataFrame(
            {
                "mol_idx": range(len(mols)),
                "frac_atoms_outside_box": 0.0,
                "pass_search_box": True,
            }
        )

    # Prefer explicit center/size if present, otherwise fall back to autobox reference.
    center = docking_cfg.get("center")
    size = docking_cfg.get("size")
    if not (
        isinstance(center, (list, tuple))
        and len(center) >= 3
        and isinstance(size, (list, tuple))
        and len(size) >= 3
    ):
        box = _autobox_from_reference_ligand(base_folder, docking_cfg)
        if box is None:
            logger.warning(
                "Search-box filter enabled but no box could be resolved; passing all poses."
            )
            return pd.DataFrame(
                {
                    "mol_idx": range(len(mols)),
                    "frac_atoms_outside_box": 0.0,
                    "pass_search_box": True,
                }
            )
        center, size = box

    cx, cy, cz = float(center[0]), float(center[1]), float(center[2])
    sx, sy, sz = float(size[0]), float(size[1]), float(size[2])
    hx, hy, hz = sx / 2.0, sy / 2.0, sz / 2.0

    max_outside = float(config.get("max_outside_fraction", 0.0) or 0.0)

    frac_outside: list[float] = []
    passed: list[bool] = []

    for mol in mols:
        if mol is None or mol.GetNumConformers() == 0:
            frac_outside.append(1.0)
            passed.append(False)
            continue
        conf = mol.GetConformer()
        n_atoms = mol.GetNumAtoms()
        outside = 0
        for i in range(n_atoms):
            p = conf.GetAtomPosition(i)
            x, y, z = float(p.x), float(p.y), float(p.z)
            if (
                x < cx - hx
                or x > cx + hx
                or y < cy - hy
                or y > cy + hy
                or z < cz - hz
                or z > cz + hz
            ):
                outside += 1
        frac = outside / n_atoms if n_atoms else 1.0
        frac_outside.append(frac)
        passed.append(frac <= max_outside)

    df = pd.DataFrame(
        {
            "mol_idx": range(len(mols)),
            "frac_atoms_outside_box": frac_outside,
            "pass_search_box": passed,
        }
    )
    logger.info(
        "Search-box filter: %d/%d passed (max_outside_fraction=%.3f)",
        int(df["pass_search_box"].sum()),
        len(df),
        max_outside,
    )
    return df


def load_molecules_from_sdf(sdf_path: str | Path) -> list[Chem.Mol]:
    """Load molecules from SDF file."""
    sdf_path = Path(sdf_path)
    if not sdf_path.exists():
        raise FileNotFoundError(f"SDF file not found: {sdf_path}")

    suppl = Chem.SDMolSupplier(str(sdf_path))
    mols = [mol for mol in suppl if mol is not None]
    logger.info(f"Loaded {len(mols)} molecules from {sdf_path}")
    return mols


def apply_pose_quality_filter(
    mols: list[Chem.Mol],
    protein_pdb: str | Path,
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Apply pose quality filter using PoseCheck.

    Checks for steric clashes and ligand strain energy.

    Args:
        mols: List of RDKit molecules with 3D coordinates
        protein_pdb: Path to protein PDB file
        config: Filter configuration dict

    Returns:
        DataFrame with columns: mol_idx, clashes, strain_energy, pass
    """
    from posecheck import PoseCheck

    max_clashes = config.get("max_clashes", 2)
    max_strain = config.get("max_strain_energy", 10.0)

    logger.info(
        f"Running pose quality filter (max_clashes={max_clashes}, max_strain={max_strain})"
    )

    pc = PoseCheck()
    pc.load_protein_from_pdb(str(protein_pdb))
    pc.ligands = mols

    # Calculate metrics
    clashes = pc.calculate_clashes()
    strain_energies = pc.calculate_strain_energy()

    # Build results DataFrame
    df = pd.DataFrame(
        {
            "mol_idx": range(len(mols)),
            "clashes": clashes,
            "strain_energy": strain_energies,
        }
    )
    df["pass_pose_quality"] = (df["clashes"] <= max_clashes) & (
        df["strain_energy"] <= max_strain
    )
    logger.info(
        f"Pose quality filter: {df['pass_pose_quality'].sum()}/{len(df)} passed"
    )
    return df


def apply_interaction_filter(
    mols: list[Chem.Mol],
    protein_pdb: str | Path,
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Apply interaction filter using ProLIF.

    Checks for required protein-ligand interactions.

    Args:
        mols: List of RDKit molecules with 3D coordinates
        protein_pdb: Path to protein PDB file
        config: Filter configuration dict

    Returns:
        DataFrame with columns: mol_idx, n_hbonds, interactions, pass
    """
    import MDAnalysis as mda
    import prolif as plf

    min_hbonds = config.get("min_hbonds", 0)
    required_residues = config.get("required_residues", [])
    forbidden_residues = config.get("forbidden_residues", [])
    interaction_types = config.get(
        "interaction_types", ["HBDonor", "HBAcceptor", "Hydrophobic", "VdWContact"]
    )

    logger.info(f"Running interaction filter (min_hbonds={min_hbonds})")

    # Load protein
    u = mda.Universe(str(protein_pdb))
    protein = plf.Molecule.from_mda(u)

    # Convert RDKit mols to ProLIF molecules
    plf_mols = []
    for mol in mols:
        mol_h = Chem.AddHs(mol, addCoords=True)
        plf_mols.append(plf.Molecule.from_rdkit(mol_h))

    # Calculate interactions (exclude FaceToFace due to edge case errors)
    safe_interactions = [it for it in interaction_types if it != "FaceToFace"]
    fp = plf.Fingerprint(interactions=safe_interactions)
    fp.run_from_iterable(plf_mols, protein, n_jobs=1)

    df_interactions = fp.to_dataframe(drop_empty=False)

    # Build results DataFrame
    results = []
    for i in range(len(mols)):
        if i < len(df_interactions):
            row = df_interactions.iloc[i]
            # Count H-bonds
            hbond_cols = [
                c for c in row.index if "HBDonor" in str(c) or "HBAcceptor" in str(c)
            ]
            n_hbonds = sum(row[hbond_cols].values) if hbond_cols else 0

            # Check required residues
            has_required = True
            for res in required_residues:
                res_cols = [c for c in row.index if res in str(c)]
                if res_cols and not any(row[res_cols].values):
                    has_required = False
                    break

            # Check forbidden residues
            has_forbidden = False
            for res in forbidden_residues:
                res_cols = [c for c in row.index if res in str(c)]
                if res_cols and any(row[res_cols].values):
                    has_forbidden = True
                    break

            passed = n_hbonds >= min_hbonds and has_required and not has_forbidden
            interactions_str = ",".join([str(c) for c in row.index if row[c]])
        else:
            n_hbonds = 0
            passed = min_hbonds == 0 and not required_residues
            interactions_str = ""

        results.append(
            {
                "mol_idx": i,
                "n_hbonds": n_hbonds,
                "interactions": interactions_str,
                "pass_interactions": passed,
            }
        )

    df = pd.DataFrame(results)
    logger.info(f"Interaction filter: {df['pass_interactions'].sum()}/{len(df)} passed")
    return df


def apply_shepherd_score_filter(
    mols: list[Chem.Mol],
    reference_mol: Chem.Mol,
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Apply Shepherd-Score filter for 3D shape similarity.

    Args:
        mols: List of RDKit molecules with 3D coordinates
        reference_mol: Reference molecule for comparison
        config: Filter configuration dict

    Returns:
        DataFrame with columns: mol_idx, shape_score, pass
    """
    import torch
    from shepherd_score.score.gaussian_overlap import shape_tanimoto

    min_shape_score = config.get("min_shape_score", 0.5)
    alpha = config.get("alpha", 0.81)

    logger.info(f"Running Shepherd-Score filter (min_shape_score={min_shape_score})")

    # Get reference coordinates
    ref_h = Chem.AddHs(reference_mol, addCoords=True)
    ref_conf = ref_h.GetConformer()
    ref_positions = torch.tensor(ref_conf.GetPositions(), dtype=torch.float32)

    results = []
    for i, mol in enumerate(mols):
        try:
            mol_h = Chem.AddHs(mol, addCoords=True)
            mol_conf = mol_h.GetConformer()
            mol_positions = torch.tensor(mol_conf.GetPositions(), dtype=torch.float32)

            # Calculate shape Tanimoto
            score = shape_tanimoto(ref_positions, mol_positions, alpha=alpha).item()
            passed = score >= min_shape_score
        except Exception as e:
            logger.warning(f"Shepherd-Score failed for mol {i}: {e}")
            score = 0.0
            passed = False

        results.append(
            {
                "mol_idx": i,
                "shape_score": score,
                "pass_shepherd_score": passed,
            }
        )

    df = pd.DataFrame(results)
    logger.info(
        f"Shepherd-Score filter: {df['pass_shepherd_score'].sum()}/{len(df)} passed"
    )
    return df


def apply_conformer_deviation_filter(
    mols: list[Chem.Mol],
    config: dict[str, Any],
) -> pd.DataFrame:
    """
    Apply conformer deviation filter.

    Checks if docked pose is geometrically plausible by comparing
    to generated conformers.

    Args:
        mols: List of RDKit molecules with 3D coordinates
        config: Filter configuration dict

    Returns:
        DataFrame with columns: mol_idx, min_rmsd, n_conformers, pass
    """
    num_conformers = config.get("num_conformers", 50)
    max_rmsd = config.get("max_rmsd_to_conformer", 3.0)
    random_seed = config.get("random_seed", 42)
    method = config.get("conformer_method", "ETKDGv3")

    logger.info(
        f"Running conformer deviation filter (max_rmsd={max_rmsd}, n_confs={num_conformers})"
    )

    # Select embedding method
    if method == "ETKDGv3":
        params = AllChem.ETKDGv3()
    elif method == "ETKDGv2":
        params = AllChem.ETKDGv2()
    else:
        params = AllChem.ETKDG()
    params.randomSeed = random_seed

    results = []
    for i, mol in enumerate(mols):
        try:
            # Add hydrogens
            mol_h = Chem.AddHs(mol, addCoords=True)

            # Create copy for conformer generation
            mol_confs = Chem.Mol(mol_h)
            mol_confs.RemoveAllConformers()

            # Generate conformers
            AllChem.EmbedMultipleConfs(
                mol_confs, numConfs=num_conformers, params=params
            )
            n_generated = mol_confs.GetNumConformers()

            if n_generated == 0:
                logger.warning(f"No conformers generated for mol {i}")
                min_rmsd = float("inf")
                passed = False
            else:
                # Calculate RMSD to each conformer
                rmsd_values = []
                for conf_id in range(n_generated):
                    # Create molecule with single conformer
                    conf_mol = Chem.Mol(mol_confs)
                    conf_mol.RemoveAllConformers()
                    conf_mol.AddConformer(
                        mol_confs.GetConformer(conf_id), assignId=True
                    )

                    # Align and get RMSD
                    rmsd = AllChem.AlignMol(conf_mol, mol_h)
                    rmsd_values.append(rmsd)

                min_rmsd = min(rmsd_values)
                passed = min_rmsd <= max_rmsd

        except Exception as e:
            logger.warning(f"Conformer deviation failed for mol {i}: {e}")
            min_rmsd = float("inf")
            n_generated = 0
            passed = False

        results.append(
            {
                "mol_idx": i,
                "min_conformer_rmsd": min_rmsd,
                "n_conformers_generated": n_generated,
                "pass_conformer_deviation": passed,
            }
        )

    df = pd.DataFrame(results)
    logger.info(
        f"Conformer deviation filter: {df['pass_conformer_deviation'].sum()}/{len(df)} passed"
    )
    return df
