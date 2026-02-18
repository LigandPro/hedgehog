"""Utility functions for docking pose filtering."""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from hedgehog.configs.logger import logger
from hedgehog.utils.parallel import parallel_map, resolve_n_jobs

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
    logger.info("Loaded %d molecules from %s", len(mols), sdf_path)
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
        "Running pose quality filter (max_clashes=%d, max_strain=%.1f)",
        max_clashes,
        max_strain,
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
        "Pose quality filter: %d/%d passed",
        int(df["pass_pose_quality"].sum()),
        len(df),
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

    logger.info("Running interaction filter (min_hbonds=%d)", min_hbonds)

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
    n_jobs = resolve_n_jobs(config)
    fp.run_from_iterable(plf_mols, protein, n_jobs=n_jobs)

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
    logger.info(
        "Interaction filter: %d/%d passed",
        int(df["pass_interactions"].sum()),
        len(df),
    )
    return df


def _check_shepherd_score_single(mol: Chem.Mol) -> dict[str, Any]:
    """Compute shape Tanimoto score for a single molecule.

    Args:
        mol: RDKit molecule.

    Returns:
        Dict with keys: score, passed, error.
    """
    import torch
    from shepherd_score.score.gaussian_overlap import shape_tanimoto

    ref_positions_list = _SHEPHERD_REF_POSITIONS_LIST
    alpha = _SHEPHERD_ALPHA
    min_score = _SHEPHERD_MIN_SCORE
    try:
        if ref_positions_list is None:
            raise RuntimeError("Shepherd-Score worker configuration was not initialized")
        ref_positions = torch.tensor(ref_positions_list, dtype=torch.float32)

        mol_h = Chem.AddHs(mol, addCoords=True)
        mol_conf = mol_h.GetConformer()
        mol_positions = torch.tensor(mol_conf.GetPositions(), dtype=torch.float32)

        score = shape_tanimoto(ref_positions, mol_positions, alpha=alpha).item()
        passed = score >= min_score
        return {"score": score, "passed": passed, "error": None}
    except Exception as e:
        return {"score": 0.0, "passed": False, "error": str(e)}


_SHEPHERD_REF_POSITIONS_LIST: list[list[float]] | None = None
_SHEPHERD_ALPHA: float = 0.81
_SHEPHERD_MIN_SCORE: float = 0.5


def _init_shepherd_score_worker(
    ref_positions_list: list[list[float]],
    alpha: float,
    min_score: float,
) -> None:
    global _SHEPHERD_REF_POSITIONS_LIST, _SHEPHERD_ALPHA, _SHEPHERD_MIN_SCORE
    _SHEPHERD_REF_POSITIONS_LIST = ref_positions_list
    _SHEPHERD_ALPHA = float(alpha)
    _SHEPHERD_MIN_SCORE = float(min_score)


def apply_shepherd_score_filter(
    mols: list[Chem.Mol],
    reference_mol: Chem.Mol,
    config: dict[str, Any],
    progress_cb=None,
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
    min_shape_score = config.get("min_shape_score", 0.5)
    alpha = config.get("alpha", 0.81)

    logger.info("Running Shepherd-Score filter (min_shape_score=%.2f)", min_shape_score)

    # Get reference coordinates as plain Python list for picklability
    ref_h = Chem.AddHs(reference_mol, addCoords=True)
    ref_conf = ref_h.GetConformer()
    ref_positions_list = ref_conf.GetPositions().tolist()

    n_jobs = resolve_n_jobs(config)
    raw_results = parallel_map(
        _check_shepherd_score_single,
        mols,
        n_jobs,
        progress=progress_cb,
        initializer=_init_shepherd_score_worker,
        initargs=(ref_positions_list, alpha, min_shape_score),
    )

    results = []
    for i, res in enumerate(raw_results):
        if res["error"]:
            logger.warning("Shepherd-Score failed for mol %d: %s", i, res["error"])
        results.append(
            {
                "mol_idx": i,
                "shape_score": res["score"],
                "pass_shepherd_score": res["passed"],
            }
        )

    df = pd.DataFrame(results)
    logger.info(
        "Shepherd-Score filter: %d/%d passed",
        int(df["pass_shepherd_score"].sum()),
        len(df),
    )
    return df


def _check_posebusters_fast_single(args: tuple) -> dict[str, Any]:
    """Run posecheck-fast checks for a single molecule.

    Args:
        args: Tuple of (mol, protein_coords, protein_atom_names, config).

    Returns:
        Dict with keys: no_clashes, no_volume_clash, not_too_far_away,
        no_internal_clash, passed, error.
    """
    from posecheck_fast import check_intermolecular_distance

    mol, protein_coords, protein_atom_names, config = args
    try:
        if mol is None or mol.GetNumConformers() == 0:
            return {
                "no_clashes": False,
                "no_volume_clash": False,
                "not_too_far_away": False,
                "no_internal_clash": False,
                "passed": False,
                "error": "no conformer",
            }

        conf = mol.GetConformer()
        pos_pred = np.array([conf.GetPositions()])  # (1, n_atoms, 3)
        atom_names_pred = np.array([atom.GetSymbol() for atom in mol.GetAtoms()])

        # check_intermolecular_distance runs all checks in one call:
        # clashes, volume overlap, distance to protein, and internal geometry.
        # Returns {"results": {"no_clashes": [...], ...}} with one bool per pose.
        raw = check_intermolecular_distance(
            mol_orig=mol,
            pos_pred=pos_pred,
            pos_cond=protein_coords,
            atom_names_pred=atom_names_pred,
            atom_names_cond=protein_atom_names,
            clash_cutoff=config.get("clash_cutoff", 0.75),
            clash_cutoff_volume=config.get("volume_clash_cutoff", 0.075),
            max_distance=config.get("max_distance", 5.0),
        )
        r = raw["results"]

        no_clashes = r["no_clashes"][0]
        no_volume_clash = r["no_volume_clash"][0]
        not_too_far = r["not_too_far_away"][0]
        no_internal_clash = r["no_internal_clash"][0]

        passed = no_clashes and no_volume_clash and not_too_far and no_internal_clash
        return {
            "no_clashes": no_clashes,
            "no_volume_clash": no_volume_clash,
            "not_too_far_away": not_too_far,
            "no_internal_clash": no_internal_clash,
            "passed": passed,
            "error": None,
        }
    except Exception as e:
        return {
            "no_clashes": False,
            "no_volume_clash": False,
            "not_too_far_away": False,
            "no_internal_clash": False,
            "passed": False,
            "error": str(e),
        }


def apply_posebusters_fast_filter(
    mols: list[Chem.Mol],
    protein_pdb: str | Path,
    config: dict[str, Any],
    progress_cb=None,
) -> pd.DataFrame:
    """Apply pose quality filter using posecheck-fast (~100x faster than PoseCheck).

    Checks for steric clashes (VDW), volume overlap (ShapeTverskyIndex),
    distance to protein, and internal geometry (bond lengths/angles).

    Args:
        mols: List of RDKit molecules with 3D coordinates.
        protein_pdb: Path to protein PDB file.
        config: Filter configuration dict with keys:
            clash_cutoff (float): Relative VDW distance threshold (default 0.75).
            volume_clash_cutoff (float): ShapeTverskyIndex threshold (default 0.075).
            max_distance (float): Max min-distance to protein in Angstrom (default 5.0).

    Returns:
        DataFrame with columns: mol_idx, no_clashes, no_volume_clash,
        not_too_far_away, no_internal_clash, pass_pose_quality.
    """
    logger.info(
        "Running posecheck-fast pose quality filter "
        "(clash_cutoff=%.2f, volume_cutoff=%.3f, max_dist=%.1f)",
        config.get("clash_cutoff", 0.75),
        config.get("volume_clash_cutoff", 0.075),
        config.get("max_distance", 5.0),
    )

    # Load protein once and extract coordinates + atom names.
    # Explicit RemoveHs is needed because MolFromPDBFile(removeHs=True) often
    # retains explicit H from prepared PDB files, causing dimension mismatches
    # inside posecheck-fast's ShapeTverskyIndex call.
    protein_mol = Chem.MolFromPDBFile(str(protein_pdb), removeHs=True, sanitize=False)
    if protein_mol is not None:
        protein_mol = Chem.RemoveHs(protein_mol, sanitize=False)
    if protein_mol is None:
        logger.error("Failed to load protein PDB: %s", protein_pdb)
        return pd.DataFrame(
            {
                "mol_idx": range(len(mols)),
                "pass_pose_quality": False,
            }
        )

    protein_conf = protein_mol.GetConformer()
    protein_coords = np.array(protein_conf.GetPositions())  # (n_atoms, 3)
    protein_atom_names = np.array([atom.GetSymbol() for atom in protein_mol.GetAtoms()])

    items = [(mol, protein_coords, protein_atom_names, config) for mol in mols]
    # Force sequential for posecheck-fast (uses torch internally, forking unsafe)
    raw_results = parallel_map(
        _check_posebusters_fast_single, items, n_jobs=1, progress=progress_cb
    )

    results = []
    for i, res in enumerate(raw_results):
        if res["error"]:
            logger.warning("posecheck-fast failed for mol %d: %s", i, res["error"])
        results.append(
            {
                "mol_idx": i,
                "no_clashes": res["no_clashes"],
                "no_volume_clash": res["no_volume_clash"],
                "not_too_far_away": res["not_too_far_away"],
                "no_internal_clash": res["no_internal_clash"],
                "pass_pose_quality": res["passed"],
            }
        )

    df = pd.DataFrame(results)
    logger.info(
        "posecheck-fast filter: %d/%d passed",
        int(df["pass_pose_quality"].sum()),
        len(df),
    )
    return df


def apply_posecheck_fast_filter(
    mols: list[Chem.Mol],
    protein_pdb: str | Path,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Backward-compatible alias for posecheck-fast backend."""
    return apply_posebusters_fast_filter(mols, protein_pdb, config)


def _check_symmetry_rmsd_single(args: tuple) -> dict[str, Any]:
    """Compute symmetry-corrected conformer RMSD for a single molecule.

    Uses RDKit's GetBestRMS which considers molecular symmetry (via substructure
    match) when computing RMSD — same correctness as spyrmsd but ~15x faster.

    Args:
        args: Tuple of (mol, num_conformers, max_rmsd, method_name, random_seed).

    Returns:
        Dict with keys: min_rmsd, n_conformers_generated, passed, error.
    """
    from rdkit.Chem import rdMolAlign

    (
        mol,
        num_conformers,
        max_rmsd,
        method_name,
        random_seed,
        include_hydrogens,
        max_matches,
        early_stop_on_pass,
    ) = args
    try:
        if include_hydrogens:
            reference_mol = Chem.AddHs(mol, addCoords=True)
        else:
            reference_mol = Chem.RemoveHs(Chem.Mol(mol), sanitize=False)

        # Select embedding method
        if method_name == "ETKDGv3":
            params = AllChem.ETKDGv3()
        elif method_name == "ETKDGv2":
            params = AllChem.ETKDGv2()
        else:
            params = AllChem.ETKDG()
        params.randomSeed = random_seed

        mol_confs = Chem.Mol(reference_mol)
        mol_confs.RemoveAllConformers()

        AllChem.EmbedMultipleConfs(mol_confs, numConfs=num_conformers, params=params)
        n_generated = mol_confs.GetNumConformers()

        if n_generated == 0:
            return {
                "min_rmsd": float("inf"),
                "n_conformers_generated": 0,
                "passed": False,
                "error": "no conformers generated",
            }

        # GetBestRMS: Kabsch alignment + symmetry-aware atom permutation.
        # NOTE: hydrogens can cause combinatorial explosion in symmetry matching.
        rmsd_values = []
        for conf_id in range(n_generated):
            rmsd = rdMolAlign.GetBestRMS(
                mol_confs,
                reference_mol,
                prbId=conf_id,
                refId=0,
                maxMatches=max_matches,
            )
            rmsd_values.append(rmsd)
            if early_stop_on_pass and rmsd <= max_rmsd:
                break

        min_rmsd = min(rmsd_values)
        return {
            "min_rmsd": min_rmsd,
            "n_conformers_generated": n_generated,
            "passed": min_rmsd <= max_rmsd,
            "error": None,
        }
    except Exception as e:
        return {
            "min_rmsd": float("inf"),
            "n_conformers_generated": 0,
            "passed": False,
            "error": str(e),
        }


def apply_symmetry_rmsd_filter(
    mols: list[Chem.Mol],
    config: dict[str, Any],
    progress_cb=None,
) -> pd.DataFrame:
    """Apply conformer deviation filter using symmetry-corrected RMSD.

    Uses symmetry-aware RDKit GetBestRMS instead of
    naive AllChem.AlignMol. This gives chemically correct RMSD for
    symmetric molecules (e.g., benzene, biphenyl).

    Args:
        mols: List of RDKit molecules with 3D coordinates.
        config: Filter configuration dict (same keys as conformer_deviation).

    Returns:
        DataFrame with columns: mol_idx, min_conformer_rmsd,
        n_conformers_generated, pass_conformer_deviation.
    """
    num_conformers = config.get("num_conformers", 50)
    max_rmsd = config.get("max_rmsd_to_conformer", 3.0)
    random_seed = config.get("random_seed", 42)
    method = config.get("conformer_method", "ETKDGv3")
    include_hydrogens = bool(config.get("include_hydrogens", False))
    max_matches = int(config.get("max_matches", 10000))
    early_stop_on_pass = bool(config.get("early_stop_on_pass", True))

    logger.info(
        "Running symmetry-RMSD conformer filter (max_rmsd=%.1f, n_confs=%d)",
        max_rmsd,
        num_conformers,
    )

    n_jobs = resolve_n_jobs(config)
    items = [
        (
            mol,
            num_conformers,
            max_rmsd,
            method,
            random_seed,
            include_hydrogens,
            max_matches,
            early_stop_on_pass,
        )
        for mol in mols
    ]
    raw_results = parallel_map(
        _check_symmetry_rmsd_single, items, n_jobs, progress=progress_cb
    )

    results = []
    for i, res in enumerate(raw_results):
        if res["error"]:
            logger.warning("Symmetry-RMSD failed for mol %d: %s", i, res["error"])
        results.append(
            {
                "mol_idx": i,
                "min_conformer_rmsd": res["min_rmsd"],
                "n_conformers_generated": res["n_conformers_generated"],
                "pass_conformer_deviation": res["passed"],
            }
        )

    df = pd.DataFrame(results)
    logger.info(
        "Symmetry-RMSD filter: %d/%d passed",
        int(df["pass_conformer_deviation"].sum()),
        len(df),
    )
    return df


def _check_conformer_deviation_single(args: tuple) -> dict[str, Any]:
    """Check conformer RMSD deviation for a single molecule.

    Args:
        args: Tuple of (mol, num_conformers, max_rmsd, method_name, random_seed).

    Returns:
        Dict with keys: min_rmsd, n_conformers_generated, passed, error.
    """
    mol, num_conformers, max_rmsd, method_name, random_seed = args
    try:
        mol_h = Chem.AddHs(mol, addCoords=True)

        # Select embedding method
        if method_name == "ETKDGv3":
            params = AllChem.ETKDGv3()
        elif method_name == "ETKDGv2":
            params = AllChem.ETKDGv2()
        else:
            params = AllChem.ETKDG()
        params.randomSeed = random_seed

        mol_confs = Chem.Mol(mol_h)
        mol_confs.RemoveAllConformers()

        AllChem.EmbedMultipleConfs(mol_confs, numConfs=num_conformers, params=params)
        n_generated = mol_confs.GetNumConformers()

        if n_generated == 0:
            return {
                "min_rmsd": float("inf"),
                "n_conformers_generated": 0,
                "passed": False,
                "error": "no conformers generated",
            }

        rmsd_values = []
        for conf_id in range(n_generated):
            conf_mol = Chem.Mol(mol_confs)
            conf_mol.RemoveAllConformers()
            conf_mol.AddConformer(mol_confs.GetConformer(conf_id), assignId=True)
            rmsd = AllChem.AlignMol(conf_mol, mol_h)
            rmsd_values.append(rmsd)

        min_rmsd = min(rmsd_values)
        return {
            "min_rmsd": min_rmsd,
            "n_conformers_generated": n_generated,
            "passed": min_rmsd <= max_rmsd,
            "error": None,
        }
    except Exception as e:
        return {
            "min_rmsd": float("inf"),
            "n_conformers_generated": 0,
            "passed": False,
            "error": str(e),
        }


def apply_conformer_deviation_filter(
    mols: list[Chem.Mol],
    config: dict[str, Any],
    progress_cb=None,
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
        "Running conformer deviation filter (max_rmsd=%.1f, n_confs=%d)",
        max_rmsd,
        num_conformers,
    )

    n_jobs = resolve_n_jobs(config)
    items = [(mol, num_conformers, max_rmsd, method, random_seed) for mol in mols]
    raw_results = parallel_map(
        _check_conformer_deviation_single, items, n_jobs, progress=progress_cb
    )

    results = []
    for i, res in enumerate(raw_results):
        if res["error"]:
            logger.warning("Conformer deviation failed for mol %d: %s", i, res["error"])
        results.append(
            {
                "mol_idx": i,
                "min_conformer_rmsd": res["min_rmsd"],
                "n_conformers_generated": res["n_conformers_generated"],
                "pass_conformer_deviation": res["passed"],
            }
        )

    df = pd.DataFrame(results)
    logger.info(
        "Conformer deviation filter: %d/%d passed",
        int(df["pass_conformer_deviation"].sum()),
        len(df),
    )
    return df
