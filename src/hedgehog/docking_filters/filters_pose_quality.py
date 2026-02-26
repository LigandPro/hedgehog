"""Pose quality filters (PoseCheck legacy and posecheck-fast)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from rdkit import Chem

from hedgehog.configs.logger import logger
from hedgehog.docking_filters.types import (
    _POSEBUSTERS_ERROR,
    _aggregate_filter_results,
    _error_result,
    _fail_result,
)
from hedgehog.utils.parallel import parallel_map

# ---------------------------------------------------------------------------
# Pose quality filter (PoseCheck legacy)
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# PoseBusters / posecheck-fast filter
# ---------------------------------------------------------------------------


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
            return _fail_result("no conformer", **_POSEBUSTERS_ERROR)

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
        return _error_result(
            no_clashes=no_clashes,
            no_volume_clash=no_volume_clash,
            not_too_far_away=not_too_far,
            no_internal_clash=no_internal_clash,
            passed=passed,
        )
    except Exception as e:
        return _fail_result(str(e), **_POSEBUSTERS_ERROR)


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

    return _aggregate_filter_results(
        raw_results,
        col_mapping={
            "no_clashes": "no_clashes",
            "no_volume_clash": "no_volume_clash",
            "not_too_far_away": "not_too_far_away",
            "no_internal_clash": "no_internal_clash",
        },
        pass_col="pass_pose_quality",
        filter_name="posecheck-fast filter",
    )


def apply_posecheck_fast_filter(
    mols: list[Chem.Mol],
    protein_pdb: str | Path,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Backward-compatible alias for posecheck-fast backend."""
    return apply_posebusters_fast_filter(mols, protein_pdb, config)
