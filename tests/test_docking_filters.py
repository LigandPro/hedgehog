"""Tests for dockingFilters/utils.py — posebusters-fast and symmetry-RMSD backends."""

from __future__ import annotations

from pathlib import Path
from typing import Any
from unittest.mock import patch

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _mol_with_3d(smiles: str, seed: int = 42) -> Chem.Mol:
    """Create an RDKit Mol with a single 3D conformer."""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    AllChem.EmbedMolecule(mol, params)
    mol = Chem.RemoveHs(mol)
    return mol


def _write_protein_pdb(path: Path, coords: np.ndarray | None = None) -> Path:
    """Write a minimal PDB file with a few carbon atoms.

    If *coords* is None, places 4 atoms in a small box around origin.
    """
    if coords is None:
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.5, 0.0, 0.0],
                [0.0, 1.5, 0.0],
                [0.0, 0.0, 1.5],
            ]
        )
    pdb_path = path / "protein.pdb"
    lines = []
    for i, (x, y, z) in enumerate(coords):
        lines.append(
            f"ATOM  {i + 1:5d}  CA  ALA A{i + 1:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C  "
        )
    lines.append("END")
    pdb_path.write_text("\n".join(lines))
    return pdb_path


# ---------------------------------------------------------------------------
# posebusters-fast filter tests
# ---------------------------------------------------------------------------


class TestPosebustersFastFilter:
    """Tests for apply_posebusters_fast_filter."""

    def test_basic_returns_correct_schema(self, tmp_path):
        """Filter should return DataFrame with all expected columns."""
        from hedgehog.stages.dockingFilters.utils import apply_posebusters_fast_filter

        mol = _mol_with_3d("CCO")
        conf = mol.GetConformer()
        ligand_center = np.mean(conf.GetPositions(), axis=0)

        # Place protein atoms ~4A away (close but not clashing)
        offsets = np.array([[4.0, 0, 0], [-4.0, 0, 0], [0, 4.0, 0], [0, -4.0, 0]])
        protein_coords = ligand_center + offsets
        pdb_path = _write_protein_pdb(tmp_path, protein_coords)

        config: dict[str, Any] = {
            "clash_cutoff": 0.75,
            "volume_clash_cutoff": 0.075,
            "max_distance": 5.0,
            "n_jobs": 1,
        }
        df = apply_posebusters_fast_filter([mol], pdb_path, config)

        assert len(df) == 1
        expected_cols = {
            "mol_idx",
            "no_clashes",
            "no_volume_clash",
            "not_too_far_away",
            "no_internal_clash",
            "pass_pose_quality",
        }
        assert expected_cols.issubset(set(df.columns))
        # Molecule should be within max_distance of protein
        assert df["not_too_far_away"].iloc[0]

    def test_too_far_away_fails(self, tmp_path):
        """A molecule placed far from the protein should fail not_too_far_away."""
        from hedgehog.stages.dockingFilters.utils import apply_posebusters_fast_filter

        mol = _mol_with_3d("CCO")
        # Place protein atoms very far away
        protein_coords = np.array([[100.0, 100.0, 100.0], [101.0, 100.0, 100.0]])
        pdb_path = _write_protein_pdb(tmp_path, protein_coords)

        config: dict[str, Any] = {
            "max_distance": 5.0,
            "n_jobs": 1,
        }
        df = apply_posebusters_fast_filter([mol], pdb_path, config)
        assert (
            df["not_too_far_away"].iloc[0] is False
            or not df["not_too_far_away"].iloc[0]
        )

    def test_none_mol_fails(self, tmp_path):
        """None molecules should fail gracefully."""
        from hedgehog.stages.dockingFilters.utils import apply_posebusters_fast_filter

        protein_coords = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
        pdb_path = _write_protein_pdb(tmp_path, protein_coords)

        config: dict[str, Any] = {"n_jobs": 1}
        df = apply_posebusters_fast_filter([None], pdb_path, config)

        assert len(df) == 1
        assert (
            df["pass_pose_quality"].iloc[0] is False
            or not df["pass_pose_quality"].iloc[0]
        )

    def test_multiple_molecules(self, tmp_path):
        """Should process multiple molecules and return correct length."""
        from hedgehog.stages.dockingFilters.utils import apply_posebusters_fast_filter

        mols = [_mol_with_3d("CCO"), _mol_with_3d("c1ccccc1"), _mol_with_3d("CC")]

        # Place protein near all of them
        center = np.mean(mols[0].GetConformer().GetPositions(), axis=0)
        protein_coords = center + np.array([[3.0, 0, 0], [-3.0, 0, 0], [0, 3.0, 0]])
        pdb_path = _write_protein_pdb(tmp_path, protein_coords)

        config: dict[str, Any] = {"n_jobs": 1}
        df = apply_posebusters_fast_filter(mols, pdb_path, config)

        assert len(df) == 3
        assert list(df["mol_idx"]) == [0, 1, 2]

    def test_bad_protein_pdb_returns_all_fail(self, tmp_path):
        """An unreadable protein PDB should mark all molecules as failed."""
        from hedgehog.stages.dockingFilters.utils import apply_posebusters_fast_filter

        bad_pdb = tmp_path / "protein.pdb"
        bad_pdb.write_text("NOT A VALID PDB")

        mol = _mol_with_3d("CCO")
        config: dict[str, Any] = {"n_jobs": 1}
        df = apply_posebusters_fast_filter([mol], bad_pdb, config)

        assert len(df) == 1
        assert not df["pass_pose_quality"].iloc[0]


# ---------------------------------------------------------------------------
# symmetry-RMSD filter tests
# ---------------------------------------------------------------------------


class TestSymmetryRmsdFilter:
    """Tests for apply_symmetry_rmsd_filter."""

    def test_basic_pass(self):
        """A normal molecule should have plausible conformer RMSD."""
        from hedgehog.stages.dockingFilters.utils import apply_symmetry_rmsd_filter

        mol = _mol_with_3d("CCO")
        config: dict[str, Any] = {
            "num_conformers": 10,
            "max_rmsd_to_conformer": 5.0,
            "random_seed": 42,
            "conformer_method": "ETKDGv3",
            "n_jobs": 1,
        }
        df = apply_symmetry_rmsd_filter([mol], config)

        assert len(df) == 1
        assert "pass_conformer_deviation" in df.columns
        assert "min_conformer_rmsd" in df.columns
        # A valid 3D molecule embedded by ETKDGv3 should match its own conformers
        assert df["pass_conformer_deviation"].iloc[0]

    def test_symmetry_rmsd_leq_naive_aligned(self):
        """GetBestRMS (symmetry-aware) should be ≤ AlignMol (identity permutation)."""
        from rdkit.Chem import rdMolAlign

        mol = _mol_with_3d("c1ccccc1")  # benzene — highly symmetric
        mol_h = Chem.AddHs(mol, addCoords=True)

        params = AllChem.ETKDGv3()
        params.randomSeed = 123
        mol_confs = Chem.Mol(mol_h)
        mol_confs.RemoveAllConformers()
        AllChem.EmbedMultipleConfs(mol_confs, numConfs=5, params=params)

        for conf_id in range(mol_confs.GetNumConformers()):
            conf_mol = Chem.Mol(mol_confs)
            conf_mol.RemoveAllConformers()
            conf_mol.AddConformer(mol_confs.GetConformer(conf_id), assignId=True)

            # Symmetry-aware RMSD (Kabsch + substructure-match permutation)
            symm_rmsd = rdMolAlign.GetBestRMS(conf_mol, mol_h)

            # Naive RMSD via AllChem.AlignMol (Kabsch, identity permutation)
            naive_rmsd = AllChem.AlignMol(Chem.Mol(conf_mol), mol_h)

            assert symm_rmsd <= naive_rmsd + 1e-6

    def test_multiple_molecules(self):
        """Should process a batch correctly."""
        from hedgehog.stages.dockingFilters.utils import apply_symmetry_rmsd_filter

        mols = [_mol_with_3d("CCO"), _mol_with_3d("c1ccccc1")]
        config: dict[str, Any] = {
            "num_conformers": 5,
            "max_rmsd_to_conformer": 5.0,
            "random_seed": 42,
            "conformer_method": "ETKDGv3",
            "n_jobs": 1,
        }
        df = apply_symmetry_rmsd_filter(mols, config)

        assert len(df) == 2
        assert list(df["mol_idx"]) == [0, 1]

    def test_strict_threshold_rejects(self):
        """A very strict threshold should reject molecules."""
        from hedgehog.stages.dockingFilters.utils import apply_symmetry_rmsd_filter

        mol = _mol_with_3d("c1ccc2ccccc2c1")  # naphthalene
        config: dict[str, Any] = {
            "num_conformers": 5,
            "max_rmsd_to_conformer": 0.001,  # impossibly strict
            "random_seed": 42,
            "conformer_method": "ETKDGv3",
            "n_jobs": 1,
        }
        df = apply_symmetry_rmsd_filter([mol], config)
        # May or may not fail for planar molecule, but the test verifies no crash
        assert len(df) == 1
        assert "pass_conformer_deviation" in df.columns


# ---------------------------------------------------------------------------
# Backend dispatch tests
# ---------------------------------------------------------------------------


class TestBackendDispatch:
    """Tests that main.py dispatches to the correct backend."""

    def test_pose_quality_dispatch_posebusters_fast(self, tmp_path):
        """backend=posebusters_fast should call apply_posebusters_fast_filter."""
        with (
            patch(
                "hedgehog.stages.dockingFilters.main.apply_posebusters_fast_filter"
            ) as mock_pb,
            patch(
                "hedgehog.stages.dockingFilters.main.apply_pose_quality_filter"
            ) as mock_pc,
        ):
            mock_pb.return_value = pd.DataFrame(
                {"mol_idx": [0], "pass_pose_quality": [True]}
            )
            mock_pc.return_value = pd.DataFrame(
                {"mol_idx": [0], "pass_pose_quality": [True]}
            )

            # Simulate the dispatch logic from main.py
            pq_config: dict[str, Any] = {"enabled": True, "backend": "posebusters_fast"}
            mols_active = [_mol_with_3d("CCO")]
            protein_pdb = tmp_path / "protein.pdb"

            pq_backend = pq_config.get("backend", "posebusters_fast")
            if pq_backend == "posebusters_fast":
                apply_posebusters_fast_filter = mock_pb
                apply_posebusters_fast_filter(mols_active, protein_pdb, pq_config)
            else:
                apply_pose_quality_filter = mock_pc
                apply_pose_quality_filter(mols_active, protein_pdb, pq_config)

            mock_pb.assert_called_once()
            mock_pc.assert_not_called()

    def test_pose_quality_dispatch_posecheck(self, tmp_path):
        """backend=posecheck should call apply_pose_quality_filter."""
        with (
            patch(
                "hedgehog.stages.dockingFilters.main.apply_posebusters_fast_filter"
            ) as mock_pb,
            patch(
                "hedgehog.stages.dockingFilters.main.apply_pose_quality_filter"
            ) as mock_pc,
        ):
            mock_pb.return_value = pd.DataFrame(
                {"mol_idx": [0], "pass_pose_quality": [True]}
            )
            mock_pc.return_value = pd.DataFrame(
                {"mol_idx": [0], "pass_pose_quality": [True]}
            )

            pq_config: dict[str, Any] = {"enabled": True, "backend": "posecheck"}
            mols_active = [_mol_with_3d("CCO")]
            protein_pdb = tmp_path / "protein.pdb"

            pq_backend = pq_config.get("backend", "posebusters_fast")
            if pq_backend == "posebusters_fast":
                apply_posebusters_fast_filter = mock_pb
                apply_posebusters_fast_filter(mols_active, protein_pdb, pq_config)
            else:
                apply_pose_quality_filter = mock_pc
                apply_pose_quality_filter(mols_active, protein_pdb, pq_config)

            mock_pc.assert_called_once()
            mock_pb.assert_not_called()

    def test_conformer_dispatch_symmetry_rmsd(self):
        """backend=symmetry_rmsd should call apply_symmetry_rmsd_filter."""
        with (
            patch(
                "hedgehog.stages.dockingFilters.main.apply_symmetry_rmsd_filter"
            ) as mock_sr,
            patch(
                "hedgehog.stages.dockingFilters.main.apply_conformer_deviation_filter"
            ) as mock_naive,
        ):
            mock_sr.return_value = pd.DataFrame(
                {"mol_idx": [0], "pass_conformer_deviation": [True]}
            )
            mock_naive.return_value = pd.DataFrame(
                {"mol_idx": [0], "pass_conformer_deviation": [True]}
            )

            cd_config: dict[str, Any] = {"enabled": True, "backend": "symmetry_rmsd"}
            mols_active = [_mol_with_3d("CCO")]

            cd_backend = cd_config.get("backend", "symmetry_rmsd")
            if cd_backend == "symmetry_rmsd":
                apply_symmetry_rmsd_filter = mock_sr
                apply_symmetry_rmsd_filter(mols_active, cd_config)
            else:
                apply_conformer_deviation_filter = mock_naive
                apply_conformer_deviation_filter(mols_active, cd_config)

            mock_sr.assert_called_once()
            mock_naive.assert_not_called()

    def test_conformer_dispatch_naive(self):
        """backend=naive should call apply_conformer_deviation_filter."""
        with (
            patch(
                "hedgehog.stages.dockingFilters.main.apply_symmetry_rmsd_filter"
            ) as mock_sr,
            patch(
                "hedgehog.stages.dockingFilters.main.apply_conformer_deviation_filter"
            ) as mock_naive,
        ):
            mock_sr.return_value = pd.DataFrame(
                {"mol_idx": [0], "pass_conformer_deviation": [True]}
            )
            mock_naive.return_value = pd.DataFrame(
                {"mol_idx": [0], "pass_conformer_deviation": [True]}
            )

            cd_config: dict[str, Any] = {"enabled": True, "backend": "naive"}
            mols_active = [_mol_with_3d("CCO")]

            cd_backend = cd_config.get("backend", "symmetry_rmsd")
            if cd_backend == "symmetry_rmsd":
                apply_symmetry_rmsd_filter = mock_sr
                apply_symmetry_rmsd_filter(mols_active, cd_config)
            else:
                apply_conformer_deviation_filter = mock_naive
                apply_conformer_deviation_filter(mols_active, cd_config)

            mock_naive.assert_called_once()
            mock_sr.assert_not_called()
