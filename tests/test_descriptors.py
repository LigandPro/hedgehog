"""Tests for descriptors/utils.py."""

import pandas as pd
from rdkit import Chem

from hedgehog.stages.descriptors.utils import (
    _compute_single_molecule_descriptors,
    _get_border_values,
    compute_metrics,
    drop_false_rows,
    order_identity_columns,
)


class TestComputeSingleMoleculeDescriptors:
    """Tests for _compute_single_molecule_descriptors function."""

    def test_benzene_descriptors(self):
        """Compute descriptors for benzene."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        result = _compute_single_molecule_descriptors(mol, "test_model", "test-0")

        assert result is not None
        assert result["model_name"] == "test_model"
        assert result["mol_idx"] == "test-0"
        assert result["n_atoms"] == 12  # 6C + 6H
        assert result["n_heavy_atoms"] == 6
        assert result["n_aroma_rings"] == 1
        assert result["n_rings"] == 1

    def test_ethanol_descriptors(self):
        """Compute descriptors for ethanol."""
        mol = Chem.MolFromSmiles("CCO")
        result = _compute_single_molecule_descriptors(mol, "model", "idx-1")

        assert result is not None
        assert result["n_heavy_atoms"] == 3  # 2C + 1O
        assert result["n_aroma_rings"] == 0
        assert result["n_rings"] == 0
        assert result["hbd"] >= 1  # OH is H-bond donor

    def test_aspirin_descriptors(self):
        """Compute descriptors for aspirin."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = _compute_single_molecule_descriptors(mol, "test", "asp-0")

        assert result is not None
        assert result["n_aroma_rings"] == 1
        assert result["molWt"] > 150  # aspirin MW ~180
        assert result["molWt"] < 200

    def test_caffeine_descriptors(self):
        """Compute descriptors for caffeine."""
        mol = Chem.MolFromSmiles("CN1C=NC2=C1C(=O)N(C)C(=O)N2C")
        result = _compute_single_molecule_descriptors(mol, "test", "caf-0")

        assert result is not None
        assert result["n_N_atoms"] == 4  # caffeine has 4 nitrogens
        assert result["n_rings"] >= 2

    def test_descriptor_keys_present(self):
        """Ensure all expected descriptor keys are present."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        result = _compute_single_molecule_descriptors(mol, "test", "idx-0")

        expected_keys = [
            "model_name",
            "mol_idx",
            "chars",
            "n_atoms",
            "n_heavy_atoms",
            "n_het_atoms",
            "n_N_atoms",
            "fN_atoms",
            "charged_mol",
            "molWt",
            "logP",
            "clogP",
            "sw",
            "ring_size",
            "n_rings",
            "n_aroma_rings",
            "n_fused_aromatic_rings",
            "n_rigid_bonds",
            "n_rot_bonds",
            "hbd",
            "hba",
            "fsp3",
            "tpsa",
            "qed",
        ]
        for key in expected_keys:
            assert key in result, f"Missing key: {key}"


class TestOrderIdentityColumns:
    """Tests for order_identity_columns function."""

    def test_reorder_columns(self):
        """Identity columns should be first."""
        df = pd.DataFrame(
            {
                "other": [1, 2],
                "mol_idx": ["a", "b"],
                "smiles": ["CCO", "CC"],
                "model_name": ["m1", "m2"],
            }
        )
        result = order_identity_columns(df)

        columns = list(result.columns)
        assert columns[0] == "smiles"
        assert columns[1] == "model_name"
        assert columns[2] == "mol_idx"

    def test_missing_identity_columns(self):
        """Should handle missing identity columns by only reordering existing ones."""
        df = pd.DataFrame(
            {
                "smiles": ["CCO"],
                "model_name": ["test"],
                "mol_idx": ["idx-0"],
                "other": [1],
            }
        )
        result = order_identity_columns(df)
        columns = list(result.columns)
        # Identity columns should come first
        assert columns[0] == "smiles"
        assert columns[1] == "model_name"
        assert columns[2] == "mol_idx"


class TestDropFalseRows:
    """Tests for drop_false_rows function."""

    def test_all_pass(self):
        """All rows pass - should return all rows."""
        df = pd.DataFrame(
            {
                "smiles": ["CCO", "CC", "CCC"],
                "molWt_pass": [True, True, True],
                "logP_pass": [True, True, True],
            }
        )
        borders = {}
        result = drop_false_rows(df, borders)
        assert len(result) == 3

    def test_some_fail(self):
        """Some rows fail - should return only passing rows."""
        df = pd.DataFrame(
            {
                "smiles": ["CCO", "CC", "CCC"],
                "molWt_pass": [True, False, True],
                "logP_pass": [True, True, True],
            }
        )
        borders = {}
        result = drop_false_rows(df, borders)
        assert len(result) == 2

    def test_all_fail(self):
        """All rows fail - should return empty dataframe."""
        df = pd.DataFrame(
            {
                "smiles": ["CCO", "CC"],
                "molWt_pass": [False, False],
            }
        )
        borders = {}
        result = drop_false_rows(df, borders)
        assert len(result) == 0


class TestGetBorderValues:
    """Tests for _get_border_values function."""

    def test_get_min_max_borders(self):
        """Should extract min and max values for a column."""
        borders = {
            "molWt_min": 100,
            "molWt_max": 500,
        }
        min_val, max_val = _get_border_values("molWt", borders)
        assert min_val == 100
        assert max_val == 500

    def test_missing_borders(self):
        """Should return None for missing borders."""
        borders = {}
        min_val, max_val = _get_border_values("molWt", borders)
        assert min_val is None
        assert max_val is None

    def test_only_min_border(self):
        """Should handle only min border."""
        borders = {"logP_min": -5}
        min_val, max_val = _get_border_values("logP", borders)
        assert min_val == -5
        assert max_val is None


class TestComputeMetrics:
    """Tests for compute_metrics function."""

    def test_compute_metrics_basic(self, tmp_path):
        """Compute metrics for basic molecules."""
        df = pd.DataFrame(
            {
                "smiles": ["c1ccccc1", "CCO", "CC(=O)O"],
                "model_name": ["test", "test", "test"],
                "mol_idx": ["t-0", "t-1", "t-2"],
            }
        )

        result = compute_metrics(df, str(tmp_path) + "/")

        assert len(result) == 3
        assert "molWt" in result.columns
        assert "logP" in result.columns
        assert "qed" in result.columns

    def test_compute_metrics_with_invalid_smiles(self, tmp_path):
        """Invalid SMILES should be skipped and logged."""
        df = pd.DataFrame(
            {
                "smiles": ["c1ccccc1", "invalid", "CCO"],
                "model_name": ["test", "test", "test"],
                "mol_idx": ["t-0", "t-1", "t-2"],
            }
        )

        result = compute_metrics(df, str(tmp_path) + "/")

        assert len(result) == 2  # Invalid skipped
        # Check skipped file was created
        assert (tmp_path / "skipped_molecules.csv").exists()

    def test_compute_metrics_preserves_model_name(self, tmp_path):
        """Model name should be preserved in output."""
        df = pd.DataFrame(
            {
                "smiles": ["c1ccccc1", "CCO"],
                "model_name": ["model_a", "model_b"],
                "mol_idx": ["a-0", "b-0"],
            }
        )

        result = compute_metrics(df, str(tmp_path) + "/")

        assert "model_name" in result.columns
        assert set(result["model_name"]) == {"model_a", "model_b"}

    def test_compute_metrics_empty_df(self, tmp_path):
        """Empty DataFrame should return empty result."""
        df = pd.DataFrame({"smiles": [], "model_name": [], "mol_idx": []})
        result = compute_metrics(df, str(tmp_path) + "/")

        assert len(result) == 0

    def test_compute_metrics_creates_output_file(self, tmp_path):
        """Output file should be created."""
        df = pd.DataFrame(
            {
                "smiles": ["c1ccccc1"],
                "model_name": ["test"],
                "mol_idx": ["t-0"],
            }
        )

        compute_metrics(df, str(tmp_path) + "/")

        assert (tmp_path / "descriptors_all.csv").exists()


class TestDescriptorValues:
    """Tests for specific descriptor value calculations."""

    def test_molecular_weight_benzene(self):
        """Benzene molecular weight should be ~78."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        result = _compute_single_molecule_descriptors(mol, "test", "idx")

        assert 77 < result["molWt"] < 79

    def test_molecular_weight_aspirin(self):
        """Aspirin molecular weight should be ~180."""
        mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
        result = _compute_single_molecule_descriptors(mol, "test", "idx")

        assert 179 < result["molWt"] < 181

    def test_rotatable_bonds_ethane(self):
        """Ethane should have 0 rotatable bonds."""
        mol = Chem.MolFromSmiles("CC")
        result = _compute_single_molecule_descriptors(mol, "test", "idx")

        assert result["n_rot_bonds"] == 0

    def test_rotatable_bonds_butane(self):
        """Butane should have 1 rotatable bond."""
        mol = Chem.MolFromSmiles("CCCC")
        result = _compute_single_molecule_descriptors(mol, "test", "idx")

        assert result["n_rot_bonds"] == 1

    def test_hydrogen_donors_phenol(self):
        """Phenol should have 1 H-bond donor."""
        mol = Chem.MolFromSmiles("c1ccc(O)cc1")
        result = _compute_single_molecule_descriptors(mol, "test", "idx")

        assert result["hbd"] == 1

    def test_charged_molecule_detection(self):
        """Charged molecule should be detected."""
        # Uncharged benzene
        mol = Chem.MolFromSmiles("c1ccccc1")
        result = _compute_single_molecule_descriptors(mol, "test", "idx")
        assert result["charged_mol"] is True  # No formal charges

        # Carboxylate ion (charged)
        mol = Chem.MolFromSmiles("CC(=O)[O-]")
        result = _compute_single_molecule_descriptors(mol, "test", "idx")
        assert result["charged_mol"] is False  # Has formal charge

    def test_qed_drug_like(self):
        """Drug-like molecules should have QED > 0.3."""
        mol = Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1")  # paracetamol
        result = _compute_single_molecule_descriptors(mol, "test", "idx")

        assert result["qed"] > 0.3


class TestDropFalseRowsAdvanced:
    """Additional tests for drop_false_rows function."""

    def test_charged_mol_filtering_enabled(self):
        """When filter_charged_mol is True, should filter by charged_mol_pass."""
        df = pd.DataFrame(
            {
                "smiles": ["CCO", "CC", "CCC"],
                "molWt_pass": [True, True, True],
                "charged_mol_pass": [True, False, True],
            }
        )
        borders = {"filter_charged_mol": True}
        result = drop_false_rows(df, borders)

        assert len(result) == 2

    def test_charged_mol_filtering_disabled(self):
        """When filter_charged_mol is False, should not filter by charged_mol_pass."""
        df = pd.DataFrame(
            {
                "smiles": ["CCO", "CC", "CCC"],
                "molWt_pass": [True, True, True],
                "charged_mol_pass": [True, False, True],
            }
        )
        borders = {"filter_charged_mol": False}
        result = drop_false_rows(df, borders)

        assert len(result) == 3  # All pass because charged_mol not considered
