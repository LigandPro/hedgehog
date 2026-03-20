"""Tests for descriptors stage modules."""

import pandas as pd
from rdkit import Chem

from hedgehog.descriptors.compute import (
    _compute_single_molecule_descriptors,
    compute_metrics,
)
from hedgehog.descriptors.filtering import (
    _get_border_values,
    drop_false_rows,
    filter_molecules,
)
from hedgehog.descriptors.io import order_identity_columns
from tests.constants import (
    COL_MODEL_NAME,
    COL_MOL_IDX,
    COL_SMILES,
    FILE_FILTERED_MOLECULES,
    MODEL_TEST,
    SMILES_ASPIRIN,
    SMILES_BENZENE,
    SMILES_ETHANOL,
)


class TestComputeSingleMoleculeDescriptors:
    """Tests for _compute_single_molecule_descriptors function."""

    def test_benzene_descriptors(self):
        """Compute descriptors for benzene."""
        mol = Chem.MolFromSmiles(SMILES_BENZENE)
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
        mol = Chem.MolFromSmiles(SMILES_ETHANOL)
        result = _compute_single_molecule_descriptors(mol, "model", "idx-1")

        assert result is not None
        assert result["n_heavy_atoms"] == 3  # 2C + 1O
        assert result["n_aroma_rings"] == 0
        assert result["n_rings"] == 0
        assert result["hbd"] >= 1  # OH is H-bond donor

    def test_aspirin_descriptors(self):
        """Compute descriptors for aspirin."""
        mol = Chem.MolFromSmiles(SMILES_ASPIRIN)
        result = _compute_single_molecule_descriptors(mol, MODEL_TEST, "asp-0")

        assert result is not None
        assert result["n_aroma_rings"] == 1
        assert result["molWt"] > 150  # aspirin MW ~180
        assert result["molWt"] < 200

    def test_caffeine_descriptors(self):
        """Compute descriptors for caffeine."""
        mol = Chem.MolFromSmiles("CN1C=NC2=C1C(=O)N(C)C(=O)N2C")
        result = _compute_single_molecule_descriptors(mol, MODEL_TEST, "caf-0")

        assert result is not None
        assert result["n_N_atoms"] == 4  # caffeine has 4 nitrogens
        assert result["n_rings"] >= 2

    def test_descriptor_keys_present(self):
        """Ensure all expected descriptor keys are present."""
        mol = Chem.MolFromSmiles(SMILES_BENZENE)
        result = _compute_single_molecule_descriptors(mol, MODEL_TEST, "idx-0")

        expected_keys = [
            "model_name",
            "mol_idx",
            "chars",
            "n_atoms",
            "n_heavy_atoms",
            "n_het_atoms",
            "n_N_atoms",
            "n_O_atoms",
            "n_S_atoms",
            "n_NO_atoms",
            "fN_atoms",
            "fNS_atoms",
            "charged_mol",
            "molWt",
            "logP",
            "clogP",
            "sw",
            "ring_size",
            "n_rings",
            "n_small_rings_3_4",
            "max_acyclic_chain_length",
            "n_aroma_rings",
            "n_fused_aromatic_rings",
            "n_rigid_bonds",
            "n_rot_bonds",
            "hbd",
            "hba",
            "fsp3",
            "tpsa",
            "qed",
            ".=O",
            "C2r",
            "C3r",
            "Car",
            "Cs2",
            "Cs3",
            "Csp",
            "Nac",
            "Nd+",
            "Nd0",
            "O_a",
            "O_d",
            "SO2",
            "Sul",
            "Hal",
        ]
        for key in expected_keys:
            assert key in result, f"Missing key: {key}"

    def test_fns_atoms_fraction(self):
        """fNS_atoms should count N+S atoms over heavy atoms."""
        mol = Chem.MolFromSmiles("NS")
        result = _compute_single_molecule_descriptors(mol, MODEL_TEST, "ns-0")

        assert result["n_heavy_atoms"] == 2
        assert result["n_N_atoms"] == 1
        assert result["fN_atoms"] == 0.5
        assert result["fNS_atoms"] == 1.0

    def test_structural_count_descriptors(self):
        """Structural count descriptors should be computed."""
        mol = Chem.MolFromSmiles("CC(=O)N")
        result = _compute_single_molecule_descriptors(mol, MODEL_TEST, "sc-0")

        assert result["n_N_atoms"] == 1
        assert result["n_O_atoms"] == 1
        assert result["n_NO_atoms"] == 2
        assert result["n_small_rings_3_4"] == 0

    def test_max_acyclic_chain_length(self):
        """Maximum acyclic chain length should be computed from non-ring chains."""
        linear = Chem.MolFromSmiles("CCCC")
        ring = Chem.MolFromSmiles("c1ccccc1")

        linear_result = _compute_single_molecule_descriptors(linear, MODEL_TEST, "l-0")
        ring_result = _compute_single_molecule_descriptors(ring, MODEL_TEST, "r-0")

        assert linear_result["max_acyclic_chain_length"] == 4
        assert ring_result["max_acyclic_chain_length"] == 0

    def test_type_alias_counts(self):
        """Configured type aliases should have expected positive counts."""
        cases = [
            ("CC(=O)C", ".=O"),
            ("C1=CCCCC1", "C2r"),
            ("C1CCCCC1", "C3r"),
            (SMILES_BENZENE, "Car"),
            ("C=C", "Cs2"),
            ("CCC", "Cs3"),
            ("C#C", "Csp"),
            ("n1ccccc1", "Nac"),
            ("C[NH3+]", "Nd+"),
            ("CN", "Nd0"),
            ("COC", "O_a"),
            (SMILES_ETHANOL, "O_d"),
            ("CS(=O)(=O)C", "SO2"),
            ("CSC", "Sul"),
            ("CCCl", "Hal"),
        ]

        for smiles, alias in cases:
            mol = Chem.MolFromSmiles(smiles)
            result = _compute_single_molecule_descriptors(mol, MODEL_TEST, f"{alias}-0")
            assert result[alias] >= 1, (
                f"Expected positive count for {alias} in {smiles}"
            )


class TestOrderIdentityColumns:
    """Tests for order_identity_columns function."""

    def test_reorder_columns(self):
        """Identity columns should be first."""
        df = pd.DataFrame(
            {
                "other": [1, 2],
                COL_MOL_IDX: ["a", "b"],
                COL_SMILES: [SMILES_ETHANOL, "CC"],
                COL_MODEL_NAME: ["m1", "m2"],
            }
        )
        result = order_identity_columns(df)

        columns = list(result.columns)
        assert columns[0] == COL_SMILES
        assert columns[1] == COL_MODEL_NAME
        assert columns[2] == COL_MOL_IDX

    def test_missing_identity_columns(self):
        """Should handle missing identity columns by only reordering existing ones."""
        df = pd.DataFrame(
            {
                COL_SMILES: [SMILES_ETHANOL],
                COL_MODEL_NAME: [MODEL_TEST],
                COL_MOL_IDX: ["idx-0"],
                "other": [1],
            }
        )
        result = order_identity_columns(df)
        columns = list(result.columns)
        # Identity columns should come first
        assert columns[0] == COL_SMILES
        assert columns[1] == COL_MODEL_NAME
        assert columns[2] == COL_MOL_IDX


class TestDropFalseRows:
    """Tests for drop_false_rows function."""

    def test_all_pass(self):
        """All rows pass - should return all rows."""
        df = pd.DataFrame(
            {
                COL_SMILES: [SMILES_ETHANOL, "CC", "CCC"],
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
                COL_SMILES: [SMILES_ETHANOL, "CC", "CCC"],
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
                COL_SMILES: [SMILES_ETHANOL, "CC"],
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
                COL_SMILES: [SMILES_BENZENE, SMILES_ETHANOL, "CC(=O)O"],
                COL_MODEL_NAME: [MODEL_TEST, MODEL_TEST, MODEL_TEST],
                COL_MOL_IDX: ["t-0", "t-1", "t-2"],
            }
        )

        result = compute_metrics(df, str(tmp_path) + "/")

        assert len(result) == 3
        assert "molWt" in result.columns

    def test_compute_metrics_with_invalid_smiles(self, tmp_path):
        """Invalid SMILES should be skipped and logged."""
        df = pd.DataFrame(
            {
                COL_SMILES: [SMILES_BENZENE, "invalid", SMILES_ETHANOL],
                COL_MODEL_NAME: [MODEL_TEST, MODEL_TEST, MODEL_TEST],
                COL_MOL_IDX: ["t-0", "t-1", "t-2"],
            }
        )

        result = compute_metrics(df, str(tmp_path) + "/")

        assert len(result) == 2  # Invalid skipped
        # Check skipped file was created
        assert (tmp_path / "skipped_molecules.csv").exists()

    def test_compute_metrics_remove_stereochemistry(self, tmp_path):
        """When enabled, stereochemistry should be removed from output SMILES."""
        df = pd.DataFrame(
            {
                COL_SMILES: ["F[C@H](Cl)Br"],
                COL_MODEL_NAME: [MODEL_TEST],
                COL_MOL_IDX: ["t-0"],
            }
        )
        cfg = {"preprocess": {"remove_stereochemistry": True}}
        result = compute_metrics(df, str(tmp_path) + "/", config_descriptors=cfg)

        assert len(result) == 1
        assert "@" not in result.iloc[0][COL_SMILES]

    def test_compute_metrics_remove_charges_best_effort(self, tmp_path):
        """When enabled, charged molecules should be neutralized or skipped."""
        df = pd.DataFrame(
            {
                COL_SMILES: ["C[NH3+]"],
                COL_MODEL_NAME: [MODEL_TEST],
                COL_MOL_IDX: ["t-0"],
            }
        )
        cfg = {"preprocess": {"remove_charges": True}}
        result = compute_metrics(df, str(tmp_path) + "/", config_descriptors=cfg)

        # Best-effort behavior: if neutralization succeeds, '+' is removed.
        if len(result) == 1:
            assert "+" not in result.iloc[0][COL_SMILES]
        else:
            assert (tmp_path / "skipped_molecules.csv").exists()

    def test_compute_metrics_reject_radicals(self, tmp_path):
        """When enabled, radicals should be skipped."""
        df = pd.DataFrame(
            {
                COL_SMILES: ["[CH3]"],
                COL_MODEL_NAME: [MODEL_TEST],
                COL_MOL_IDX: ["t-0"],
            }
        )
        cfg = {"preprocess": {"remove_radicals": True}}
        result = compute_metrics(df, str(tmp_path) + "/", config_descriptors=cfg)

        assert len(result) == 0
        assert (tmp_path / "skipped_molecules.csv").exists()

    def test_compute_metrics_preserves_model_name(self, tmp_path):
        """Model name should be preserved in output."""
        df = pd.DataFrame(
            {
                COL_SMILES: [SMILES_BENZENE, SMILES_ETHANOL],
                COL_MODEL_NAME: ["model_a", "model_b"],
                COL_MOL_IDX: ["a-0", "b-0"],
            }
        )

        result = compute_metrics(df, str(tmp_path) + "/")

        assert COL_MODEL_NAME in result.columns
        assert set(result[COL_MODEL_NAME]) == {"model_a", "model_b"}

    def test_compute_metrics_empty_df(self, tmp_path):
        """Empty DataFrame should return empty result."""
        df = pd.DataFrame({COL_SMILES: [], COL_MODEL_NAME: [], COL_MOL_IDX: []})
        result = compute_metrics(df, str(tmp_path) + "/")

        assert len(result) == 0

    def test_compute_metrics_creates_output_file(self, tmp_path):
        """Output file should be created."""
        df = pd.DataFrame(
            {
                COL_SMILES: [SMILES_BENZENE],
                COL_MODEL_NAME: [MODEL_TEST],
                COL_MOL_IDX: ["t-0"],
            }
        )

        compute_metrics(df, str(tmp_path) + "/")

        assert (tmp_path / "descriptors_all.csv").exists()


class TestFilterMolecules:
    """Tests for descriptor-based filtering."""

    def test_allowed_chars_filters_out_salts(self, tmp_path):
        mol_ok = Chem.MolFromSmiles(SMILES_ETHANOL)
        mol_salt = Chem.MolFromSmiles("CCO.[Na+]")

        row_ok = _compute_single_molecule_descriptors(mol_ok, "m", "ok-0")
        row_ok[COL_SMILES] = SMILES_ETHANOL
        row_salt = _compute_single_molecule_descriptors(mol_salt, "m", "salt-0")
        row_salt[COL_SMILES] = "CCO.[Na+]"

        df = pd.DataFrame([row_ok, row_salt])

        borders = {
            "allowed_chars": ["C", "N", "O", "S", "F", "Cl", "Br", "H"],
        }
        filter_molecules(df, borders, str(tmp_path))

        passed = pd.read_csv(tmp_path / FILE_FILTERED_MOLECULES)
        assert len(passed) == 1
        assert passed.iloc[0][COL_SMILES] == SMILES_ETHANOL

        flags = pd.read_csv(tmp_path / "pass_flags.csv")
        assert "chars_pass" in flags.columns

    def test_charged_mol_filter_applies_when_enabled(self, tmp_path):
        mol_neutral = Chem.MolFromSmiles(SMILES_ETHANOL)
        mol_charged = Chem.MolFromSmiles("CC(=O)[O-]")

        row_neutral = _compute_single_molecule_descriptors(mol_neutral, "m", "n-0")
        row_neutral[COL_SMILES] = SMILES_ETHANOL
        row_charged = _compute_single_molecule_descriptors(mol_charged, "m", "c-0")
        row_charged[COL_SMILES] = "CC(=O)[O-]"

        df = pd.DataFrame([row_neutral, row_charged])

        borders = {
            "charged_mol_allowed": False,
            "filter_charged_mol": True,
        }
        filter_molecules(df, borders, str(tmp_path))

        passed = pd.read_csv(tmp_path / FILE_FILTERED_MOLECULES)
        assert len(passed) == 1
        assert passed.iloc[0][COL_SMILES] == SMILES_ETHANOL

        flags = pd.read_csv(tmp_path / "pass_flags.csv")
        assert "charged_mol_pass" in flags.columns

    def test_structural_constraints_nested_block(self, tmp_path):
        """Nested structural_constraints should be applied via generic pass flags."""
        mol_pass = Chem.MolFromSmiles(SMILES_ETHANOL)  # O_d=1, no SO2
        mol_fail = Chem.MolFromSmiles("CS(=O)(=O)C")  # SO2=1, n_O_atoms=2

        row_pass = _compute_single_molecule_descriptors(mol_pass, "m", "p-0")
        row_pass[COL_SMILES] = SMILES_ETHANOL
        row_fail = _compute_single_molecule_descriptors(mol_fail, "m", "f-0")
        row_fail[COL_SMILES] = "CS(=O)(=O)C"

        df = pd.DataFrame([row_pass, row_fail])

        config_block = {
            "borders": {"molWt_min": 0},
            "structural_constraints": {
                "enabled": True,
                "type_limits": {
                    "SO2": 0,
                    "O_d": 1,
                },
                "element_limits": {
                    "O": 1,
                    "S": 0,
                },
                "max_n_or_o_atoms": 1,
                "max_small_rings_3_4": 0,
                "max_acyclic_chain_length": 3,
            },
        }
        filter_molecules(df, config_block, str(tmp_path))

        passed = pd.read_csv(tmp_path / FILE_FILTERED_MOLECULES)
        assert len(passed) == 1
        assert passed.iloc[0][COL_SMILES] == SMILES_ETHANOL

        flags = pd.read_csv(tmp_path / "pass_flags.csv")
        assert "SO2_pass" in flags.columns
        assert "n_O_atoms_pass" in flags.columns
        assert "n_NO_atoms_pass" in flags.columns
        assert "max_acyclic_chain_length_pass" in flags.columns

    def test_structural_constraints_disabled(self, tmp_path):
        """Disabled structural_constraints should not affect filtering."""
        mol_a = Chem.MolFromSmiles(SMILES_ETHANOL)
        mol_b = Chem.MolFromSmiles("CS(=O)(=O)C")

        row_a = _compute_single_molecule_descriptors(mol_a, "m", "a-0")
        row_a[COL_SMILES] = SMILES_ETHANOL
        row_b = _compute_single_molecule_descriptors(mol_b, "m", "b-0")
        row_b[COL_SMILES] = "CS(=O)(=O)C"

        df = pd.DataFrame([row_a, row_b])

        config_block = {
            "borders": {},
            "structural_constraints": {
                "enabled": False,
                "type_limits": {"SO2": 0},
                "max_n_or_o_atoms": 0,
            },
        }
        filter_molecules(df, config_block, str(tmp_path))

        passed = pd.read_csv(tmp_path / FILE_FILTERED_MOLECULES)
        assert len(passed) == 2


class TestDescriptorValues:
    """Tests for specific descriptor value calculations."""

    def test_molecular_weight_benzene(self):
        """Benzene molecular weight should be ~78."""
        mol = Chem.MolFromSmiles(SMILES_BENZENE)
        result = _compute_single_molecule_descriptors(mol, MODEL_TEST, "idx")

        assert 77 < result["molWt"] < 79

    def test_molecular_weight_aspirin(self):
        """Aspirin molecular weight should be ~180."""
        mol = Chem.MolFromSmiles(SMILES_ASPIRIN)
        result = _compute_single_molecule_descriptors(mol, MODEL_TEST, "idx")

        assert 179 < result["molWt"] < 181

    def test_rotatable_bonds_ethane(self):
        """Ethane should have 0 rotatable bonds."""
        mol = Chem.MolFromSmiles("CC")
        result = _compute_single_molecule_descriptors(mol, MODEL_TEST, "idx")

        assert result["n_rot_bonds"] == 0

    def test_rotatable_bonds_butane(self):
        """Butane should have 1 rotatable bond."""
        mol = Chem.MolFromSmiles("CCCC")
        result = _compute_single_molecule_descriptors(mol, MODEL_TEST, "idx")

        assert result["n_rot_bonds"] == 1

    def test_hydrogen_donors_phenol(self):
        """Phenol should have 1 H-bond donor."""
        mol = Chem.MolFromSmiles("c1ccc(O)cc1")
        result = _compute_single_molecule_descriptors(mol, MODEL_TEST, "idx")

        assert result["hbd"] == 1

    def test_charged_molecule_detection(self):
        """Charged molecule should be detected."""
        # Uncharged benzene
        mol = Chem.MolFromSmiles("c1ccccc1")
        result = _compute_single_molecule_descriptors(mol, MODEL_TEST, "idx")
        assert result["charged_mol"] is True  # No formal charges

        # Carboxylate ion (charged)
        mol = Chem.MolFromSmiles("CC(=O)[O-]")
        result = _compute_single_molecule_descriptors(mol, MODEL_TEST, "idx")
        assert result["charged_mol"] is False  # Has formal charge

    def test_qed_drug_like(self):
        """Drug-like molecules should have QED > 0.3."""
        mol = Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1")  # paracetamol
        result = _compute_single_molecule_descriptors(mol, MODEL_TEST, "idx")

        assert result["qed"] > 0.3


class TestDropFalseRowsAdvanced:
    """Additional tests for drop_false_rows function."""

    def test_charged_mol_filtering_enabled(self):
        """When filter_charged_mol is True, should filter by charged_mol_pass."""
        df = pd.DataFrame(
            {
                COL_SMILES: [SMILES_ETHANOL, "CC", "CCC"],
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
                COL_SMILES: [SMILES_ETHANOL, "CC", "CCC"],
                "molWt_pass": [True, True, True],
                "charged_mol_pass": [True, False, True],
            }
        )
        borders = {"filter_charged_mol": False}
        result = drop_false_rows(df, borders)

        assert len(result) == 3  # All pass because charged_mol not considered
