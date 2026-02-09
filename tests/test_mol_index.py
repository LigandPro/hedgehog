"""Tests for mol_index.py utilities."""

import pandas as pd

from hedgehog.utils.mol_index import assign_mol_idx


class TestAssignMolIdx:
    """Tests for assign_mol_idx function."""

    def test_basic_assignment(self, tmp_path):
        """Basic mol_idx assignment."""
        df = pd.DataFrame(
            {
                "smiles": ["CCO", "CC", "CCC"],
                "model_name": ["test", "test", "test"],
            }
        )

        result = assign_mol_idx(df, tmp_path)

        assert "mol_idx" in result.columns
        assert len(result["mol_idx"].unique()) == 3

    def test_format_lp_prefix(self, tmp_path):
        """mol_idx should have LP- prefix."""
        df = pd.DataFrame(
            {
                "smiles": ["CCO"],
                "model_name": ["test"],
            }
        )

        result = assign_mol_idx(df, tmp_path)

        assert result["mol_idx"].iloc[0].startswith("LP-")

    def test_unique_per_model(self, tmp_path):
        """Each molecule in each model should have unique index."""
        df = pd.DataFrame(
            {
                "smiles": ["CCO", "CC", "CCC", "CCCC"],
                "model_name": ["a", "a", "b", "b"],
            }
        )

        result = assign_mol_idx(df, tmp_path)

        assert len(result["mol_idx"].unique()) == 4

    def test_model_number_different(self, tmp_path):
        """Different models should have different model numbers in index."""
        df = pd.DataFrame(
            {
                "smiles": ["CCO", "CC"],
                "model_name": ["model_a", "model_b"],
            }
        )

        result = assign_mol_idx(df, tmp_path)

        idx1 = result["mol_idx"].iloc[0]
        idx2 = result["mol_idx"].iloc[1]

        # Model numbers (second part) should be different
        model_num1 = idx1.split("-")[1]
        model_num2 = idx2.split("-")[1]
        assert model_num1 != model_num2

    def test_counter_increments(self, tmp_path):
        """Counter should increment for molecules in same model."""
        df = pd.DataFrame(
            {
                "smiles": ["CCO", "CC", "CCC"],
                "model_name": ["test", "test", "test"],
            }
        )

        result = assign_mol_idx(df, tmp_path)

        counters = [idx.split("-")[2] for idx in result["mol_idx"]]
        assert counters == ["00001", "00002", "00003"]

    def test_preserves_original_columns(self, tmp_path):
        """Original columns should be preserved."""
        df = pd.DataFrame(
            {
                "smiles": ["CCO"],
                "model_name": ["test"],
                "extra_col": [42],
            }
        )

        result = assign_mol_idx(df, tmp_path)

        assert "smiles" in result.columns
        assert "model_name" in result.columns
        assert "extra_col" in result.columns
        assert result["extra_col"].iloc[0] == 42

    def test_no_temp_columns(self, tmp_path):
        """Temporary columns should be removed from output."""
        df = pd.DataFrame(
            {
                "smiles": ["CCO"],
                "model_name": ["test"],
            }
        )

        result = assign_mol_idx(df, tmp_path)

        assert "_model_number" not in result.columns
        assert "_per_model_counter" not in result.columns

    def test_missing_model_name_column(self, tmp_path):
        """Missing model_name should default to 'single'."""
        df = pd.DataFrame(
            {
                "smiles": ["CCO", "CC"],
            }
        )

        result = assign_mol_idx(df, tmp_path)

        assert "model_name" in result.columns
        assert result["model_name"].iloc[0] == "single"

    def test_saves_model_index_map(self, tmp_path):
        """Model index map should be saved to file."""
        df = pd.DataFrame(
            {
                "smiles": ["CCO", "CC"],
                "model_name": ["model_a", "model_b"],
            }
        )

        assign_mol_idx(df, tmp_path)

        map_file = tmp_path / "configs" / "model_index_map.json"
        assert map_file.exists()

    def test_does_not_modify_input(self, tmp_path):
        """Input dataframe should not be modified."""
        df = pd.DataFrame(
            {
                "smiles": ["CCO"],
                "model_name": ["test"],
            }
        )
        original_cols = list(df.columns)

        assign_mol_idx(df, tmp_path)

        assert list(df.columns) == original_cols
        assert "mol_idx" not in df.columns


class TestMolIdxFormat:
    """Tests for mol_idx format validation."""

    def test_format_structure(self, tmp_path):
        """mol_idx should have format LP-XXXX-XXXXX."""
        df = pd.DataFrame(
            {
                "smiles": ["CCO"],
                "model_name": ["test"],
            }
        )

        result = assign_mol_idx(df, tmp_path)
        mol_idx = result["mol_idx"].iloc[0]

        parts = mol_idx.split("-")
        assert len(parts) == 3
        assert parts[0] == "LP"
        assert len(parts[1]) == 4  # model number
        assert len(parts[2]) == 5  # counter

    def test_zero_padding(self, tmp_path):
        """Numbers should be zero-padded."""
        df = pd.DataFrame(
            {
                "smiles": ["CCO"],
                "model_name": ["test"],
            }
        )

        result = assign_mol_idx(df, tmp_path)
        mol_idx = result["mol_idx"].iloc[0]

        parts = mol_idx.split("-")
        # First molecule should be 0001-00001
        assert parts[1] == "0001"
        assert parts[2] == "00001"
