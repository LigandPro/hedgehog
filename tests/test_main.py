"""Tests for main.py utilities."""

from pathlib import Path

import pandas as pd

from hedgehog.main import (
    Stage,
    _canonicalize_smiles,
    _folder_is_empty,
    _get_input_format_flag,
    _get_unique_results_folder,
    _validate_input_path,
    preprocess_input_with_rdkit,
)


class TestCanonicalizeSmiles:
    """Tests for _canonicalize_smiles function."""

    def test_valid_smiles_benzene(self):
        """Valid benzene SMILES - should return canonical form."""
        result = _canonicalize_smiles("c1ccccc1")
        assert result is not None
        assert result == "c1ccccc1"

    def test_valid_smiles_ethanol(self):
        """Valid ethanol SMILES - should return canonical form."""
        result = _canonicalize_smiles("CCO")
        assert result is not None
        assert result == "CCO"

    def test_valid_smiles_aspirin(self):
        """Valid aspirin SMILES - should return canonical form."""
        result = _canonicalize_smiles("CC(=O)Oc1ccccc1C(=O)O")
        assert result is not None

    def test_invalid_smiles(self):
        """Invalid SMILES string - should return None."""
        result = _canonicalize_smiles("invalid")
        assert result is None

    def test_empty_smiles(self):
        """Empty SMILES string - should return empty string."""
        result = _canonicalize_smiles("")
        assert result == ""

    def test_smiles_with_stereochemistry(self):
        """SMILES with stereochemistry should be preserved."""
        result = _canonicalize_smiles("[C@@H](O)(F)Cl")
        assert result is not None

    def test_malformed_smiles(self):
        """Malformed SMILES - should return None."""
        result = _canonicalize_smiles("C(C)(C)(C)(C)C")  # invalid valence
        assert result is None


class TestFolderIsEmpty:
    """Tests for _folder_is_empty function."""

    def test_nonexistent_folder(self, tmp_path):
        """Nonexistent folder - should return True."""
        nonexistent = tmp_path / "nonexistent"
        assert _folder_is_empty(nonexistent) is True

    def test_empty_folder(self, tmp_path):
        """Empty folder - should return True."""
        empty_folder = tmp_path / "empty"
        empty_folder.mkdir()
        assert _folder_is_empty(empty_folder) is True

    def test_folder_with_file(self, tmp_path):
        """Folder with file - should return False."""
        folder = tmp_path / "with_file"
        folder.mkdir()
        (folder / "file.txt").touch()
        assert _folder_is_empty(folder) is False

    def test_folder_with_subdirectory(self, tmp_path):
        """Folder with subdirectory - should return False."""
        folder = tmp_path / "with_subdir"
        folder.mkdir()
        (folder / "subdir").mkdir()
        assert _folder_is_empty(folder) is False


class TestGetUniqueResultsFolder:
    """Tests for _get_unique_results_folder function."""

    def test_empty_folder(self, tmp_path):
        """If no numbered folders exist, return base_1."""
        result = _get_unique_results_folder(tmp_path / "results")
        assert result == tmp_path / "results_1"

    def test_nonexistent_folder(self, tmp_path):
        """If parent folder doesn't exist, return base_1."""
        result = _get_unique_results_folder(tmp_path / "new_results")
        assert result == tmp_path / "new_results_1"

    def test_folder_exists_with_content(self, tmp_path):
        """If folder exists with content but no numbered folders, return _1."""
        results = tmp_path / "results"
        results.mkdir()
        (results / "file.txt").touch()

        result = _get_unique_results_folder(results)
        assert result == tmp_path / "results_1"

    def test_multiple_existing_folders(self, tmp_path):
        """If multiple numbered folders exist, return next available."""
        # Create results_1 and results_2
        (tmp_path / "results_1").mkdir()
        (tmp_path / "results_2").mkdir()

        result = _get_unique_results_folder(tmp_path / "results")
        assert result == tmp_path / "results_3"

    def test_folder_already_numbered(self, tmp_path):
        """Sequential numbering continues from highest existing number."""
        # Create results_5
        (tmp_path / "results_5").mkdir()

        result = _get_unique_results_folder(tmp_path / "results")
        assert result == tmp_path / "results_6"


class TestValidateInputPath:
    """Tests for _validate_input_path function."""

    def test_existing_file(self, tmp_path):
        """Existing file should return Path object."""
        test_file = tmp_path / "test.csv"
        test_file.write_text("smiles\nCCO")

        result = _validate_input_path(str(test_file))
        assert result is not None
        assert result.exists()

    def test_nonexistent_file(self, tmp_path):
        """Nonexistent file should return None."""
        result = _validate_input_path(str(tmp_path / "nonexistent.csv"))
        assert result is None

    def test_glob_pattern(self):
        """Glob pattern should return None."""
        result = _validate_input_path("/path/to/*.csv")
        assert result is None

    def test_question_mark_pattern(self):
        """Pattern with ? should return None."""
        result = _validate_input_path("/path/to/file?.csv")
        assert result is None


class TestGetInputFormatFlag:
    """Tests for _get_input_format_flag function."""

    def test_csv_extension(self):
        """CSV extension should return -icsv."""
        assert _get_input_format_flag("csv") == "-icsv"
        assert _get_input_format_flag(".csv") == "-icsv"
        assert _get_input_format_flag("CSV") == "-icsv"

    def test_smi_extension(self):
        """SMI extension should return -ismi."""
        assert _get_input_format_flag("smi") == "-ismi"
        assert _get_input_format_flag("ismi") == "-ismi"
        assert _get_input_format_flag("txt") == "-ismi"

    def test_unsupported_extension(self):
        """Unsupported extension should return None."""
        assert _get_input_format_flag("pdf") is None
        assert _get_input_format_flag("xyz") is None


class TestPreprocessInputWithRdkit:
    """Tests for preprocess_input_with_rdkit function."""

    def test_valid_input(self, tmp_path, mock_logger):
        """Preprocess valid input file."""
        input_file = tmp_path / "input.csv"
        input_file.write_text("smiles,model_name\nCCO,test\nc1ccccc1,test")

        output_folder = tmp_path / "output"
        result = preprocess_input_with_rdkit(
            str(input_file), output_folder, mock_logger
        )

        assert result is not None
        assert Path(result).exists()

    def test_removes_duplicates(self, tmp_path, mock_logger):
        """Should remove duplicate SMILES within models."""
        input_file = tmp_path / "input.csv"
        input_file.write_text("smiles,model_name\nCCO,test\nCCO,test\nc1ccccc1,test")

        output_folder = tmp_path / "output"
        result = preprocess_input_with_rdkit(
            str(input_file), output_folder, mock_logger
        )

        output_df = pd.read_csv(result)
        assert len(output_df) == 2  # Duplicates removed

    def test_glob_pattern_returns_none(self, tmp_path, mock_logger):
        """Glob pattern input should return None."""
        result = preprocess_input_with_rdkit("*.csv", tmp_path, mock_logger)
        assert result is None

    def test_invalid_smiles_filtered(self, tmp_path, mock_logger):
        """Invalid SMILES should be filtered out."""
        input_file = tmp_path / "input.csv"
        input_file.write_text(
            "smiles,model_name\nCCO,test\ninvalid,test\nc1ccccc1,test"
        )

        output_folder = tmp_path / "output"
        result = preprocess_input_with_rdkit(
            str(input_file), output_folder, mock_logger
        )

        output_df = pd.read_csv(result)
        assert len(output_df) == 2  # Invalid SMILES filtered


class TestStageEnum:
    """Tests for Stage enum."""

    def test_stage_values(self):
        """Test stage enum values."""
        assert Stage.descriptors.value == "descriptors"
        assert Stage.struct_filters.value == "struct_filters"
        assert Stage.synthesis.value == "synthesis"
        assert Stage.docking.value == "docking"

    def test_stage_description(self):
        """Test stage descriptions."""
        assert "descriptors" in Stage.descriptors.description.lower()
        assert "filter" in Stage.struct_filters.description.lower()
        assert "synth" in Stage.synthesis.description.lower()
        assert "docking" in Stage.docking.description.lower()
