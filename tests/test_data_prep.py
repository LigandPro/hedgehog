"""Tests for data_prep.py utilities."""

import pandas as pd
import pytest

from hedgehog.utils.data_prep import (
    MODE_MULTI,
    MODE_SINGLE,
    _apply_sampling,
    _detect_mode_and_paths,
    _extract_model_name_from_path,
    _file_has_multiple_models,
    _find_column_case_insensitive,
    _normalize_model_name_column,
    _normalize_smiles_column,
    _read_csv_with_fallback,
    _remove_duplicates,
)


class TestFindColumnCaseInsensitive:
    """Tests for _find_column_case_insensitive function."""

    def test_exact_match(self):
        """Exact match - should return the column name."""
        df = pd.DataFrame({"smiles": [], "Name": []})
        assert _find_column_case_insensitive(df, "smiles") == "smiles"

    def test_case_insensitive_match(self):
        """Case-insensitive match - should return original column name."""
        df = pd.DataFrame({"SMILES": [], "Name": []})
        assert _find_column_case_insensitive(df, "smiles") == "SMILES"

    def test_mixed_case_match(self):
        """Mixed case match - should return original column name."""
        df = pd.DataFrame({"SmIlEs": [], "Name": []})
        assert _find_column_case_insensitive(df, "smiles") == "SmIlEs"

    def test_missing_column(self):
        """Column not found - should return None."""
        df = pd.DataFrame({"col1": [], "col2": []})
        assert _find_column_case_insensitive(df, "smiles") is None

    def test_empty_dataframe(self):
        """Empty DataFrame - should return None."""
        df = pd.DataFrame()
        assert _find_column_case_insensitive(df, "smiles") is None

    def test_model_name_column(self):
        """Test finding model_name column."""
        df = pd.DataFrame({"smiles": [], "MODEL_NAME": []})
        assert _find_column_case_insensitive(df, "model_name") == "MODEL_NAME"


class TestApplySampling:
    """Tests for _apply_sampling function."""

    def test_sampling_smaller_than_data(self):
        """Sample size smaller than data - should return sampled data."""
        df = pd.DataFrame({"a": range(100)})
        result, warning = _apply_sampling(df, 10, model_name="test")
        assert len(result) == 10
        assert warning is None

    def test_sampling_larger_than_data(self):
        """Sample size larger than data - should return all data with warning."""
        df = pd.DataFrame({"a": range(5)})
        result, warning = _apply_sampling(df, 10, model_name="test")
        assert len(result) == 5
        assert warning is not None
        assert warning["requested"] == 10
        assert warning["available"] == 5

    def test_no_sampling(self):
        """No sample size - should return original data."""
        df = pd.DataFrame({"a": range(10)})
        result, warning = _apply_sampling(df, None)
        assert len(result) == 10
        assert warning is None

    def test_sample_size_equals_data_size(self):
        """Sample size equals data size - should return all data."""
        df = pd.DataFrame({"a": range(10)})
        result, warning = _apply_sampling(df, 10)
        assert len(result) == 10
        assert warning is None

    def test_reproducibility(self):
        """Sampling should be reproducible (random_state=42)."""
        df = pd.DataFrame({"a": range(100)})
        result1, _ = _apply_sampling(df, 10)
        result2, _ = _apply_sampling(df, 10)
        pd.testing.assert_frame_equal(
            result1.reset_index(drop=True), result2.reset_index(drop=True)
        )


class TestNormalizeSmilesColumn:
    """Tests for _normalize_smiles_column function."""

    def test_already_normalized(self):
        """Column already named 'smiles' - should return unchanged."""
        df = pd.DataFrame({"smiles": ["CCO", "CC"]})
        result = _normalize_smiles_column(df)
        assert "smiles" in result.columns

    def test_uppercase_smiles(self):
        """SMILES column - should be renamed to 'smiles'."""
        df = pd.DataFrame({"SMILES": ["CCO", "CC"]})
        result = _normalize_smiles_column(df)
        assert "smiles" in result.columns
        assert "SMILES" not in result.columns

    def test_no_smiles_column(self):
        """No smiles column - should return unchanged."""
        df = pd.DataFrame({"other": ["a", "b"]})
        result = _normalize_smiles_column(df)
        assert "smiles" not in result.columns


class TestNormalizeModelNameColumn:
    """Tests for _normalize_model_name_column function."""

    def test_already_normalized(self):
        """Column already named 'model_name' - should return unchanged."""
        df = pd.DataFrame({"smiles": ["CCO"], "model_name": ["test"]})
        result = _normalize_model_name_column(df, "/path/to/file.csv")
        assert "model_name" in result.columns

    def test_uppercase_model_name(self):
        """MODEL_NAME column - should be renamed."""
        df = pd.DataFrame({"smiles": ["CCO"], "MODEL_NAME": ["test"]})
        result = _normalize_model_name_column(df, "/path/to/file.csv")
        assert "model_name" in result.columns

    def test_name_column(self):
        """name column - should be renamed to model_name."""
        df = pd.DataFrame({"smiles": ["CCO"], "name": ["test"]})
        result = _normalize_model_name_column(df, "/path/to/file.csv")
        assert "model_name" in result.columns

    def test_no_model_name_column(self):
        """No model_name column - should extract from path."""
        df = pd.DataFrame({"smiles": ["CCO"]})
        result = _normalize_model_name_column(df, "/path/to/my_model.csv")
        assert "model_name" in result.columns
        assert result["model_name"].iloc[0] == "my_model"


class TestDetectModeAndPaths:
    """Tests for _detect_mode_and_paths function."""

    def test_single_file_mode(self, tmp_path):
        """Single file should be detected as single comparison."""
        test_file = tmp_path / "test.csv"
        test_file.write_text("smiles\nCCO\nCC")

        mode, paths = _detect_mode_and_paths(str(test_file))

        assert mode == MODE_SINGLE
        assert len(paths) == 1

    def test_glob_pattern_multiple_files(self, tmp_path):
        """Glob pattern matching multiple files should be multi comparison."""
        (tmp_path / "model_a.csv").write_text("smiles\nCCO")
        (tmp_path / "model_b.csv").write_text("smiles\nCC")

        mode, paths = _detect_mode_and_paths(str(tmp_path / "model_*.csv"))

        assert mode == MODE_MULTI
        assert len(paths) == 2

    def test_file_with_multiple_models(self, tmp_path):
        """Single file with multiple models should be multi comparison."""
        test_file = tmp_path / "test.csv"
        test_file.write_text("smiles,model_name\nCCO,model_a\nCC,model_b")

        mode, paths = _detect_mode_and_paths(str(test_file))

        assert mode == MODE_MULTI
        assert len(paths) == 1

    def test_nonexistent_file_raises_error(self, tmp_path):
        """Nonexistent file should raise FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            _detect_mode_and_paths(str(tmp_path / "nonexistent.csv"))

    def test_glob_no_matches_raises_error(self, tmp_path):
        """Glob pattern with no matches should raise FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            _detect_mode_and_paths(str(tmp_path / "*.xyz"))


class TestFileHasMultipleModels:
    """Tests for _file_has_multiple_models function."""

    def test_multiple_models(self, tmp_path):
        """File with multiple model names should return True."""
        test_file = tmp_path / "test.csv"
        test_file.write_text("smiles,model_name\nCCO,model_a\nCC,model_b\nCCC,model_c")

        assert _file_has_multiple_models(str(test_file)) is True

    def test_single_model(self, tmp_path):
        """File with single model name should return False."""
        test_file = tmp_path / "test.csv"
        test_file.write_text("smiles,model_name\nCCO,model_a\nCC,model_a")

        assert _file_has_multiple_models(str(test_file)) is False

    def test_no_model_column(self, tmp_path):
        """File without model_name column should return False."""
        test_file = tmp_path / "test.csv"
        test_file.write_text("smiles\nCCO\nCC")

        assert _file_has_multiple_models(str(test_file)) is False


class TestRemoveDuplicates:
    """Tests for _remove_duplicates function."""

    def test_removes_duplicates_with_model(self, mock_logger):
        """Removes duplicates when model_name column present."""
        df = pd.DataFrame(
            {
                "smiles": ["CCO", "CCO", "CC"],
                "model_name": ["a", "a", "a"],
            }
        )
        result = _remove_duplicates(df, mock_logger)

        assert len(result) == 2

    def test_removes_duplicates_without_model(self, mock_logger):
        """Removes duplicates when model_name column absent."""
        df = pd.DataFrame({"smiles": ["CCO", "CCO", "CC"]})
        result = _remove_duplicates(df, mock_logger)

        assert len(result) == 2

    def test_keeps_different_model_duplicates(self, mock_logger):
        """Same SMILES in different models should be kept."""
        df = pd.DataFrame(
            {
                "smiles": ["CCO", "CCO", "CC"],
                "model_name": ["a", "b", "a"],
            }
        )
        result = _remove_duplicates(df, mock_logger)

        assert len(result) == 3  # All kept because different models


class TestExtractModelNameFromPath:
    """Tests for _extract_model_name_from_path function."""

    def test_simple_path(self):
        """Extract model name from simple path."""
        result = _extract_model_name_from_path("/path/to/my_model.csv")
        assert result == "my_model"

    def test_nested_path(self):
        """Extract model name from nested path."""
        result = _extract_model_name_from_path("/a/b/c/d/model_name.csv")
        assert result == "model_name"

    def test_multiple_dots(self):
        """Handle filename with multiple dots."""
        result = _extract_model_name_from_path("/path/to/my.model.csv")
        assert result == "my"


class TestReadCsvWithFallback:
    """Tests for _read_csv_with_fallback function."""

    def test_standard_csv(self, tmp_path):
        """Read standard CSV with header."""
        test_file = tmp_path / "test.csv"
        test_file.write_text("smiles,model_name\nCCO,test")

        result = _read_csv_with_fallback(str(test_file))

        assert "smiles" in result.columns
        assert "model_name" in result.columns

    def test_headerless_csv(self, tmp_path):
        """Read headerless CSV - first row becomes header unless parsing fails."""
        test_file = tmp_path / "test.csv"
        test_file.write_text("CCO\nCC\nCCC")

        result = _read_csv_with_fallback(str(test_file))

        # First row 'CCO' becomes column name
        assert "CCO" in result.columns
        assert len(result) == 2  # CC and CCC are data rows
