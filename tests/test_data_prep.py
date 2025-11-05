"""Unit tests for data preparation module."""

from unittest.mock import Mock

import pandas as pd
import pytest

from hedge.utils.data_prep import (
    _apply_sampling,
    _detect_mode_and_paths,
    _extract_model_name_from_path,
    _find_smiles_column,
    _load_multi_comparison_data,
    _load_single_comparison_data,
    _normalize_model_name_column,
    _normalize_smiles_column,
    prepare_input_data,
)


@pytest.mark.unit
def test_find_smiles_column_lowercase():
    """Test finding SMILES column with lowercase name."""
    df = pd.DataFrame({"smiles": ["CCO", "CCN"], "other": [1, 2]})
    result = _find_smiles_column(df)
    assert result == "smiles"


@pytest.mark.unit
def test_find_smiles_column_uppercase():
    """Test finding SMILES column with uppercase name."""
    df = pd.DataFrame({"SMILES": ["CCO", "CCN"], "other": [1, 2]})
    result = _find_smiles_column(df)
    assert result == "SMILES"


@pytest.mark.unit
def test_find_smiles_column_not_found():
    """Test when SMILES column doesn't exist."""
    df = pd.DataFrame({"molecules": ["CCO", "CCN"], "other": [1, 2]})
    result = _find_smiles_column(df)
    assert result is None


@pytest.mark.unit
def test_normalize_smiles_column():
    """Test normalizing SMILES column name."""
    df = pd.DataFrame({"SMILES": ["CCO", "CCN"], "other": [1, 2]})
    result = _normalize_smiles_column(df, "test.csv")
    assert "smiles" in result.columns
    assert "SMILES" not in result.columns


@pytest.mark.unit
def test_extract_model_name_from_path():
    """Test extracting model name from file path."""
    path = "/path/to/data/model_123.csv"
    result = _extract_model_name_from_path(path)
    assert result == "model_123"


@pytest.mark.unit
def test_extract_model_name_from_simple_path():
    """Test extracting model name from simple path."""
    path = "data.csv"
    result = _extract_model_name_from_path(path)
    assert result == "data"


@pytest.mark.unit
def test_normalize_model_name_column_with_model_name():
    """Test normalizing when model_name column exists."""
    df = pd.DataFrame({"smiles": ["CCO"], "model_name": ["test"]})
    result = _normalize_model_name_column(df, "test.csv")
    assert "model_name" in result.columns
    assert result["model_name"].iloc[0] == "test"


@pytest.mark.unit
def test_normalize_model_name_column_with_name():
    """Test normalizing when 'name' column exists."""
    df = pd.DataFrame({"smiles": ["CCO"], "name": ["test"]})
    result = _normalize_model_name_column(df, "test.csv")
    assert "model_name" in result.columns
    assert result["model_name"].iloc[0] == "test"


@pytest.mark.unit
def test_normalize_model_name_column_missing():
    """Test normalizing when model_name column is missing."""
    df = pd.DataFrame({"smiles": ["CCO"], "other": [1]})
    result = _normalize_model_name_column(df, "/path/to/mymodel.csv")
    assert "model_name" in result.columns
    assert result["model_name"].iloc[0] == "mymodel"


@pytest.mark.unit
def test_apply_sampling_exact():
    """Test applying exact sampling."""
    df = pd.DataFrame({"smiles": ["CCO", "CCN", "CCC", "CCF", "CCCl"]})
    mock_logger = Mock()
    result = _apply_sampling(df, sample_size=3, logger=mock_logger)
    assert len(result) == 3


@pytest.mark.unit
def test_apply_sampling_exceeds_data():
    """Test sampling when sample_size exceeds data size."""
    df = pd.DataFrame({"smiles": ["CCO", "CCN"]})
    mock_logger = Mock()
    result = _apply_sampling(df, sample_size=10, logger=mock_logger)
    assert len(result) == 2
    mock_logger.warning.assert_called_once()


@pytest.mark.unit
def test_apply_sampling_none():
    """Test when sample_size is None."""
    df = pd.DataFrame({"smiles": ["CCO", "CCN", "CCC"]})
    mock_logger = Mock()
    result = _apply_sampling(df, sample_size=None, logger=mock_logger)
    assert len(result) == 3


@pytest.mark.unit
def test_detect_mode_single_file(sample_csv_file):
    """Test detecting single file mode."""
    mode, paths = _detect_mode_and_paths(str(sample_csv_file))
    assert mode == "single_comparison"
    assert len(paths) == 1


@pytest.mark.unit
def test_detect_mode_multi_file_glob(tmp_path):
    """Test detecting multi-file mode with glob pattern."""
    # Create multiple files
    for i in range(3):
        file_path = tmp_path / f"model_{i}.csv"
        df = pd.DataFrame({"smiles": ["CCO", "CCN"], "model_name": [f"model_{i}"] * 2})
        df.to_csv(file_path, index=False)

    pattern = str(tmp_path / "model_*.csv")
    mode, paths = _detect_mode_and_paths(pattern)
    assert mode == "multi_comparison"
    assert len(paths) == 3


@pytest.mark.unit
def test_detect_mode_file_not_found():
    """Test when file pattern doesn't match any files."""
    with pytest.raises(FileNotFoundError):
        _detect_mode_and_paths("/nonexistent/path/*.csv")


@pytest.mark.unit
def test_load_single_comparison_data(sample_csv_file):
    """Test loading single comparison data."""
    mock_logger = Mock()
    result = _load_single_comparison_data(str(sample_csv_file), sample_size=None, logger=mock_logger)
    assert len(result) > 0
    assert "smiles" in result.columns
    assert "model_name" in result.columns


@pytest.mark.unit
def test_load_multi_comparison_data(tmp_path):
    """Test loading multi-comparison data from multiple files."""
    files = []
    for i in range(2):
        file_path = tmp_path / f"model_{i}.csv"
        df = pd.DataFrame({"smiles": ["CCO", "CCN"], "model_name": [f"model_{i}"] * 2})
        df.to_csv(file_path, index=False)
        files.append(str(file_path))

    mock_logger = Mock()
    result = _load_multi_comparison_data(files, sample_size=None, logger=mock_logger)
    assert len(result) == 4  # 2 rows per file * 2 files
    assert "smiles" in result.columns
    assert "model_name" in result.columns


@pytest.mark.unit
def test_prepare_input_data_single_mode(sample_csv_file, temp_output_dir):
    """Test preparing input data in single mode."""
    mock_logger = Mock()
    config = {
        "generated_mols_path": str(sample_csv_file),
        "folder_to_save": str(temp_output_dir),
        "save_sampled_mols": False,
    }
    result = prepare_input_data(config, mock_logger)
    assert len(result) > 0
    assert "smiles" in result.columns
    assert "model_name" in result.columns


@pytest.mark.unit
def test_prepare_input_data_removes_duplicates(tmp_path, temp_output_dir):
    """Test that prepare_input_data removes duplicates."""
    # Create CSV with duplicates
    csv_file = tmp_path / "duplicates.csv"
    df = pd.DataFrame({
        "smiles": ["CCO", "CCO", "CCN", "CCN"],
        "model_name": ["test", "test", "test", "test"],
    })
    df.to_csv(csv_file, index=False)

    mock_logger = Mock()
    config = {
        "generated_mols_path": str(csv_file),
        "folder_to_save": str(temp_output_dir),
        "save_sampled_mols": False,
    }
    result = prepare_input_data(config, mock_logger)
    # Should have 2 unique molecules
    assert len(result) == 2
