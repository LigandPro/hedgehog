"""Unit tests for CLI module."""

from pathlib import Path
from unittest.mock import Mock, patch

import pandas as pd
import pytest

from hedge.main import (
    _validate_input_path,
    preprocess_input_with_rdkit,
    preprocess_input_with_tool,
)


@pytest.mark.unit
def test_validate_input_path_valid(sample_csv_file):
    """Test validating a valid input path."""
    result = _validate_input_path(str(sample_csv_file))
    assert result is not None
    assert isinstance(result, Path)
    assert result.exists()


@pytest.mark.unit
def test_validate_input_path_with_glob():
    """Test validating a path with glob pattern."""
    result = _validate_input_path("data/*.csv")
    assert result is None


@pytest.mark.unit
def test_validate_input_path_not_exists():
    """Test validating a non-existent path."""
    result = _validate_input_path("/nonexistent/file.csv")
    assert result is None


@pytest.mark.unit
def test_preprocess_input_with_rdkit_success(tmp_path):
    """Test successful preprocessing with RDKit."""
    # Create input file
    input_file = tmp_path / "input.csv"
    df = pd.DataFrame({
        "smiles": ["CCO", "c1ccccc1", "CCN"],
        "model_name": ["test", "test", "test"],
    })
    df.to_csv(input_file, index=False)

    mock_logger = Mock()
    result = preprocess_input_with_rdkit(str(input_file), tmp_path, mock_logger)

    assert result is not None
    assert Path(result).exists()
    # Check output file
    output_df = pd.read_csv(result)
    assert len(output_df) == 3
    assert "smiles" in output_df.columns


@pytest.mark.unit
def test_preprocess_input_with_rdkit_invalid_smiles(tmp_path):
    """Test preprocessing with invalid SMILES."""
    # Create input file with some invalid SMILES
    input_file = tmp_path / "input.csv"
    df = pd.DataFrame({
        "smiles": ["CCO", "INVALID_SMILES", "CCN"],
        "model_name": ["test", "test", "test"],
    })
    df.to_csv(input_file, index=False)

    mock_logger = Mock()
    result = preprocess_input_with_rdkit(str(input_file), tmp_path, mock_logger)

    assert result is not None
    # Invalid SMILES should be filtered out
    output_df = pd.read_csv(result)
    assert len(output_df) == 2
    assert "INVALID_SMILES" not in output_df["smiles"].values


@pytest.mark.unit
def test_preprocess_input_with_rdkit_removes_duplicates(tmp_path):
    """Test that preprocessing removes duplicates."""
    # Create input file with duplicates
    input_file = tmp_path / "input.csv"
    df = pd.DataFrame({
        "smiles": ["CCO", "CCO", "CCN", "CCN"],
        "model_name": ["test", "test", "test", "test"],
    })
    df.to_csv(input_file, index=False)

    mock_logger = Mock()
    result = preprocess_input_with_rdkit(str(input_file), tmp_path, mock_logger)

    assert result is not None
    output_df = pd.read_csv(result)
    # Should have only 2 unique molecules
    assert len(output_df) == 2
    mock_logger.info.assert_called()


@pytest.mark.unit
def test_preprocess_input_with_rdkit_glob_pattern(tmp_path):
    """Test preprocessing with glob pattern returns None."""
    mock_logger = Mock()
    result = preprocess_input_with_rdkit("data/*.csv", tmp_path, mock_logger)
    assert result is None


@pytest.mark.unit
def test_preprocess_input_with_tool_csv(tmp_path):
    """Test preprocessing with external tool for CSV file."""
    input_file = tmp_path / "input.csv"
    df = pd.DataFrame({"smiles": ["CCO", "CCN"], "model_name": ["test", "test"]})
    df.to_csv(input_file, index=False)

    mock_logger = Mock()
    with patch("subprocess.run") as mock_run:
        mock_run.return_value = Mock(returncode=0)
        # Create expected output file
        output_file = tmp_path / "prepared_input.csv"
        df.to_csv(output_file, index=False)

        result = preprocess_input_with_tool(
            str(input_file),
            "fake_tool",
            tmp_path,
            mock_logger,
        )

        # Check that result is returned and subprocess was called
        assert result is not None
        mock_run.assert_called_once()
        args = mock_run.call_args[0][0]
        assert "fake_tool" in args
        assert "-icsv" in args


@pytest.mark.unit
def test_preprocess_input_with_tool_smi(tmp_path):
    """Test preprocessing with external tool for SMI file."""
    input_file = tmp_path / "input.smi"
    input_file.write_text("CCO\nCCN\n")

    mock_logger = Mock()
    with patch("subprocess.run") as mock_run:
        mock_run.return_value = Mock(returncode=0)
        output_file = tmp_path / "prepared_input.csv"
        pd.DataFrame({"smiles": ["CCO", "CCN"]}).to_csv(output_file, index=False)

        result = preprocess_input_with_tool(
            str(input_file),
            "fake_tool",
            tmp_path,
            mock_logger,
        )

        # Check that result is returned and subprocess was called with -ismi
        assert result is not None
        args = mock_run.call_args[0][0]
        assert "-ismi" in args


@pytest.mark.unit
def test_preprocess_input_with_tool_unsupported_format(tmp_path):
    """Test preprocessing with unsupported file format."""
    input_file = tmp_path / "input.xyz"
    input_file.write_text("some data")

    mock_logger = Mock()
    result = preprocess_input_with_tool(str(input_file), "fake_tool", tmp_path, mock_logger)
    assert result is None


@pytest.mark.unit
def test_preprocess_input_with_tool_fails(tmp_path):
    """Test preprocessing when external tool fails."""
    input_file = tmp_path / "input.csv"
    df = pd.DataFrame({"smiles": ["CCO"]})
    df.to_csv(input_file, index=False)

    mock_logger = Mock()
    with patch("subprocess.run") as mock_run:
        mock_run.side_effect = Exception("Tool failed")

        result = preprocess_input_with_tool(
            str(input_file),
            "fake_tool",
            tmp_path,
            mock_logger,
        )
        assert result is None
