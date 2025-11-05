"""Unit tests for pipeline module."""

from pathlib import Path
from unittest.mock import Mock, patch

import pandas as pd
import pytest

from hedge.pipeline import (
    DataChecker,
    MolecularAnalysisPipeline,
    PipelineStage,
    calculate_metrics,
)


@pytest.mark.unit
def test_pipeline_stage_initialization():
    """Test PipelineStage initialization."""
    stage = PipelineStage("test_stage", "config_test", "test_dir")
    assert stage.name == "test_stage"
    assert stage.config_key == "config_test"
    assert stage.directory == "test_dir"
    assert stage.enabled is False
    assert stage.completed is False


@pytest.mark.unit
def test_data_checker_initialization(sample_config_dict):
    """Test DataChecker initialization."""
    checker = DataChecker(sample_config_dict)
    assert checker.config == sample_config_dict
    assert isinstance(checker.base_path, Path)


@pytest.mark.unit
def test_data_checker_stage_data_not_exists(sample_config_dict):
    """Test checking for non-existent stage data."""
    checker = DataChecker(sample_config_dict)
    result = checker.check_stage_data("Descriptors")
    assert result is False


@pytest.mark.unit
def test_data_checker_stage_data_exists(sample_config_dict, temp_output_dir):
    """Test checking for existing stage data."""
    # Create output file
    descriptors_dir = temp_output_dir / "Descriptors"
    descriptors_dir.mkdir(parents=True, exist_ok=True)
    output_file = descriptors_dir / "passDescriptorsSMILES.csv"
    pd.DataFrame({"smiles": ["CCO"], "model_name": ["test"]}).to_csv(output_file, index=False)

    checker = DataChecker(sample_config_dict)
    result = checker.check_stage_data("Descriptors")
    assert result is True


@pytest.mark.unit
def test_molecular_analysis_pipeline_initialization(sample_config_dict):
    """Test MolecularAnalysisPipeline initialization."""
    pipeline = MolecularAnalysisPipeline(sample_config_dict)
    assert pipeline.config == sample_config_dict
    assert isinstance(pipeline.data_checker, DataChecker)
    assert len(pipeline.stages) == 6  # 6 stages in total


@pytest.mark.unit
def test_pipeline_get_latest_data_no_data(sample_config_dict):
    """Test getting latest data when no data exists."""
    pipeline = MolecularAnalysisPipeline(sample_config_dict)
    pipeline.current_data = pd.DataFrame({"smiles": ["CCO"], "model_name": ["test"]})
    result = pipeline.get_latest_data()
    # Should return current_data when no stage data exists
    assert result is not None
    assert len(result) == 1


@pytest.mark.unit
def test_pipeline_get_latest_data_with_data(sample_config_dict, temp_output_dir):
    """Test getting latest data when stage data exists."""
    # Create output file
    descriptors_dir = temp_output_dir / "Descriptors"
    descriptors_dir.mkdir(parents=True, exist_ok=True)
    output_file = descriptors_dir / "passDescriptorsSMILES.csv"
    test_data = pd.DataFrame({"smiles": ["CCO", "CCN"], "model_name": ["test", "test"]})
    test_data.to_csv(output_file, index=False)

    pipeline = MolecularAnalysisPipeline(sample_config_dict)
    result = pipeline.get_latest_data()
    assert result is not None
    assert len(result) == 2


@pytest.mark.unit
def test_calculate_metrics_with_mock(sample_config_dict, sample_smiles_data):
    """Test calculate_metrics function with mocked pipeline."""
    with patch("hedge.pipeline.MolecularAnalysisPipeline") as mock_pipeline_class:
        mock_pipeline = Mock()
        mock_pipeline.run_pipeline.return_value = True
        mock_pipeline_class.return_value = mock_pipeline

        result = calculate_metrics(sample_smiles_data, sample_config_dict)
        assert result is True
        mock_pipeline.run_pipeline.assert_called_once_with(sample_smiles_data)


@pytest.mark.unit
def test_calculate_metrics_error_handling(sample_config_dict, sample_smiles_data):
    """Test calculate_metrics error handling."""
    with patch("hedge.pipeline.MolecularAnalysisPipeline") as mock_pipeline_class:
        mock_pipeline_class.side_effect = Exception("Test error")

        result = calculate_metrics(sample_smiles_data, sample_config_dict)
        assert result is False


@pytest.mark.integration
def test_pipeline_stages_disabled_by_default(sample_config_dict):
    """Test that certain stages are disabled in test config."""
    pipeline = MolecularAnalysisPipeline(sample_config_dict)

    # Check that struct_filters is disabled (as per test config)
    struct_filter_stages = [s for s in pipeline.stages if "struct" in s.name.lower()]
    # At least one should be disabled
    disabled_count = sum(1 for s in struct_filter_stages if not s.enabled)
    assert disabled_count > 0


@pytest.mark.integration
def test_pipeline_descriptors_enabled_by_default(sample_config_dict):
    """Test that descriptors stage is enabled in test config."""
    pipeline = MolecularAnalysisPipeline(sample_config_dict)

    descriptor_stages = [s for s in pipeline.stages if "descriptor" in s.name.lower()]
    # At least one descriptors stage should be enabled
    enabled_count = sum(1 for s in descriptor_stages if s.enabled)
    assert enabled_count > 0
