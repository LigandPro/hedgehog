"""Shared fixtures for pytest tests."""

import os
import shutil
from pathlib import Path

import pandas as pd
import pytest


@pytest.fixture
def fixtures_dir():
    """Return path to fixtures directory."""
    return Path(__file__).parent / "fixtures"


@pytest.fixture
def test_molecules_csv(fixtures_dir):
    """Return path to test molecules CSV file."""
    return fixtures_dir / "test_molecules.csv"


@pytest.fixture
def test_config_path(fixtures_dir):
    """Return path to test config YAML file."""
    return fixtures_dir / "test_config.yml"


@pytest.fixture
def temp_output_dir(tmp_path):
    """Create temporary output directory for tests."""
    output_dir = tmp_path / "test_output"
    output_dir.mkdir(parents=True, exist_ok=True)
    yield output_dir
    # Cleanup after test
    if output_dir.exists():
        shutil.rmtree(output_dir)


@pytest.fixture
def sample_smiles_data():
    """Return sample SMILES data as pandas DataFrame."""
    data = {
        "smiles": [
            "CCO",
            "CC(C)O",
            "c1ccccc1",
            "CC(=O)O",
            "CCN",
        ],
        "model_name": [
            "test_model",
            "test_model",
            "test_model",
            "test_model",
            "test_model",
        ],
    }
    return pd.DataFrame(data)


@pytest.fixture
def sample_multi_model_data():
    """Return sample data with multiple models."""
    data = {
        "smiles": [
            "CCO",
            "CC(C)O",
            "c1ccccc1",
            "CC(=O)O",
            "CCN",
            "c1ccc(cc1)O",
        ],
        "model_name": [
            "model_a",
            "model_a",
            "model_a",
            "model_b",
            "model_b",
            "model_b",
        ],
    }
    return pd.DataFrame(data)


@pytest.fixture
def sample_csv_file(tmp_path, sample_smiles_data):
    """Create temporary CSV file with sample data."""
    csv_path = tmp_path / "sample.csv"
    sample_smiles_data.to_csv(csv_path, index=False)
    return csv_path


@pytest.fixture
def sample_multi_model_csv(tmp_path, sample_multi_model_data):
    """Create temporary CSV file with multi-model data."""
    csv_path = tmp_path / "multi_model.csv"
    sample_multi_model_data.to_csv(csv_path, index=False)
    return csv_path


@pytest.fixture
def sample_config_dict(temp_output_dir, test_molecules_csv, fixtures_dir):
    """Return sample configuration dictionary."""
    return {
        "generated_mols_path": str(test_molecules_csv),
        "folder_to_save": str(temp_output_dir),
        "n_jobs": 2,
        "sample_size": 10,
        "save_sampled_mols": True,
        "config_descriptors": str(fixtures_dir / "test_config_descriptors.yml"),
        "config_structFilters": str(fixtures_dir / "test_config_structFilters.yml"),
        "config_synthesis": str(fixtures_dir / "test_config_synthesis.yml"),
        "config_docking": str(fixtures_dir / "test_config_docking.yml"),
    }


@pytest.fixture(autouse=True)
def cleanup_test_outputs():
    """Cleanup test output directories after each test."""
    yield
    # Cleanup common test output directories
    test_dirs = ["tests/test_output", "results/test"]
    for dir_path in test_dirs:
        if os.path.exists(dir_path):
            shutil.rmtree(dir_path)
