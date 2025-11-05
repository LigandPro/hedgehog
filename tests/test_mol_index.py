"""Unit tests for molecular index assignment module."""

import json

import pandas as pd
import pytest

from hedge.utils.mol_index import assign_mol_idx


@pytest.mark.unit
def test_assign_mol_idx_format(sample_smiles_data, temp_output_dir):
    """Test that mol_idx has correct format LP-XXXX-XXXXX."""
    result = assign_mol_idx(sample_smiles_data, temp_output_dir)
    assert "mol_idx" in result.columns

    for idx in result["mol_idx"]:
        assert idx.startswith("LP-")
        parts = idx.split("-")
        assert len(parts) == 3
        assert len(parts[1]) == 4  # Model number is 4 digits
        assert len(parts[2]) == 5  # Counter is 5 digits


@pytest.mark.unit
def test_assign_mol_idx_single_model(sample_smiles_data, temp_output_dir):
    """Test mol_idx assignment for single model."""
    result = assign_mol_idx(sample_smiles_data, temp_output_dir)

    # All should have same model number (0001)
    model_numbers = {idx.split("-")[1] for idx in result["mol_idx"]}
    assert len(model_numbers) == 1
    assert "0001" in model_numbers

    # Counters should be sequential
    counters = [int(idx.split("-")[2]) for idx in result["mol_idx"]]
    assert counters == list(range(1, len(result) + 1))


@pytest.mark.unit
def test_assign_mol_idx_multiple_models(sample_multi_model_data, temp_output_dir):
    """Test mol_idx assignment for multiple models."""
    result = assign_mol_idx(sample_multi_model_data, temp_output_dir)

    # Should have different model numbers for different models
    model_a = result[result["model_name"] == "model_a"]
    model_b = result[result["model_name"] == "model_b"]

    model_a_numbers = {idx.split("-")[1] for idx in model_a["mol_idx"]}
    model_b_numbers = {idx.split("-")[1] for idx in model_b["mol_idx"]}

    assert len(model_a_numbers) == 1
    assert len(model_b_numbers) == 1
    assert model_a_numbers != model_b_numbers


@pytest.mark.unit
def test_assign_mol_idx_uniqueness(sample_multi_model_data, temp_output_dir):
    """Test that all mol_idx values are unique."""
    result = assign_mol_idx(sample_multi_model_data, temp_output_dir)
    assert result["mol_idx"].nunique() == len(result)


@pytest.mark.unit
def test_assign_mol_idx_preserves_data(sample_smiles_data, temp_output_dir):
    """Test that assign_mol_idx preserves original data."""
    original_smiles = sample_smiles_data["smiles"].tolist()
    result = assign_mol_idx(sample_smiles_data, temp_output_dir)

    assert result["smiles"].tolist() == original_smiles
    assert len(result) == len(sample_smiles_data)


@pytest.mark.unit
def test_assign_mol_idx_missing_model_name(temp_output_dir):
    """Test mol_idx assignment when model_name is missing."""
    df = pd.DataFrame({"smiles": ["CCO", "CCN", "CCC"]})
    result = assign_mol_idx(df, temp_output_dir)

    assert "model_name" in result.columns
    assert "mol_idx" in result.columns
    # Should add default model name
    assert all(result["model_name"] == "single")


@pytest.mark.unit
def test_assign_mol_idx_saves_mapping(sample_smiles_data, temp_output_dir):
    """Test that model index mapping is saved to JSON."""
    assign_mol_idx(sample_smiles_data, temp_output_dir)

    mapping_file = temp_output_dir / "run_configs" / "model_index_map.json"
    assert mapping_file.exists()

    with open(mapping_file) as f:
        mapping = json.load(f)

    assert "test_model" in mapping
    assert mapping["test_model"] == 1


@pytest.mark.unit
def test_assign_mol_idx_per_model_counter(sample_multi_model_data, temp_output_dir):
    """Test that counter is per-model and starts from 1."""
    result = assign_mol_idx(sample_multi_model_data, temp_output_dir)

    for model in result["model_name"].unique():
        model_data = result[result["model_name"] == model]
        counters = [int(idx.split("-")[2]) for idx in model_data["mol_idx"]]
        # Counters should start from 1 and be sequential
        assert counters == list(range(1, len(model_data) + 1))
