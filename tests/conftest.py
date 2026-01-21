"""Common fixtures for HEDGEHOG tests."""

import logging
from pathlib import Path

import pandas as pd
import pytest
import yaml


@pytest.fixture
def sample_smiles():
    """10 simple valid SMILES strings."""
    return [
        "c1ccccc1",  # benzene
        "CC(C)C",  # isobutane
        "CCO",  # ethanol
        "CC=O",  # acetaldehyde
        "c1ccc(O)cc1",  # phenol
        "CC(=O)O",  # acetic acid
        "c1ccc(N)cc1",  # aniline
        "CCCC",  # butane
        "C1CCCCC1",  # cyclohexane
        "CC(C)(C)C",  # neopentane
    ]


@pytest.fixture
def sample_df(sample_smiles):
    """DataFrame with test molecules."""
    return pd.DataFrame(
        {
            "smiles": sample_smiles,
            "model_name": ["test"] * len(sample_smiles),
            "mol_idx": [f"test-{i}" for i in range(len(sample_smiles))],
        }
    )


@pytest.fixture
def multimodel_df():
    """DataFrame from data/test/generated_mols.csv."""
    test_data_path = (
        Path(__file__).parent.parent / "data" / "test" / "generated_mols.csv"
    )
    if test_data_path.exists():
        return pd.read_csv(test_data_path)
    # Fallback test data
    return pd.DataFrame(
        {
            "SMILES": [
                "c1cc2c(nc1N1CCN(C[C@@H]3CCCO3)CC1)OCO2",
                "COc1cccc(-c2nccnc2C2CN(C(=O)c3cc4ccccc4[nH]3)C2)c1",
                "O=C(NCCN1CCOCC1)c1ccc(-c2cnc3ccc(O)cc3c2)cc1",
            ],
            "model_name": ["test_model_a", "test_model_a", "test_model_b"],
        }
    )


@pytest.fixture
def sample_config(tmp_path):
    """Minimal config for tests."""
    return {
        "folder_to_save": str(tmp_path),
        "generated_mols_path": str(tmp_path / "generated_mols.csv"),
        "n_jobs": 1,
        "save_sampled_mols": False,
        "sample_size": None,
        "config_descriptors": str(tmp_path / "config_descriptors.yml"),
        "config_structFilters": str(tmp_path / "config_structFilters.yml"),
        "config_synthesis": str(tmp_path / "config_synthesis.yml"),
        "config_docking": str(tmp_path / "config_docking.yml"),
    }


@pytest.fixture
def mock_logger():
    """Mock logger for testing."""
    logger = logging.getLogger("test")
    logger.setLevel(logging.DEBUG)
    return logger


@pytest.fixture
def drug_like_smiles():
    """SMILES representing drug-like molecules for descriptor tests."""
    return [
        "CC(=O)Nc1ccc(O)cc1",  # paracetamol
        "CC(=O)Oc1ccccc1C(=O)O",  # aspirin
        "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",  # caffeine
        "c1ccccc1",  # benzene
        "CCO",  # ethanol
    ]


@pytest.fixture
def invalid_smiles():
    """List of invalid SMILES strings."""
    return [
        "invalid",
        "not_a_smiles",
        "???",
        "",
        "C(C)(C)(C)(C)C",  # invalid valence
    ]


@pytest.fixture
def minimal_pipeline_config(tmp_path):
    """Minimal config for pipeline integration tests."""
    config_names = [
        "config_descriptors",
        "config_structFilters",
        "config_synthesis",
        "config_docking",
    ]
    config_paths = {}
    disabled_config = yaml.dump({"run": False})

    for name in config_names:
        cfg_path = tmp_path / f"{name}.yml"
        cfg_path.write_text(disabled_config)
        config_paths[name] = str(cfg_path)

    return {
        "folder_to_save": str(tmp_path),
        "generated_mols_path": str(tmp_path / "generated_mols.csv"),
        "n_jobs": 1,
        "save_sampled_mols": False,
        "sample_size": None,
        **config_paths,
    }


@pytest.fixture
def stage_directories(tmp_path):
    """Create standard stage directory structure."""
    stages = tmp_path / "stages"
    stages.mkdir()

    dirs = {
        "descriptors": stages / "01_descriptors_initial",
        "struct_pre": stages / "02_structural_filters_pre",
        "struct_post": stages / "03_structural_filters_post",
        "synthesis": stages / "04_synthesis",
        "docking": stages / "05_docking",
    }
    for d in dirs.values():
        d.mkdir(parents=True)

    return dirs


@pytest.fixture
def test_csv_content():
    """Return a function to create test CSV content."""

    def _make_csv(smiles_list, model_names=None):
        if model_names is None:
            model_names = ["test"] * len(smiles_list)
        lines = ["smiles,model_name"]
        for smi, model in zip(smiles_list, model_names):
            lines.append(f"{smi},{model}")
        return "\n".join(lines)

    return _make_csv
