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


class _FakeProgress:
    """Minimal Progress replacement for testing CLI event rendering."""

    instances: list["_FakeProgress"] = []

    def __init__(self, *args, **kwargs):
        self.add_calls: list[dict] = []
        self.update_calls: list[dict] = []
        _FakeProgress.instances.append(self)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def add_task(self, description: str, **kwargs) -> int:
        task_id = len(self.add_calls) + 1
        self.add_calls.append({"task_id": task_id, "description": description, **kwargs})
        return task_id

    def update(self, task_id: int, **kwargs) -> None:
        self.update_calls.append({"task_id": task_id, **kwargs})


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
        """If base folder is available, keep it without suffix."""
        result = _get_unique_results_folder(tmp_path / "results")
        assert result == tmp_path / "results"

    def test_nonexistent_folder(self, tmp_path):
        """If base folder doesn't exist, keep the base name."""
        result = _get_unique_results_folder(tmp_path / "new_results")
        assert result == tmp_path / "new_results"

    def test_folder_exists_with_content(self, tmp_path):
        """If folder exists with content but no numbered folders, return _1."""
        results = tmp_path / "results"
        results.mkdir()
        (results / "file.txt").touch()

        result = _get_unique_results_folder(results)
        assert result == tmp_path / "results_1"

    def test_multiple_existing_folders(self, tmp_path):
        """If base folder is free, use it even when numbered folders exist."""
        # Create results_1 and results_2
        (tmp_path / "results_1").mkdir()
        (tmp_path / "results_2").mkdir()

        result = _get_unique_results_folder(tmp_path / "results")
        assert result == tmp_path / "results"

    def test_folder_already_numbered(self, tmp_path):
        """Base folder is still preferred when available."""
        # Create results_5
        (tmp_path / "results_5").mkdir()

        result = _get_unique_results_folder(tmp_path / "results")
        assert result == tmp_path / "results"


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

    def test_invalid_smiles_preserved_for_molprep_stage(self, tmp_path, mock_logger):
        """Invalid SMILES are preserved; MolPrep handles structural validation later."""
        input_file = tmp_path / "input.csv"
        input_file.write_text(
            "smiles,model_name\nCCO,test\ninvalid,test\nc1ccccc1,test"
        )

        output_folder = tmp_path / "output"
        result = preprocess_input_with_rdkit(
            str(input_file), output_folder, mock_logger
        )

        output_df = pd.read_csv(result)
        assert len(output_df) == 3


class TestStageEnum:
    """Tests for Stage enum."""

    def test_stage_values(self):
        """Test stage enum values."""
        assert Stage.mol_prep.value == "mol_prep"
        assert Stage.descriptors.value == "descriptors"
        assert Stage.struct_filters.value == "struct_filters"
        assert Stage.synthesis.value == "synthesis"
        assert Stage.docking.value == "docking"

    def test_stage_description(self):
        """Test stage descriptions."""
        assert "standard" in Stage.mol_prep.description.lower()
        assert "descriptors" in Stage.descriptors.description.lower()
        assert "filter" in Stage.struct_filters.description.lower()
        assert "synth" in Stage.synthesis.description.lower()
        assert "docking" in Stage.docking.description.lower()


def test_run_uses_single_progress_task_and_consistent_stage_numbers(tmp_path, monkeypatch):
    """CLI progress should reuse one task and keep stage numbering monotonic."""
    import hedgehog.main as main_mod

    _FakeProgress.instances.clear()

    results_dir = tmp_path / "results"
    input_path = tmp_path / "input.csv"
    input_path.write_text("smiles,model_name,mol_idx\nCCO,m1,0\n", encoding="utf-8")

    base_config = {
        "folder_to_save": str(results_dir),
        "generated_mols_path": str(input_path),
        "save_sampled_mols": False,
    }
    prepared_df = pd.DataFrame({"smiles": ["CCO"], "model_name": ["m1"], "mol_idx": [0]})

    monkeypatch.setattr(main_mod, "_display_banner", lambda: None)
    monkeypatch.setattr(main_mod, "_plain_output_enabled", lambda: False)
    monkeypatch.setattr(main_mod, "load_config", lambda *args, **kwargs: base_config.copy())
    monkeypatch.setattr(main_mod, "_apply_cli_overrides", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        main_mod,
        "_resolve_output_folder",
        lambda *args, **kwargs: Path(base_config["folder_to_save"]),
    )
    monkeypatch.setattr(
        main_mod.LoggerSingleton,
        "configure_log_directory",
        lambda self, folder_to_save: None,
    )
    monkeypatch.setattr(main_mod, "_preprocess_input", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        main_mod,
        "prepare_input_data",
        lambda *args, **kwargs: prepared_df.copy(),
    )
    monkeypatch.setattr(main_mod, "_save_sampled_molecules", lambda *args, **kwargs: None)
    monkeypatch.setattr(main_mod, "Progress", _FakeProgress)

    def _fake_calculate_metrics(data, config, progress_callback):
        progress_callback(
            {
                "type": "stage_start",
                "stage": "mol_prep",
                "stage_index": 1,
                "total_stages": 2,
                "message": "Mol Prep",
            }
        )
        progress_callback(
            {
                "type": "stage_progress",
                "stage": "mol_prep",
                "stage_index": 1,
                "total_stages": 2,
                "current": 10,
                "total": 10,
                "message": "Mol Prep",
            }
        )
        progress_callback(
            {
                "type": "stage_complete",
                "stage": "mol_prep",
                "stage_index": 1,
                "total_stages": 2,
                "message": "Mol Prep complete",
            }
        )
        progress_callback(
            {
                "type": "stage_start",
                "stage": "descriptors",
                "stage_index": 2,
                "total_stages": 2,
                "message": "Descriptors",
            }
        )
        progress_callback(
            {
                "type": "stage_progress",
                "stage": "descriptors",
                "stage_index": 2,
                "total_stages": 2,
                "current": 5,
                "total": 10,
                "message": "Descriptors",
            }
        )
        progress_callback(
            {
                "type": "stage_complete",
                "stage": "descriptors",
                "stage_index": 2,
                "total_stages": 2,
                "message": "Descriptors complete",
            }
        )
        return True

    monkeypatch.setattr(main_mod, "calculate_metrics", _fake_calculate_metrics)

    main_mod.run(
        config_path="unused.yml",
        generated_mols_path=None,
        out_dir=None,
        stage=None,
        reuse_folder=False,
        force_new_folder=False,
        auto_install=False,
    )

    progress_instance = _FakeProgress.instances[-1]
    assert len(progress_instance.add_calls) == 1

    descriptions = [
        call["description"]
        for call in progress_instance.update_calls
        if "description" in call
    ]
    assert any(desc.startswith("1/2 - Prep") for desc in descriptions)
    assert any(desc.startswith("2/2 - Descriptors") for desc in descriptions)
    assert not any(desc.startswith("0/2 - ") for desc in descriptions)
