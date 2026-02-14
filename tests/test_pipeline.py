"""Tests for pipeline.py utilities."""

from pathlib import Path

import pandas as pd
import pytest

from hedgehog.pipeline import (
    DIR_MOL_PREP,
    DIR_DESCRIPTORS_INITIAL,
    DIR_SYNTHESIS,
    DOCKING_SCORE_COLUMNS,
    DataChecker,
    MolecularAnalysisPipeline,
    PipelineStage,
    PipelineStageRunner,
    _directory_has_files,
    _file_exists_and_not_empty,
)


class TestFileExistsAndNotEmpty:
    """Tests for _file_exists_and_not_empty function."""

    def test_file_with_content(self, tmp_path):
        """File exists and has content - should return True."""
        f = tmp_path / "test.txt"
        f.write_text("content")
        assert _file_exists_and_not_empty(f) is True

    def test_empty_file(self, tmp_path):
        """File exists but is empty - should return False."""
        f = tmp_path / "empty.txt"
        f.touch()
        assert _file_exists_and_not_empty(f) is False

    def test_missing_file(self):
        """File does not exist - should return False."""
        assert _file_exists_and_not_empty(Path("/nonexistent/path/file.txt")) is False

    def test_with_string_path(self, tmp_path):
        """Function should also work with string paths converted to Path."""
        f = tmp_path / "test.txt"
        f.write_text("content")
        assert _file_exists_and_not_empty(f) is True

    def test_binary_file(self, tmp_path):
        """Binary file with content - should return True."""
        f = tmp_path / "binary.bin"
        f.write_bytes(b"\x00\x01\x02")
        assert _file_exists_and_not_empty(f) is True


class TestDirectoryHasFiles:
    """Tests for _directory_has_files function."""

    def test_directory_with_files(self, tmp_path):
        """Directory exists and contains files - should return True."""
        (tmp_path / "file.csv").touch()
        assert _directory_has_files(tmp_path) is True

    def test_empty_directory(self, tmp_path):
        """Directory exists but is empty - should return False."""
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()
        assert _directory_has_files(empty_dir) is False

    def test_missing_directory(self):
        """Directory does not exist - should return False."""
        assert _directory_has_files(Path("/nonexistent/directory")) is False

    def test_directory_with_subdirectory_only(self, tmp_path):
        """Directory contains only subdirectories, no files - should return False."""
        subdir = tmp_path / "subdir"
        subdir.mkdir()
        assert _directory_has_files(tmp_path) is False

    def test_directory_with_files_and_subdirs(self, tmp_path):
        """Directory contains both files and subdirectories - should return True."""
        (tmp_path / "subdir").mkdir()
        (tmp_path / "file.txt").touch()
        assert _directory_has_files(tmp_path) is True

    def test_directory_with_multiple_files(self, tmp_path):
        """Directory contains multiple files - should return True."""
        (tmp_path / "file1.csv").touch()
        (tmp_path / "file2.csv").touch()
        (tmp_path / "file3.csv").touch()
        assert _directory_has_files(tmp_path) is True


class TestDataChecker:
    """Tests for DataChecker class."""

    def test_check_stage_data_exists(self, tmp_path):
        """Check stage data when file exists with content."""
        # Create synthesis output
        synthesis_dir = tmp_path / "stages" / "04_synthesis"
        synthesis_dir.mkdir(parents=True)
        output_file = synthesis_dir / "filtered_molecules.csv"
        output_file.write_text("smiles,model_name\nCCO,test")

        config = {"folder_to_save": str(tmp_path)}
        checker = DataChecker(config)

        assert checker.check_stage_data(DIR_SYNTHESIS) is True

    def test_check_stage_data_missing(self, tmp_path):
        """Check stage data when file doesn't exist."""
        config = {"folder_to_save": str(tmp_path)}
        checker = DataChecker(config)

        assert checker.check_stage_data(DIR_SYNTHESIS) is False

    def test_check_stage_data_empty_file(self, tmp_path):
        """Check stage data when file exists but is empty."""
        synthesis_dir = tmp_path / "stages" / "04_synthesis"
        synthesis_dir.mkdir(parents=True)
        output_file = synthesis_dir / "filtered_molecules.csv"
        output_file.touch()  # Empty file

        config = {"folder_to_save": str(tmp_path)}
        checker = DataChecker(config)

        assert checker.check_stage_data(DIR_SYNTHESIS) is False

    def test_stage_has_molecules_header_only_csv(self, tmp_path):
        """Header-only stage CSV should be treated as having zero molecules."""
        synthesis_dir = tmp_path / "stages" / "04_synthesis"
        synthesis_dir.mkdir(parents=True)
        output_file = synthesis_dir / "filtered_molecules.csv"
        output_file.write_text("smiles,model_name,mol_idx\n")

        config = {"folder_to_save": str(tmp_path)}
        checker = DataChecker(config)

        assert checker.check_stage_data(DIR_SYNTHESIS) is True
        assert checker.stage_has_molecules(DIR_SYNTHESIS) is False

    def test_check_stage_data_legacy_structure(self, tmp_path):
        """Check stage data for legacy directory structure."""
        legacy_dir = tmp_path / "Synthesis"
        legacy_dir.mkdir()
        output_file = legacy_dir / "passSynthesisSMILES.csv"
        output_file.write_text("smiles\nCCO")

        config = {"folder_to_save": str(tmp_path)}
        checker = DataChecker(config)

        assert checker.check_stage_data("Synthesis") is True

    def test_check_stage_data_descriptors(self, tmp_path):
        """Check stage data for descriptors output."""
        desc_dir = tmp_path / "stages" / "01_descriptors_initial" / "filtered"
        desc_dir.mkdir(parents=True)
        output_file = desc_dir / "filtered_molecules.csv"
        output_file.write_text("smiles,model_name\nCCO,test")

        config = {"folder_to_save": str(tmp_path)}
        checker = DataChecker(config)

        assert checker.check_stage_data(DIR_DESCRIPTORS_INITIAL) is True

    def test_check_stage_data_mol_prep(self, tmp_path):
        """Check stage data for MolPrep output."""
        prep_dir = tmp_path / "stages" / "00_mol_prep"
        prep_dir.mkdir(parents=True)
        (prep_dir / "filtered_molecules.csv").write_text("smiles,model_name\nCCO,test")

        config = {"folder_to_save": str(tmp_path)}
        checker = DataChecker(config)

        assert checker.check_stage_data(DIR_MOL_PREP) is True


class TestPipelineStageRunner:
    """Tests for PipelineStageRunner class."""

    def test_find_latest_data_source_with_synthesis(self, tmp_path):
        """Find latest data source when synthesis output exists."""
        synthesis_dir = tmp_path / "stages" / "04_synthesis"
        synthesis_dir.mkdir(parents=True)
        (synthesis_dir / "filtered_molecules.csv").write_text("smiles\nCCO")

        config = {"folder_to_save": str(tmp_path)}
        checker = DataChecker(config)
        runner = PipelineStageRunner(config, checker)

        result = runner.find_latest_data_source()
        assert result == DIR_SYNTHESIS

    def test_find_latest_data_source_with_descriptors(self, tmp_path):
        """Find latest data source when only descriptors output exists."""
        desc_dir = tmp_path / "stages" / "01_descriptors_initial" / "filtered"
        desc_dir.mkdir(parents=True)
        (desc_dir / "filtered_molecules.csv").write_text("smiles\nCCO")

        config = {"folder_to_save": str(tmp_path)}
        checker = DataChecker(config)
        runner = PipelineStageRunner(config, checker)

        result = runner.find_latest_data_source()
        assert result == DIR_DESCRIPTORS_INITIAL

    def test_find_latest_data_source_with_mol_prep(self, tmp_path):
        """Find latest data source when only MolPrep output exists."""
        prep_dir = tmp_path / "stages" / "00_mol_prep"
        prep_dir.mkdir(parents=True)
        (prep_dir / "filtered_molecules.csv").write_text("smiles\nCCO")

        config = {"folder_to_save": str(tmp_path)}
        checker = DataChecker(config)
        runner = PipelineStageRunner(config, checker)

        result = runner.find_latest_data_source()
        assert result == DIR_MOL_PREP

    def test_find_latest_data_source_empty(self, tmp_path):
        """Find latest data source when no output exists."""
        config = {"folder_to_save": str(tmp_path)}
        checker = DataChecker(config)
        runner = PipelineStageRunner(config, checker)

        result = runner.find_latest_data_source()
        assert result is None

    def test_find_latest_data_source_priority(self, tmp_path):
        """Synthesis should be prioritized over descriptors."""
        # Create both synthesis and descriptors outputs
        synthesis_dir = tmp_path / "stages" / "04_synthesis"
        synthesis_dir.mkdir(parents=True)
        (synthesis_dir / "filtered_molecules.csv").write_text("smiles\nCCO")

        desc_dir = tmp_path / "stages" / "01_descriptors_initial" / "filtered"
        desc_dir.mkdir(parents=True)
        (desc_dir / "filtered_molecules.csv").write_text("smiles\nCCO")

        config = {"folder_to_save": str(tmp_path)}
        checker = DataChecker(config)
        runner = PipelineStageRunner(config, checker)

        result = runner.find_latest_data_source()
        assert result == DIR_SYNTHESIS

    def test_parse_docking_tools_string(self, tmp_path):
        """Parse docking tools from comma-separated string."""
        config = {"folder_to_save": str(tmp_path)}
        checker = DataChecker(config)
        runner = PipelineStageRunner(config, checker)

        result = runner._parse_docking_tools("smina, gnina")
        assert set(result) == {"smina", "gnina"}

    def test_parse_docking_tools_list(self, tmp_path):
        """Parse docking tools from list."""
        config = {"folder_to_save": str(tmp_path)}
        checker = DataChecker(config)
        runner = PipelineStageRunner(config, checker)

        result = runner._parse_docking_tools(["smina"])
        assert result == ["smina"]

    def test_parse_docking_tools_both(self, tmp_path):
        """Parse docking tools when 'both' specified."""
        config = {"folder_to_save": str(tmp_path)}
        checker = DataChecker(config)
        runner = PipelineStageRunner(config, checker)

        result = runner._parse_docking_tools("both")
        assert set(result) == {"smina", "gnina"}


class TestPipelineStage:
    """Tests for PipelineStage class."""

    def test_stage_initialization(self):
        """Test PipelineStage initialization."""
        stage = PipelineStage("test_stage", "config_test", "test/dir")

        assert stage.name == "test_stage"
        assert stage.config_key == "config_test"
        assert stage.directory == "test/dir"
        assert stage.enabled is False
        assert stage.completed is False

    def test_stage_enabled_flag(self):
        """Test modifying stage enabled flag."""
        stage = PipelineStage("test", "config", "dir")
        stage.enabled = True

        assert stage.enabled is True


def _write_test_docking_sdf(path: Path, records: list[dict]) -> None:
    """Write a minimal docking SDF with specified properties."""
    pytest.importorskip("rdkit")
    from rdkit import Chem

    path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(path))
    for rec in records:
        mol = Chem.MolFromSmiles("CC")
        assert mol is not None
        mol.SetProp("_Name", str(rec["name"]))
        for key, value in rec.items():
            if key == "name":
                continue
            mol.SetProp(str(key), str(value))
        writer.write(mol)
    writer.close()


class TestFinalOutputWithDockingScores:
    """Tests for final output enrichment with docking scores."""

    def _pipeline(self, tmp_path: Path) -> MolecularAnalysisPipeline:
        return MolecularAnalysisPipeline({"folder_to_save": str(tmp_path)})

    def test_preserves_existing_columns_and_adds_score_columns(self, tmp_path):
        """Keep all existing columns and append docking score schema."""
        pipeline = self._pipeline(tmp_path)
        final_data = pd.DataFrame(
            {
                "smiles": ["CCO", "CCC"],
                "model_name": ["model_a", "model_b"],
                "mol_idx": ["m1", "m2"],
                "sa_score": [2.1, 3.4],
                "custom_flag": [True, False],
            }
        )

        pipeline._save_final_output(final_data, len(final_data))

        output_path = tmp_path / "output" / "final_molecules.csv"
        out = pd.read_csv(output_path)
        assert "sa_score" in out.columns
        assert "custom_flag" in out.columns
        for col in DOCKING_SCORE_COLUMNS:
            assert col in out.columns
            assert out[col].isna().all()

    def test_adds_gnina_and_smina_scores_using_best_pose(self, tmp_path):
        """For each tool, choose the best pose (minimum affinity) per molecule."""
        pipeline = self._pipeline(tmp_path)

        _write_test_docking_sdf(
            tmp_path / "stages" / "05_docking" / "gnina" / "gnina_out.sdf",
            [
                {
                    "name": "mol-1",
                    "minimizedAffinity": -7.1,
                    "CNNscore": 0.51,
                    "CNNaffinity": -6.0,
                    "CNN_VS": 0.22,
                },
                {
                    "name": "mol-1",
                    "minimizedAffinity": -8.4,
                    "CNNscore": 0.73,
                    "CNNaffinity": -6.8,
                    "CNN_VS": 0.31,
                },
                {
                    "name": "mol-2",
                    "minimizedAffinity": -6.3,
                    "CNNscore": 0.44,
                    "CNNaffinity": -5.1,
                    "CNN_VS": 0.11,
                },
            ],
        )
        _write_test_docking_sdf(
            tmp_path / "stages" / "05_docking" / "smina" / "smina_out.sdf",
            [
                {"name": "mol-1", "minimizedAffinity": -6.0},
                {"name": "mol-1", "minimizedAffinity": -9.2},
                {"name": "mol-2", "minimizedAffinity": -5.7},
            ],
        )

        final_data = pd.DataFrame(
            {
                "smiles": ["CCO", "CCC", "CCN"],
                "model_name": ["model_a", "model_a", "model_b"],
                "mol_idx": ["mol-1", "mol-2", "mol-3"],
            }
        )
        pipeline._save_final_output(final_data, len(final_data))

        output_path = tmp_path / "output" / "final_molecules.csv"
        out = pd.read_csv(output_path).assign(
            mol_idx=lambda df: df["mol_idx"].astype(str)
        )
        rows = out.set_index("mol_idx")

        assert rows.loc["mol-1", "gnina_affinity"] == pytest.approx(-8.4)
        assert rows.loc["mol-1", "gnina_cnnscore"] == pytest.approx(0.73)
        assert rows.loc["mol-1", "gnina_cnnaffinity"] == pytest.approx(-6.8)
        assert rows.loc["mol-1", "gnina_cnn_vs"] == pytest.approx(0.31)
        assert rows.loc["mol-1", "smina_affinity"] == pytest.approx(-9.2)

        assert rows.loc["mol-2", "gnina_affinity"] == pytest.approx(-6.3)
        assert rows.loc["mol-2", "smina_affinity"] == pytest.approx(-5.7)

        assert pd.isna(rows.loc["mol-3", "gnina_affinity"])
        assert pd.isna(rows.loc["mol-3", "smina_affinity"])

    def test_empty_final_output_has_stable_header(self, tmp_path):
        """Empty final output should still expose a stable score schema."""
        pipeline = self._pipeline(tmp_path)
        empty_data = pd.DataFrame(columns=["smiles", "model_name", "mol_idx"])

        pipeline._save_final_output(empty_data, 0)

        output_path = tmp_path / "output" / "final_molecules.csv"
        out = pd.read_csv(output_path)
        assert list(out.columns) == [
            "smiles",
            "model_name",
            "mol_idx",
            *DOCKING_SCORE_COLUMNS,
        ]
