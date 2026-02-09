"""Tests for pipeline.py utilities."""

from pathlib import Path

from hedgehog.pipeline import (
    DIR_DESCRIPTORS_INITIAL,
    DIR_SYNTHESIS,
    DataChecker,
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
