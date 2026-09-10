"""Tests for main.py utilities."""

import os
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
import yaml

from hedgehog.main import (
    Stage,
    _folder_is_empty,
    _get_input_format_flag,
    _get_unique_results_folder,
    _resolve_config_paths,
    _resolve_output_folder,
    _save_sampled_molecules,
    _validate_input_path,
    preprocess_input_with_rdkit,
)


class _FakeProgress:
    """Minimal Progress replacement for testing CLI event rendering."""

    instances: list["_FakeProgress"] = []

    def __init__(self, *args, **kwargs):
        self.add_calls: list[dict] = []
        self.reset_calls: list[dict] = []
        self.update_calls: list[dict] = []
        _FakeProgress.instances.append(self)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def add_task(self, description: str, **kwargs) -> int:
        task_id = len(self.add_calls) + 1
        self.add_calls.append(
            {"task_id": task_id, "description": description, **kwargs}
        )
        return task_id

    def reset(self, task_id: int, **kwargs) -> None:
        self.reset_calls.append({"task_id": task_id, **kwargs})

    def update(self, task_id: int, **kwargs) -> None:
        self.update_calls.append({"task_id": task_id, **kwargs})


def test_resolve_config_paths_relative_to_config_file(tmp_path):
    config_dir = tmp_path / "configs"
    config_dir.mkdir()
    input_file = config_dir / "mols.csv"
    input_file.write_text("smiles\nCCO\n", encoding="utf-8")
    sub_config = config_dir / "config_descriptors.yml"
    sub_config.write_text("run: true\n", encoding="utf-8")
    config_file = config_dir / "config.yml"
    config_file.write_text("", encoding="utf-8")

    config = {
        "generated_mols_path": "mols.csv",
        "config_descriptors": "config_descriptors.yml",
        "folder_to_save": "results/run",
    }

    _resolve_config_paths(config, str(config_file))

    assert config["generated_mols_path"] == str(input_file.resolve())
    assert config["config_descriptors"] == str(sub_config.resolve())
    assert config["folder_to_save"] == "results/run"


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

    def test_first_run_is_numbered_one(self, tmp_path):
        """Fresh installs always get base_name_1."""
        result = _get_unique_results_folder(tmp_path / "results")
        assert result == tmp_path / "results_1"

    def test_nonexistent_parent_still_numbers_from_one(self, tmp_path):
        """Numbering starts at 1 even when the parent directory is new."""
        result = _get_unique_results_folder(tmp_path / "new_results")
        assert result == tmp_path / "new_results_1"

    def test_unnumbered_base_folder_advances_sequence(self, tmp_path):
        """Unnumbered base folder with content advances the run suffix."""
        results = tmp_path / "results"
        results.mkdir()
        (results / "file.txt").touch()

        result = _get_unique_results_folder(results)
        assert result == tmp_path / "results_2"

    def test_increments_past_existing_numbered_folders(self, tmp_path):
        """Next folder is one above the highest existing suffix."""
        (tmp_path / "results_1").mkdir()
        (tmp_path / "results_2").mkdir()

        result = _get_unique_results_folder(tmp_path / "results")
        assert result == tmp_path / "results_3"

    def test_skips_gaps_in_numbering(self, tmp_path):
        """Highest suffix wins even when intermediate numbers are missing."""
        (tmp_path / "results_5").mkdir()

        result = _get_unique_results_folder(tmp_path / "results")
        assert result == tmp_path / "results_6"


class TestResolveOutputFolder:
    """Tests for explicit fresh-run and reuse folder behavior."""

    def test_stage_selection_creates_fresh_folder_by_default(self, tmp_path):
        """Selecting one stage must not silently reuse an existing run."""
        base = tmp_path / "run"
        base.mkdir()
        (base / "RUN_INFO.md").write_text("existing run\n", encoding="utf-8")

        result = _resolve_output_folder(
            {"folder_to_save": str(base)},
            reuse_folder=False,
            force_new_folder=False,
        )

        assert result == tmp_path / "run_2"

    def test_reuse_is_the_only_existing_folder_mode(self, tmp_path):
        """The configured directory is returned only with explicit reuse."""
        base = tmp_path / "run"
        base.mkdir()

        result = _resolve_output_folder(
            {"folder_to_save": str(base)},
            reuse_folder=True,
            force_new_folder=False,
        )

        assert result == base


def test_continue_command_loads_saved_input_and_skips_completed_stages(
    tmp_path, monkeypatch
):
    from hedgehog import main as main_mod

    run_dir = tmp_path / "unfinished"
    configs_dir = run_dir / "configs"
    input_dir = run_dir / "input"
    configs_dir.mkdir(parents=True)
    input_dir.mkdir()
    (run_dir / ".RUN_INCOMPLETE").write_text("interrupted\n", encoding="utf-8")
    (input_dir / "sampled_molecules.csv").write_text(
        "smiles,model_name,mol_idx\nCCO,target,0\n", encoding="utf-8"
    )

    master = {"folder_to_save": str(run_dir), "generated_mols_path": "unused.csv"}
    for key in (
        "config_mol_prep",
        "config_descriptors",
        "config_structFilters",
        "config_synthesis",
        "config_docking",
        "config_docking_filters",
    ):
        config_path = configs_dir / f"{key}.yml"
        config_path.write_text("run: true\n", encoding="utf-8")
        master[key] = str(config_path)
    (configs_dir / "master_config_resolved.yml").write_text(
        yaml.safe_dump(master, sort_keys=False), encoding="utf-8"
    )

    for directory in (
        "stages/01_mol_prep",
        "stages/02_descriptors_initial",
        "stages/03_structural_filters_post",
    ):
        stage_dir = run_dir / directory
        stage_dir.mkdir(parents=True)
        (stage_dir / "filtered_molecules.csv").write_text(
            "smiles,model_name,mol_idx\nCCO,target,0\n", encoding="utf-8"
        )

    captured = {}

    def fake_calculate(data, config, progress_callback):
        captured["data"] = data
        captured["config"] = config
        captured["progress_callback"] = progress_callback
        return True

    monkeypatch.setattr(main_mod, "calculate_metrics", fake_calculate)
    main_mod._run_pipeline_command(
        config_path=None,
        generated_mols_path=None,
        out_dir=None,
        stage=None,
        reuse_folder=False,
        force_new_folder=False,
        auto_install=False,
        show_progress=False,
        large_dataset=False,
        continue_folder=str(run_dir),
    )

    assert len(captured["data"]) == 1
    assert captured["config"]["folder_to_save"] == str(run_dir.resolve())
    assert captured["config"]["_continue_mode"] is True
    assert captured["config"]["_continue_completed_stages"] == [
        "mol_prep",
        "descriptors",
        "struct_filters",
    ]
    assert captured["config"]["_run_stage_selection_override"] == [
        "synthesis",
        "docking",
        "docking_filters",
        "final_descriptors",
    ]


def test_target_alignment_preserves_cli_stage_selection_in_probe_and_candidates(
    tmp_path, monkeypatch
):
    from hedgehog import main as main_mod

    target = tmp_path / "targets.csv"
    target.write_text("smiles,mol_idx\nCCO,target-1\n", encoding="utf-8")
    target_run = tmp_path / "outer" / "target_alignment" / "calibration_target_run"
    stage_selection = ["mol_prep", "descriptors"]
    config = {
        "generated_mols_path": str(target),
        "target_mols_path": str(target),
        "folder_to_save": str(tmp_path / "outer"),
        "save_sampled_mols": True,
        "_run_stage_selection_override": stage_selection,
    }
    captured = {}

    def fake_create_probe(master, _target_path, _alignment_root):
        probe = dict(master)
        probe["folder_to_save"] = str(target_run)
        captured["probe_created"] = probe
        return probe

    def fake_calculate(data, probe, _callback):
        captured["probe_run"] = dict(probe)
        return True

    def fake_finalize(*_args, **_kwargs):
        return (
            {"folder_to_save": str(tmp_path / "discarded")},
            tmp_path / "aligned.yml",
            tmp_path / "thresholds.yml",
        )

    monkeypatch.setattr(main_mod, "create_probe_config", fake_create_probe)
    monkeypatch.setattr(main_mod, "_preprocess_input", lambda *_args: None)
    monkeypatch.setattr(
        main_mod,
        "prepare_input_data",
        lambda *_args: pd.DataFrame({"smiles": ["CCO"], "mol_idx": ["target-1"]}),
    )
    monkeypatch.setattr(
        main_mod, "set_probe_molprep_allowed_atoms", lambda *_args: None
    )
    monkeypatch.setattr(main_mod, "calculate_metrics", fake_calculate)
    monkeypatch.setattr(main_mod, "finalize_global_alignment", fake_finalize)

    aligned = main_mod._align_config_with_target_molecules(
        config, tmp_path / "outer", 95
    )

    assert captured["probe_created"]["_run_stage_selection_override"] == stage_selection
    assert captured["probe_run"]["_run_stage_selection_override"] == stage_selection
    assert aligned["_run_stage_selection_override"] == stage_selection


def test_continue_nested_alignment_then_starts_candidates(tmp_path, monkeypatch):
    """An interrupted target probe should finish alignment and resume its outer run."""
    from hedgehog import main as main_mod

    outer = tmp_path / "run"
    target_run = outer / "target_alignment" / "target_run"
    target_configs = target_run / "configs"
    target_configs.mkdir(parents=True)
    (target_run / ".RUN_INCOMPLETE").write_text("unfinished\n", encoding="utf-8")
    (target_run / "input").mkdir()
    (target_run / "input" / "sampled_molecules.csv").write_text(
        "smiles,model_name,mol_idx\nCCO,target,0\n", encoding="utf-8"
    )
    target_csv = tmp_path / "targets.csv"
    target_csv.write_text("smiles\nCCO\n", encoding="utf-8")
    candidate_csv = tmp_path / "candidates.csv"
    candidate_csv.write_text("smiles\nCCN\n", encoding="utf-8")
    docking_config = target_configs / "config_docking.yml"
    docking_config.write_text("run: true\n", encoding="utf-8")
    target_master = {
        "generated_mols_path": str(target_csv),
        "target_mols_path": str(target_csv),
        "folder_to_save": str(target_run),
        "sample_size": None,
        "save_sampled_mols": True,
        "alignment": {"enabled": True, "target_coverage_percent": 90},
        "config_docking": str(docking_config),
    }
    (target_configs / "master_config_resolved.yml").write_text(
        yaml.safe_dump(target_master, sort_keys=False), encoding="utf-8"
    )

    aligned_dir = outer / "target_alignment" / "aligned_configs"
    aligned_dir.mkdir(parents=True)
    aligned_master_path = aligned_dir / "aligned_config.yml"
    aligned_master = {
        "generated_mols_path": str(candidate_csv),
        "target_mols_path": str(target_csv),
        "folder_to_save": str(outer),
        "sample_size": None,
        "save_sampled_mols": True,
        "alignment": {"enabled": False, "target_coverage_percent": 90},
        "config_docking": str(docking_config),
    }
    aligned_master_path.write_text(
        yaml.safe_dump(aligned_master, sort_keys=False), encoding="utf-8"
    )

    calls = []

    def fake_calculate(data, config, progress_callback):
        calls.append((data.copy(), dict(config)))
        if config.get("_continue_mode"):
            assert progress_callback is not None
            progress_callback(
                {"type": "stage_complete", "stage": "docking", "ok": True}
            )
        return True

    def fake_create(*args, **kwargs):
        return dict(aligned_master), aligned_master_path, aligned_dir / "thresholds.yml"

    monkeypatch.setattr(main_mod, "calculate_metrics", fake_calculate)
    monkeypatch.setattr(main_mod, "create_aligned_stage_config", fake_create)
    monkeypatch.setattr(main_mod, "finalize_global_alignment", fake_create)

    main_mod._run_pipeline_command(
        config_path=None,
        generated_mols_path=None,
        out_dir=None,
        stage=None,
        reuse_folder=False,
        force_new_folder=False,
        auto_install=False,
        show_progress=False,
        large_dataset=False,
        continue_folder=str(outer),
    )

    assert len(calls) == 2
    assert calls[0][1]["folder_to_save"] == str(target_run.resolve())
    assert calls[1][1]["folder_to_save"] == str(outer.resolve())
    assert calls[1][0]["smiles"].tolist() == ["CCN"]


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

    def test_preserves_existing_identity_and_provenance_columns(
        self, tmp_path, mock_logger
    ):
        """Schema normalization must not replace stable upstream molecule IDs."""
        input_file = tmp_path / "input.csv"
        input_file.write_text(
            "smiles,model_name,mol_idx,source_row\n"
            "CCO,test,zinc-17,17\n"
            "c1ccccc1,test,zinc-29,29\n"
        )

        result = preprocess_input_with_rdkit(
            str(input_file), tmp_path / "output", mock_logger
        )

        output_df = pd.read_csv(result)
        assert output_df["mol_idx"].tolist() == ["zinc-17", "zinc-29"]
        assert output_df["source_row"].tolist() == [17, 29]

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

    def test_empty_smiles_rows_disable_preprocessing(self, tmp_path, mock_logger):
        """Malformed CSV with empty smiles should fall back to raw input handling."""
        input_file = tmp_path / "input.csv"
        input_file.write_text(
            "smiles,model_name,extra\nCCO,test,\n,jtvae,CCN\n,hiergraphvae,CCC"
        )

        output_folder = tmp_path / "output"
        result = preprocess_input_with_rdkit(
            str(input_file), output_folder, mock_logger
        )

        assert result is None


def test_save_sampled_molecules_writes_identity_columns_only(tmp_path):
    data = pd.DataFrame(
        {
            "smiles": ["CCO"],
            "model_name": ["m1"],
            "mol_idx": ["LP-0001-00001"],
            "extra": ["should_not_be_saved"],
        }
    )

    _save_sampled_molecules(data, tmp_path, should_save=True)

    saved = pd.read_csv(tmp_path / "input" / "sampled_molecules.csv")
    assert list(saved.columns) == ["smiles", "model_name", "mol_idx"]
    assert saved.iloc[0].to_dict() == {
        "smiles": "CCO",
        "model_name": "m1",
        "mol_idx": "LP-0001-00001",
    }


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


def test_apply_cli_overrides_supports_multiple_stage_selection():
    """Repeated --stage flags should map to an ordered stage-selection override."""
    from hedgehog import main as main_mod

    config = {"generated_mols_path": "input.csv"}

    main_mod._apply_cli_overrides(
        config,
        generated_mols_path=None,
        stages=[main_mod.Stage.descriptors, main_mod.Stage.struct_filters],
    )

    assert config[main_mod.STAGE_SELECTION_KEY] == [
        "descriptors",
        "struct_filters",
    ]
    assert main_mod.STAGE_OVERRIDE_KEY not in config


def test_resolve_cli_mols_paths_filters_shell_expanded_sdf_list(tmp_path):
    """Shell-expanded --mols globs should keep molecule files and drop README."""
    from hedgehog import main as main_mod

    drugflow = tmp_path / "drugflow.sdf"
    pocket = tmp_path / "pocket2mol.sdf"
    readme = tmp_path / "README.md"
    drugflow.write_text("x")
    pocket.write_text("y")
    readme.write_text("docs")

    primary, paths = main_mod._resolve_cli_mols_paths(
        str(drugflow),
        [str(pocket), str(readme)],
    )

    assert primary == str(tmp_path)
    assert paths == [str(drugflow.resolve()), str(pocket.resolve())]


def test_setup_aizynthfinder_auto_accepts_by_default(monkeypatch, tmp_path):
    """setup aizynthfinder should auto-accept downloads without extra flags."""
    import hedgehog.setup as setup_mod
    from hedgehog import main as main_mod

    captured: dict[str, str | None] = {"auto": None}

    def _fake_ensure(_project_root):
        captured["auto"] = os.environ.get("HEDGEHOG_AUTO_INSTALL")
        return tmp_path / "config.yml"

    monkeypatch.setattr(setup_mod, "ensure_aizynthfinder", _fake_ensure)
    monkeypatch.setattr(main_mod.console, "print", lambda *args, **kwargs: None)
    monkeypatch.delenv("HEDGEHOG_AUTO_INSTALL", raising=False)

    main_mod.setup_aizynthfinder()

    assert captured["auto"] == "1"


def test_setup_aizynthfinder_no_yes_restores_prompt(monkeypatch, tmp_path):
    """--no-yes should avoid forcing auto-install confirmations."""
    import hedgehog.setup as setup_mod
    from hedgehog import main as main_mod

    captured: dict[str, str | None] = {"auto": "unexpected"}

    def _fake_ensure(_project_root):
        captured["auto"] = os.environ.get("HEDGEHOG_AUTO_INSTALL")
        return tmp_path / "config.yml"

    monkeypatch.setattr(setup_mod, "ensure_aizynthfinder", _fake_ensure)
    monkeypatch.setattr(main_mod.console, "print", lambda *args, **kwargs: None)
    monkeypatch.setenv("HEDGEHOG_AUTO_INSTALL", "1")

    main_mod.setup_aizynthfinder(yes=False)

    assert captured["auto"] is None


def test_setup_fsscore_sets_auto_install_when_yes(monkeypatch, tmp_path):
    """setup fsscore --yes should set HEDGEHOG_AUTO_INSTALL."""
    import hedgehog.setup as setup_mod
    from hedgehog import main as main_mod

    captured: dict[str, str | None] = {"auto": None}

    def _fake_ensure(_project_root):
        captured["auto"] = os.environ.get("HEDGEHOG_AUTO_INSTALL")
        return tmp_path / "modules" / "fsscore"

    monkeypatch.setattr(setup_mod, "ensure_fsscore_checkout", _fake_ensure)
    monkeypatch.setattr(main_mod.console, "print", lambda *args, **kwargs: None)
    monkeypatch.delenv("HEDGEHOG_AUTO_INSTALL", raising=False)

    main_mod.setup_fsscore(yes=True)

    assert captured["auto"] == "1"


def test_setup_nonpher_check_success(monkeypatch):
    """setup nonpher-check should report availability and return normally."""
    import hedgehog.setup as setup_mod
    from hedgehog import main as main_mod

    printed: list[str] = []

    monkeypatch.setattr(
        setup_mod,
        "check_nonpher_runtime",
        lambda **kwargs: SimpleNamespace(available=True, detail="ok"),
    )
    monkeypatch.setattr(
        main_mod.console,
        "print",
        lambda *args, **kwargs: printed.append(" ".join(str(a) for a in args)),
    )

    main_mod.setup_nonpher_check(python_bin="/tmp/nonpher/bin/python")

    assert any("Nonpher runtime is available." in line for line in printed)


def test_setup_nonpher_check_failure_exits_with_guidance(monkeypatch):
    """setup nonpher-check should exit code 1 and print Linux setup guidance."""
    import hedgehog.setup as setup_mod
    from hedgehog import main as main_mod

    printed: list[str] = []

    monkeypatch.setattr(
        setup_mod,
        "check_nonpher_runtime",
        lambda **kwargs: SimpleNamespace(
            available=False,
            detail="No module named nonpher",
        ),
    )
    monkeypatch.setattr(
        setup_mod,
        "nonpher_lobachevsky_setup_commands",
        lambda: ["ssh lobachevsky", "conda create -n hedgehog-nonpher python=3.10 -y"],
    )
    monkeypatch.setattr(
        main_mod.console,
        "print",
        lambda *args, **kwargs: printed.append(" ".join(str(a) for a in args)),
    )

    try:
        main_mod.setup_nonpher_check()
    except main_mod.typer.Exit as exc:
        assert exc.exit_code == 1
    else:
        raise AssertionError("Expected typer.Exit(code=1)")

    assert any("ssh lobachevsky" in line for line in printed)


def test_run_uses_single_progress_task_and_consistent_stage_numbers(
    tmp_path, monkeypatch
):
    """CLI progress should reuse one task and keep stage numbering monotonic."""
    from hedgehog import main as main_mod

    _FakeProgress.instances.clear()

    results_dir = tmp_path / "results"
    input_path = tmp_path / "input.csv"
    input_path.write_text("smiles,model_name,mol_idx\nCCO,m1,0\n", encoding="utf-8")

    base_config = {
        "folder_to_save": str(results_dir),
        "generated_mols_path": str(input_path),
        "save_sampled_mols": False,
    }
    prepared_df = pd.DataFrame(
        {"smiles": ["CCO"], "model_name": ["m1"], "mol_idx": [0]}
    )

    monkeypatch.setattr(main_mod, "_display_banner", lambda: None)
    monkeypatch.setattr(main_mod, "_plain_output_enabled", lambda: False)
    monkeypatch.setattr(
        main_mod, "load_config", lambda *args, **kwargs: base_config.copy()
    )
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
    monkeypatch.setattr(
        main_mod, "_save_sampled_molecules", lambda *args, **kwargs: None
    )
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
        ctx=SimpleNamespace(invoked_subcommand=None, args=[]),
        config_path="unused.yml",
        generated_mols_path=None,
        out_dir=None,
        stage=None,
        reuse_folder=False,
        force_new_folder=False,
        auto_install=False,
        show_progress=True,
        large_dataset=False,
    )

    progress_instance = _FakeProgress.instances[-1]
    assert len(progress_instance.add_calls) == 1
    assert len(progress_instance.reset_calls) == 2
    assert all(call["start"] is True for call in progress_instance.reset_calls)
    assert all(call["completed"] == 0 for call in progress_instance.reset_calls)

    descriptions = [
        call["description"]
        for call in progress_instance.update_calls
        if "description" in call
    ]
    assert any(desc.startswith("1/2 - Prep") for desc in descriptions)
    assert any(desc.startswith("2/2 - Descriptors") for desc in descriptions)
    assert not any(desc.startswith("0/2 - ") for desc in descriptions)


def test_run_progress_strips_duplicate_stage_prefix(tmp_path, monkeypatch):
    """StructFilters progress should not duplicate stage name in description."""
    from hedgehog import main as main_mod

    _FakeProgress.instances.clear()

    results_dir = tmp_path / "results"
    input_path = tmp_path / "input.csv"
    input_path.write_text("smiles,model_name,mol_idx\nCCO,m1,0\n", encoding="utf-8")

    base_config = {
        "folder_to_save": str(results_dir),
        "generated_mols_path": str(input_path),
        "save_sampled_mols": False,
    }
    prepared_df = pd.DataFrame(
        {"smiles": ["CCO"], "model_name": ["m1"], "mol_idx": [0]}
    )

    monkeypatch.setattr(main_mod, "_display_banner", lambda: None)
    monkeypatch.setattr(main_mod, "_plain_output_enabled", lambda: False)
    monkeypatch.setattr(
        main_mod, "load_config", lambda *args, **kwargs: base_config.copy()
    )
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
    monkeypatch.setattr(
        main_mod, "_save_sampled_molecules", lambda *args, **kwargs: None
    )
    monkeypatch.setattr(main_mod, "Progress", _FakeProgress)

    def _fake_calculate_metrics(data, config, progress_callback):
        progress_callback(
            {
                "type": "stage_start",
                "stage": "struct_filters",
                "stage_index": 1,
                "total_stages": 1,
                "message": "StructFilters: common_alerts",
            }
        )
        progress_callback(
            {
                "type": "stage_progress",
                "stage": "struct_filters",
                "stage_index": 1,
                "total_stages": 1,
                "current": 1,
                "total": 10,
                "message": "StructFilters: common_alerts",
            }
        )
        return True

    monkeypatch.setattr(main_mod, "calculate_metrics", _fake_calculate_metrics)

    main_mod.run(
        ctx=SimpleNamespace(invoked_subcommand=None, args=[]),
        config_path="unused.yml",
        generated_mols_path=None,
        out_dir=None,
        stage=None,
        reuse_folder=False,
        force_new_folder=False,
        auto_install=False,
        show_progress=True,
        large_dataset=False,
    )

    progress_instance = _FakeProgress.instances[-1]
    descriptions = [
        call["description"]
        for call in progress_instance.update_calls
        if "description" in call
    ]
    assert any(
        desc.startswith("1/1 - StructFilters · common_alerts") for desc in descriptions
    )
    assert not any("StructFilters · StructFilters:" in desc for desc in descriptions)


def test_run_disables_progress_bar_by_default(tmp_path, monkeypatch):
    """CLI should not create Rich progress/task unless --progress is enabled."""
    from hedgehog import main as main_mod

    results_dir = tmp_path / "results"
    input_path = tmp_path / "input.csv"
    input_path.write_text("smiles,model_name,mol_idx\nCCO,m1,0\n", encoding="utf-8")

    base_config = {
        "folder_to_save": str(results_dir),
        "generated_mols_path": str(input_path),
        "save_sampled_mols": False,
    }
    prepared_df = pd.DataFrame(
        {"smiles": ["CCO"], "model_name": ["m1"], "mol_idx": [0]}
    )

    monkeypatch.setattr(main_mod, "_display_banner", lambda: None)
    monkeypatch.setattr(main_mod, "_plain_output_enabled", lambda: False)
    monkeypatch.setattr(
        main_mod, "load_config", lambda *args, **kwargs: base_config.copy()
    )
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
    monkeypatch.setattr(
        main_mod, "_save_sampled_molecules", lambda *args, **kwargs: None
    )

    def _progress_should_not_be_created(*args, **kwargs):
        raise AssertionError(
            "Progress should not be initialized when --progress is not set"
        )

    monkeypatch.setattr(main_mod, "Progress", _progress_should_not_be_created)

    captured: dict[str, object] = {}

    def _fake_calculate_metrics(data, config, progress_callback):
        captured["progress_callback"] = progress_callback
        return True

    monkeypatch.setattr(main_mod, "calculate_metrics", _fake_calculate_metrics)

    main_mod.run(
        ctx=SimpleNamespace(invoked_subcommand=None, args=[]),
        config_path="unused.yml",
        generated_mols_path=None,
        out_dir=None,
        stage=None,
        reuse_folder=False,
        force_new_folder=False,
        auto_install=False,
        show_progress=False,
        large_dataset=False,
    )

    assert captured["progress_callback"] is None


def test_large_dataset_defaults_to_compute_only_statistics_path(tmp_path, monkeypatch):
    """Large-dataset mode should default to the no-drop statistics stage set."""
    import hedgehog.main as main_mod

    input_path = tmp_path / "input.csv"
    input_path.write_text("smiles,model_name,mol_idx\nCCO,m1,0\n", encoding="utf-8")
    base_config = {
        "folder_to_save": str(tmp_path / "configured"),
        "generated_mols_path": str(input_path),
    }

    monkeypatch.setattr(main_mod, "_display_banner", lambda: None)
    monkeypatch.setattr(
        main_mod, "load_config", lambda *args, **kwargs: base_config.copy()
    )
    monkeypatch.setattr(main_mod, "_resolve_config_paths", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        main_mod.LoggerSingleton,
        "configure_log_directory",
        lambda self, folder_to_save: None,
    )

    captured: dict[str, object] = {}

    def _fake_calculate_metrics(data, config, progress_callback):
        captured["data"] = data
        captured["config"] = config
        captured["progress_callback"] = progress_callback
        return True

    monkeypatch.setattr(main_mod, "calculate_metrics", _fake_calculate_metrics)

    main_mod._run_pipeline_command(
        config_path="unused.yml",
        generated_mols_path=None,
        out_dir=str(tmp_path / "run"),
        stage=None,
        reuse_folder=False,
        force_new_folder=False,
        auto_install=False,
        show_progress=False,
        large_dataset=True,
    )

    config = captured["config"]
    assert captured["data"] is None
    assert config[main_mod.STAGE_SELECTION_KEY] == [
        "mol_prep",
        "descriptors",
        "struct_filters",
        "synthesis",
    ]
    assert config["large_dataset_mode"] is True
    assert config["large_dataset_filter_data"] is False
    assert config["large_dataset_enable_all_filters"] is True
    assert config["folder_to_save"] == str(tmp_path / "run_1")


def test_large_dataset_still_rejects_docking_stage(tmp_path, monkeypatch):
    """Large-dataset mode allows synthesis but not docking-heavy stages."""
    import hedgehog.main as main_mod

    monkeypatch.setattr(main_mod, "_display_banner", lambda: None)
    monkeypatch.setattr(
        main_mod,
        "load_config",
        lambda *args, **kwargs: {
            "folder_to_save": str(tmp_path / "run"),
            "generated_mols_path": str(tmp_path / "input.csv"),
        },
    )
    monkeypatch.setattr(main_mod, "_resolve_config_paths", lambda *args, **kwargs: None)

    with pytest.raises(main_mod.typer.Exit):
        main_mod._run_pipeline_command(
            config_path="unused.yml",
            generated_mols_path=None,
            out_dir=str(tmp_path / "run"),
            stage=[main_mod.Stage.docking],
            reuse_folder=False,
            force_new_folder=False,
            auto_install=False,
            show_progress=False,
            large_dataset=True,
        )


def test_run_subcommand_delegates_to_pipeline_command(monkeypatch):
    """Explicit run subcommand should delegate to pipeline command helper."""
    from hedgehog import main as main_mod

    captured: dict[str, object] = {}

    def _fake_run_pipeline_command(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(main_mod, "_run_pipeline_command", _fake_run_pipeline_command)

    main_mod.run_command(
        ctx=SimpleNamespace(args=[]),
        config_path="cfg.yml",
        generated_mols_path="input.csv",
        out_dir="results/x",
        stage=[main_mod.Stage.docking],
        reuse_folder=True,
        force_new_folder=False,
        auto_install=True,
        show_progress=True,
        large_dataset=False,
    )

    assert captured == {
        "config_path": "cfg.yml",
        "generated_mols_path": "input.csv",
        "generated_mols_paths": None,
        "out_dir": "results/x",
        "stage": [main_mod.Stage.docking],
        "reuse_folder": True,
        "force_new_folder": False,
        "auto_install": True,
        "show_progress": True,
        "large_dataset": False,
    }
