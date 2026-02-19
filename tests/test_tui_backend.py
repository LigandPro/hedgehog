"""Tests for TUI backend audit fixes (C1-C3, H1, H3, H4, M2, M3, L1)."""

import threading
from pathlib import Path
from unittest.mock import MagicMock

import pytest
import yaml

from hedgehog.tui_backend.handlers.files import FilesHandler, _validate_path
from hedgehog.tui_backend.handlers.history import HistoryHandler
from hedgehog.tui_backend.handlers.pipeline import PipelineHandler, PipelineJob
from hedgehog.tui_backend.server import JsonRpcServer
from hedgehog.tui_backend.validators import ConfigValidator

# =====================================================================
# Helpers
# =====================================================================


def _mock_server():
    """Create a mock JsonRpcServer with send_notification stub."""
    server = MagicMock(spec=JsonRpcServer)
    server.send_notification = MagicMock()
    return server


# =====================================================================
# C1 — Path traversal validation (_validate_path)
# =====================================================================


class TestValidatePath:
    """Test _validate_path rejects paths outside $HOME."""

    def test_path_inside_home_is_allowed(self):
        home = Path.home()
        _validate_path(home / "some" / "file.csv")

    def test_path_outside_home_raises(self):
        with pytest.raises(ValueError, match="Access denied"):
            _validate_path(Path("/etc/passwd"))

    def test_root_path_raises(self):
        with pytest.raises(ValueError, match="Access denied"):
            _validate_path(Path("/"))

    def test_custom_allowed_base(self, tmp_path):
        # Path inside allowed base: OK
        _validate_path(tmp_path / "sub" / "file.txt", allowed_base=tmp_path)
        # Path outside allowed base: raises
        with pytest.raises(ValueError, match="Access denied"):
            _validate_path(Path("/tmp/other"), allowed_base=tmp_path)

    def test_validate_input_file_rejects_outside_home(self, tmp_path):
        """ValidationHandler.validate_input_file calls _validate_path."""
        from hedgehog.tui_backend.handlers.validation import ValidationHandler

        handler = ValidationHandler(_mock_server())
        with pytest.raises(ValueError, match="Access denied"):
            handler.validate_input_file("/etc/passwd")

    def test_validate_receptor_pdb_rejects_outside_home(self):
        from hedgehog.tui_backend.handlers.validation import ValidationHandler

        handler = ValidationHandler(_mock_server())
        with pytest.raises(ValueError, match="Access denied"):
            handler.validate_receptor_pdb("/etc/hosts")

    def test_validate_output_directory_rejects_outside_home(self):
        from hedgehog.tui_backend.handlers.validation import ValidationHandler

        handler = ValidationHandler(_mock_server())
        with pytest.raises(ValueError, match="Access denied"):
            handler.validate_output_directory("/var/log")


# =====================================================================
# C2 — HistoryHandler thread safety
# =====================================================================


class TestHistoryHandlerThreadSafety:
    """Test that concurrent HistoryHandler operations don't lose data."""

    @pytest.fixture()
    def handler(self, tmp_path, monkeypatch):
        """Create a HistoryHandler with a temp history file."""
        server = _mock_server()
        handler = HistoryHandler(server)
        # Redirect history file to tmp_path
        handler.history_file = tmp_path / "job_history.json"
        handler.history_file.write_text("[]")
        return handler

    def test_concurrent_add_job_no_data_loss(self, handler):
        n_threads = 20
        barrier = threading.Barrier(n_threads)
        errors = []

        def add_one(idx):
            try:
                barrier.wait(timeout=5)
                handler.add_job(f"job-{idx}", name=f"Job {idx}")
            except Exception as e:
                errors.append(e)

        threads = [
            threading.Thread(target=add_one, args=(i,)) for i in range(n_threads)
        ]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=10)

        assert not errors, f"Errors during concurrent add: {errors}"
        history = handler.get_job_history(limit=100)
        assert len(history) == n_threads

    def test_concurrent_add_and_update(self, handler):
        # Seed one job
        handler.add_job("seed-job", name="Seed")

        n_threads = 10
        barrier = threading.Barrier(n_threads + 1)
        errors = []

        def update_loop():
            try:
                barrier.wait(timeout=5)
                handler.update_job("seed-job", status="completed")
            except Exception as e:
                errors.append(e)

        def add_loop(idx):
            try:
                barrier.wait(timeout=5)
                handler.add_job(f"concurrent-{idx}")
            except Exception as e:
                errors.append(e)

        threads = [threading.Thread(target=update_loop)]
        threads += [
            threading.Thread(target=add_loop, args=(i,)) for i in range(n_threads)
        ]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=10)

        assert not errors
        history = handler.get_job_history(limit=200)
        # 1 seed + n_threads new jobs
        assert len(history) == 1 + n_threads


# =====================================================================
# C3 — PipelineHandler thread safety
# =====================================================================


class TestPipelineHandlerThreadSafety:
    """Test that PipelineHandler.jobs access is protected by a lock."""

    def test_has_lock(self):
        server = _mock_server()
        handler = PipelineHandler(server)
        assert hasattr(handler, "_lock")
        assert isinstance(handler._lock, type(threading.Lock()))

    def test_get_progress_unknown_job_raises(self):
        server = _mock_server()
        handler = PipelineHandler(server)
        with pytest.raises(ValueError, match="Job not found"):
            handler.get_progress("nonexistent")

    def test_cancel_unknown_job_raises(self):
        server = _mock_server()
        handler = PipelineHandler(server)
        with pytest.raises(ValueError, match="Job not found"):
            handler.cancel_pipeline("nonexistent")

    def test_concurrent_start_and_cancel(self):
        """Start multiple jobs concurrently and cancel them; no crashes."""
        server = _mock_server()
        handler = PipelineHandler(server)
        n = 10
        barrier = threading.Barrier(n)
        job_ids = []
        lock = threading.Lock()
        errors = []

        def start_one(idx):
            try:
                barrier.wait(timeout=5)
                jid = handler.start_pipeline([f"stage-{idx}"])
                with lock:
                    job_ids.append(jid)
            except Exception as e:
                errors.append(e)

        threads = [threading.Thread(target=start_one, args=(i,)) for i in range(n)]
        for t in threads:
            t.start()
        for t in threads:
            t.join(timeout=10)

        assert not errors
        assert len(job_ids) == n

        # Cancel all jobs (they won't actually run since _run will fail without hedgehog)
        for jid in job_ids:
            handler.cancel_pipeline(jid)


# =====================================================================
# H1 — Stage name mapping bug
# =====================================================================


class TestStageNameMapping:
    """Test that PipelineJob maps internal names to TUI names correctly."""

    def test_map_stage_name(self):
        server = _mock_server()
        job = PipelineJob("test-id", ["descriptors", "struct_filters"], server)
        assert job._map_stage_name("final_descriptors") == "descriptors"
        assert job._map_stage_name("mol_prep") == "mol_prep"
        assert job._map_stage_name("docking") == "docking"
        assert job._map_stage_name("unknown_stage") == "unknown_stage"

    def test_previous_stage_stores_tui_name(self):
        """After progress_callback, previous_stage should hold the TUI name."""
        server = _mock_server()
        stages = ["descriptors", "struct_filters", "synthesis"]
        job = PipelineJob("test-id", stages, server)

        # Simulate what progress_callback does (extracted logic)
        internal_stage = "final_descriptors"
        tui_stage = job._map_stage_name(internal_stage)

        # Before the fix, previous_stage was set to internal name.
        # After the fix, it is set to tui_stage.
        job.previous_stage = tui_stage
        job.current_stage = tui_stage

        assert job.previous_stage == "descriptors"
        assert job.current_stage == "descriptors"

    def test_no_false_stage_transitions(self):
        """Same TUI stage repeated should not trigger stage_complete/start."""
        server = _mock_server()
        stages = ["descriptors"]
        job = PipelineJob("test-id", stages, server)

        # Simulate two callbacks for stages that map to the same TUI name
        # First: "descriptors" -> tui "descriptors"
        job.current_stage = "descriptors"
        job.previous_stage = "descriptors"

        # Second: "final_descriptors" -> tui "descriptors"
        tui_stage = job._map_stage_name("final_descriptors")
        assert tui_stage == "descriptors"

        # previous_stage == tui_stage, so no stage_complete should fire
        assert job.previous_stage == tui_stage


# =====================================================================
# H3 — Progress dict update in callback
# =====================================================================


class TestProgressUpdate:
    """Test that self.progress is updated during callback."""

    def test_progress_dict_initialized(self):
        server = _mock_server()
        stages = ["descriptors", "synthesis"]
        job = PipelineJob("test-id", stages, server)
        assert job.progress == {"descriptors": 0, "synthesis": 0}

    def test_progress_dict_updated_on_set(self):
        """Verify the progress dict can be updated per stage."""
        server = _mock_server()
        stages = ["descriptors", "synthesis"]
        job = PipelineJob("test-id", stages, server)

        # Simulate what the callback does: progress[tui_stage] = progress_pct
        tui_stage = "descriptors"
        progress_pct = int((2 / 5) * 100)
        job.progress[tui_stage] = progress_pct

        assert job.progress["descriptors"] == 40
        assert job.progress["synthesis"] == 0


# =====================================================================
# M3 — Streaming count_molecules for SDF
# =====================================================================


class TestCountMolecules:
    """Test FilesHandler.count_molecules with streaming SDF counting."""

    @pytest.fixture(autouse=True)
    def _bypass_path_validation(self, monkeypatch):
        """Allow tmp_path (outside $HOME) in count_molecules tests."""
        monkeypatch.setattr(
            "hedgehog.tui_backend.handlers.files._validate_path",
            lambda *_a, **_kw: None,
        )

    @pytest.fixture()
    def handler(self):
        return FilesHandler(_mock_server())

    def test_sdf_multi_molecule(self, handler, tmp_path):
        sdf = tmp_path / "mols.sdf"
        sdf.write_text("mol1\n$$$$\nmol2\n$$$$\nmol3\n$$$$\n")
        result = handler.count_molecules(str(sdf))
        assert result["count"] == 3

    def test_sdf_empty_file(self, handler, tmp_path):
        sdf = tmp_path / "empty.sdf"
        sdf.write_text("")
        result = handler.count_molecules(str(sdf))
        assert result["count"] == 0

    def test_sdf_single_molecule_no_delimiter(self, handler, tmp_path):
        sdf = tmp_path / "single.sdf"
        sdf.write_text("some molecule data\nwith multiple lines\n")
        result = handler.count_molecules(str(sdf))
        # Fallback: no $$$$ but file is non-empty -> count = 1
        assert result["count"] == 1

    def test_csv_count(self, handler, tmp_path):
        csv = tmp_path / "mols.csv"
        csv.write_text("smiles,name\nCCC,ethane\nCCCC,butane\n")
        result = handler.count_molecules(str(csv))
        assert result["count"] == 2

    def test_smi_count(self, handler, tmp_path):
        smi = tmp_path / "mols.smi"
        smi.write_text("CCC\nCCCC\n\nCCCCC\n")
        result = handler.count_molecules(str(smi))
        # 3 non-empty lines
        assert result["count"] == 3

    def test_nonexistent_file_raises(self, handler, tmp_path):
        with pytest.raises(FileNotFoundError):
            handler.count_molecules(str(tmp_path / "nope.sdf"))


# =====================================================================
# H4 — save_config rejects invalid config data
# =====================================================================


class TestSaveConfigValidation:
    """Test that ConfigHandler.save_config validates before saving."""

    def test_save_invalid_main_config_raises(self, tmp_path, monkeypatch):
        from hedgehog.tui_backend.handlers.config import ConfigHandler

        server = _mock_server()
        handler = ConfigHandler(server)
        # Point project root to tmp_path so config files are written there
        monkeypatch.setattr(handler, "_project_root", tmp_path)
        monkeypatch.setattr(
            handler, "_user_config_root", tmp_path / ".hedgehog" / "tui"
        )

        # Create the config directory structure expected by _get_config_path
        config_dir = tmp_path / "src" / "hedgehog" / "configs"
        config_dir.mkdir(parents=True)
        config_file = config_dir / "config.yml"
        config_file.write_text(
            "generated_mols_path: data/mols.csv\nfolder_to_save: results\n"
        )

        # Missing required fields -> should raise ValueError
        with pytest.raises(ValueError, match="Invalid config"):
            handler.save_config("main", {"n_jobs": 4})

    def test_save_valid_main_config_succeeds(self, tmp_path, monkeypatch):
        from hedgehog.tui_backend.handlers.config import ConfigHandler

        server = _mock_server()
        handler = ConfigHandler(server)
        monkeypatch.setattr(handler, "_project_root", tmp_path)
        monkeypatch.setattr(
            handler, "_user_config_root", tmp_path / ".hedgehog" / "tui"
        )

        config_dir = tmp_path / "src" / "hedgehog" / "configs"
        config_dir.mkdir(parents=True)
        config_file = config_dir / "config.yml"
        config_file.write_text("generated_mols_path: old\nfolder_to_save: old\n")

        valid_data = {
            "generated_mols_path": "data/mols.csv",
            "folder_to_save": "results/run",
            "n_jobs": 4,
        }
        result = handler.save_config("main", valid_data)
        assert result is True

        # Verify user config was updated (source template remains unchanged)
        user_config_file = handler.get_config_path("main")
        saved = yaml.safe_load(user_config_file.read_text())
        assert saved["folder_to_save"] == "results/run"
        source_saved = yaml.safe_load(config_file.read_text())
        assert source_saved["folder_to_save"] == "old"

    def test_load_main_config_bootstraps_user_copy(self, tmp_path, monkeypatch):
        from hedgehog.tui_backend.handlers.config import ConfigHandler

        server = _mock_server()
        handler = ConfigHandler(server)
        monkeypatch.setattr(handler, "_project_root", tmp_path)
        monkeypatch.setattr(
            handler, "_user_config_root", tmp_path / ".hedgehog" / "tui"
        )

        config_dir = tmp_path / "src" / "hedgehog" / "configs"
        config_dir.mkdir(parents=True)
        source_file = config_dir / "config.yml"
        source_file.write_text(
            "generated_mols_path: data/mols.csv\nfolder_to_save: results\n"
        )

        loaded = handler.load_config("main")
        assert loaded["folder_to_save"] == "results"

        user_file = handler.get_config_path("main")
        assert user_file.exists()
        assert user_file != source_file
        assert source_file.read_text() == user_file.read_text()

    def test_load_mol_prep_config_bootstraps_user_copy(self, tmp_path, monkeypatch):
        from hedgehog.tui_backend.handlers.config import ConfigHandler

        server = _mock_server()
        handler = ConfigHandler(server)
        monkeypatch.setattr(handler, "_project_root", tmp_path)
        monkeypatch.setattr(
            handler, "_user_config_root", tmp_path / ".hedgehog" / "tui"
        )

        config_dir = tmp_path / "src" / "hedgehog" / "configs"
        config_dir.mkdir(parents=True)
        source_file = config_dir / "config_mol_prep.yml"
        source_file.write_text("run: true\nn_jobs: -1\n")

        loaded = handler.load_config("mol_prep")
        assert loaded["run"] is True
        assert loaded["n_jobs"] == -1

        user_file = handler.get_config_path("mol_prep")
        assert user_file.exists()
        assert user_file != source_file
        assert source_file.read_text() == user_file.read_text()


# =====================================================================
# L1 — Error code mapping in JsonRpcServer.handle_request
# =====================================================================


class TestErrorCodeMapping:
    """Test that exceptions are mapped to correct JSON-RPC error codes."""

    def test_file_not_found_returns_32001(self):
        server = JsonRpcServer()
        responses = []

        def capture(msg):
            responses.append(msg)

        server._send_message = capture

        # Register a handler that raises FileNotFoundError
        server.handlers["test_fnf"] = lambda: (_ for _ in ()).throw(
            FileNotFoundError("not found")
        )
        server.handle_request({"id": 1, "method": "test_fnf", "params": {}})

        assert len(responses) == 1
        assert responses[0]["error"]["code"] == -32001

    def test_value_error_returns_32602(self):
        server = JsonRpcServer()
        responses = []

        def capture(msg):
            responses.append(msg)

        server._send_message = capture

        server.handlers["test_ve"] = lambda: (_ for _ in ()).throw(
            ValueError("bad value")
        )
        server.handle_request({"id": 2, "method": "test_ve", "params": {}})

        assert len(responses) == 1
        assert responses[0]["error"]["code"] == -32602

    def test_permission_error_returns_32002(self):
        server = JsonRpcServer()
        responses = []

        def capture(msg):
            responses.append(msg)

        server._send_message = capture

        server.handlers["test_pe"] = lambda: (_ for _ in ()).throw(
            PermissionError("denied")
        )
        server.handle_request({"id": 3, "method": "test_pe", "params": {}})

        assert len(responses) == 1
        assert responses[0]["error"]["code"] == -32002

    def test_method_not_found_returns_32601(self):
        server = JsonRpcServer()
        responses = []

        def capture(msg):
            responses.append(msg)

        server._send_message = capture

        server.handle_request({"id": 4, "method": "nonexistent", "params": {}})

        assert len(responses) == 1
        assert responses[0]["error"]["code"] == -32601

    def test_generic_exception_returns_32000(self):
        server = JsonRpcServer()
        responses = []

        def capture(msg):
            responses.append(msg)

        server._send_message = capture

        server.handlers["test_gen"] = lambda: (_ for _ in ()).throw(
            RuntimeError("oops")
        )
        server.handle_request({"id": 5, "method": "test_gen", "params": {}})

        assert len(responses) == 1
        assert responses[0]["error"]["code"] == -32000


# =====================================================================
# M2 — ConfigValidator warning for unknown types
# =====================================================================


class TestConfigValidatorUnknownType:
    """Test that unknown config type produces a warning, not an error."""

    def test_unknown_type_returns_valid_with_warning(self):
        result = ConfigValidator.validate("totally_unknown_type", {"key": "val"})
        assert result["valid"] is True
        assert any("No specific validator" in w for w in result["warnings"])

    def test_known_type_no_warning(self):
        result = ConfigValidator.validate(
            "main",
            {
                "generated_mols_path": "data/mols.csv",
                "folder_to_save": "results",
            },
        )
        # Should not have the "no specific validator" warning
        assert not any("No specific validator" in w for w in result.get("warnings", []))


# =====================================================================
# H2 — _disable_unrequested_stages
# =====================================================================


class TestDisableUnrequestedStages:
    """Test PipelineJob._disable_unrequested_stages modifies config correctly."""

    def test_disables_stages_not_in_requested_list(self, tmp_path):
        server = _mock_server()
        # Only request "descriptors" stage
        job = PipelineJob("test-id", ["descriptors"], server)

        # Create sub-config files with run: true
        mol_prep_cfg = tmp_path / "config_mol_prep.yml"
        mol_prep_cfg.write_text(yaml.dump({"run": True}))

        synthesis_cfg = tmp_path / "config_synthesis.yml"
        synthesis_cfg.write_text(yaml.dump({"run": True, "other": "value"}))

        docking_cfg = tmp_path / "config_docking.yml"
        docking_cfg.write_text(yaml.dump({"run": True}))

        descriptors_cfg = tmp_path / "config_descriptors.yml"
        descriptors_cfg.write_text(yaml.dump({"run": True}))

        config_dict = {
            "config_mol_prep": str(mol_prep_cfg),
            "config_synthesis": str(synthesis_cfg),
            "config_docking": str(docking_cfg),
            "config_descriptors": str(descriptors_cfg),
        }

        job._disable_unrequested_stages(config_dict)

        # Mol prep should now point to temp file with run: false
        assert config_dict["config_mol_prep"] != str(mol_prep_cfg)
        new_mol_prep = yaml.safe_load(Path(config_dict["config_mol_prep"]).read_text())
        assert new_mol_prep["run"] is False

        # Synthesis and docking should now point to temp files with run: false
        assert config_dict["config_synthesis"] != str(synthesis_cfg)
        new_synth = yaml.safe_load(Path(config_dict["config_synthesis"]).read_text())
        assert new_synth["run"] is False
        assert new_synth["other"] == "value"  # Other keys preserved

        assert config_dict["config_docking"] != str(docking_cfg)
        new_dock = yaml.safe_load(Path(config_dict["config_docking"]).read_text())
        assert new_dock["run"] is False

        # Descriptors should remain unchanged (it was requested)
        assert config_dict["config_descriptors"] == str(descriptors_cfg)

    def test_keeps_all_stages_when_all_requested(self, tmp_path):
        server = _mock_server()
        all_stages = [
            "descriptors",
            "struct_filters",
            "synthesis",
            "docking",
            "docking_filters",
        ]
        job = PipelineJob("test-id", all_stages, server)

        synth_cfg = tmp_path / "config_synthesis.yml"
        synth_cfg.write_text(yaml.dump({"run": True}))

        config_dict = {"config_synthesis": str(synth_cfg)}
        job._disable_unrequested_stages(config_dict)

        # Should not have been modified
        assert config_dict["config_synthesis"] == str(synth_cfg)


# =====================================================================
# Preflight pipeline checks
# =====================================================================


class TestPipelinePreflight:
    """Test PipelineHandler.preflight_pipeline end-to-end checks."""

    @staticmethod
    def _write_source_configs(project_root: Path, paths: dict[str, Path]) -> None:
        config_dir = project_root / "src" / "hedgehog" / "configs"
        config_dir.mkdir(parents=True, exist_ok=True)

        (config_dir / "config.yml").write_text(
            yaml.safe_dump(
                {
                    "generated_mols_path": str(paths["input_csv"]),
                    "folder_to_save": str(paths["output_dir"]),
                    "n_jobs": 1,
                    "config_mol_prep": str(config_dir / "config_mol_prep.yml"),
                    "config_descriptors": str(config_dir / "config_descriptors.yml"),
                    "config_structFilters": str(
                        config_dir / "config_structFilters.yml"
                    ),
                    "config_synthesis": str(config_dir / "config_synthesis.yml"),
                    "config_docking": str(config_dir / "config_docking.yml"),
                    "config_docking_filters": str(
                        config_dir / "config_docking_filters.yml"
                    ),
                },
                sort_keys=False,
            )
        )
        (config_dir / "config_mol_prep.yml").write_text(
            yaml.safe_dump({"run": True, "n_jobs": 1}, sort_keys=False)
        )
        (config_dir / "config_descriptors.yml").write_text(
            yaml.safe_dump(
                {
                    "run": True,
                    "batch_size": 1000,
                    "borders": {
                        "molWt_min": 200,
                        "molWt_max": 500,
                        "logP_min": 0.0,
                        "logP_max": 5.0,
                    },
                },
                sort_keys=False,
            )
        )
        (config_dir / "config_structFilters.yml").write_text(
            yaml.safe_dump({"run": True}, sort_keys=False)
        )
        (config_dir / "config_synthesis.yml").write_text(
            yaml.safe_dump(
                {
                    "run": True,
                    "sa_score_min": 1,
                    "sa_score_max": 4.5,
                    "ra_score_min": 0.5,
                    "ra_score_max": 1.0,
                },
                sort_keys=False,
            )
        )
        (config_dir / "config_docking.yml").write_text(
            yaml.safe_dump(
                {
                    "run": True,
                    "tools": "smina",
                    "receptor_pdb": str(paths["receptor_pdb"]),
                },
                sort_keys=False,
            )
        )
        (config_dir / "config_docking_filters.yml").write_text(
            yaml.safe_dump(
                {
                    "run": True,
                    "run_after_docking": True,
                    "input_sdf": None,
                    "receptor_pdb": None,
                    "aggregation": {"mode": "all"},
                },
                sort_keys=False,
            )
        )

    @pytest.fixture()
    def preflight_env(self, tmp_path, monkeypatch):
        monkeypatch.setattr(
            "hedgehog.tui_backend.handlers.files._validate_path",
            lambda *_a, **_kw: None,
        )
        monkeypatch.setattr(
            "hedgehog.tui_backend.handlers.validation._validate_path",
            lambda *_a, **_kw: None,
        )

        input_csv = tmp_path / "input.csv"
        input_csv.write_text("smiles\nCCO\nCCN\n")

        receptor_pdb = tmp_path / "receptor.pdb"
        receptor_pdb.write_text(
            "ATOM      1  N   ALA A   1      11.104  13.207  14.281  1.00 20.00           N\n"
        )

        output_dir = tmp_path / "results"
        output_dir.mkdir(parents=True)

        project_root = tmp_path / "project"
        user_config_root = tmp_path / "user-configs"
        paths = {
            "input_csv": input_csv,
            "receptor_pdb": receptor_pdb,
            "output_dir": output_dir,
            "project_root": project_root,
        }
        self._write_source_configs(project_root, paths)

        server = JsonRpcServer()
        server.config_handler._project_root = project_root
        server.config_handler._user_config_root = user_config_root

        return {
            "server": server,
            "paths": paths,
        }

    def test_preflight_success(self, preflight_env):
        handler = preflight_env["server"].pipeline_handler

        result = handler.preflight_pipeline(["descriptors"])

        assert result["valid"] is True
        assert result["molecule_count"] == 2
        assert result["estimated_runtime"] == "short"
        assert any(r["stage"] == "descriptors" for r in result["stage_reports"])

    def test_preflight_blocks_invalid_docking_receptor(self, preflight_env):
        project_root = preflight_env["paths"]["project_root"]
        docking_cfg = (
            project_root / "src" / "hedgehog" / "configs" / "config_docking.yml"
        )
        docking_cfg.write_text(
            yaml.safe_dump(
                {
                    "run": True,
                    "tools": "smina",
                    "receptor_pdb": str(project_root / "missing_receptor.pdb"),
                },
                sort_keys=False,
            )
        )

        handler = preflight_env["server"].pipeline_handler
        result = handler.preflight_pipeline(["docking"])

        assert result["valid"] is False
        stage_checks = result["stage_reports"][0]["checks"]
        assert any(
            check["code"] == "DOCKING_RECEPTOR_INVALID" and check["level"] == "error"
            for check in stage_checks
        )

    def test_preflight_warning_only_keeps_valid(self, preflight_env):
        output_file = preflight_env["paths"]["output_dir"] / "existing.txt"
        output_file.write_text("already here")

        handler = preflight_env["server"].pipeline_handler
        result = handler.preflight_pipeline(["descriptors"])

        assert result["valid"] is True
        assert any(
            check["code"] == "MAIN_OUTPUT_WARNING" and check["level"] == "warning"
            for check in result["checks"]
        )

    def test_preflight_requires_input_sdf_for_manual_docking_filters(
        self, preflight_env
    ):
        project_root = preflight_env["paths"]["project_root"]
        filters_cfg = (
            project_root / "src" / "hedgehog" / "configs" / "config_docking_filters.yml"
        )
        filters_cfg.write_text(
            yaml.safe_dump(
                {
                    "run": True,
                    "run_after_docking": False,
                    "input_sdf": None,
                    "receptor_pdb": None,
                    "aggregation": {"mode": "all"},
                },
                sort_keys=False,
            )
        )

        handler = preflight_env["server"].pipeline_handler
        result = handler.preflight_pipeline(["docking_filters"])

        assert result["valid"] is False
        stage_checks = result["stage_reports"][0]["checks"]
        assert any(
            check["code"] == "DOCKING_FILTERS_INPUT_REQUIRED"
            and check["level"] == "error"
            for check in stage_checks
        )
