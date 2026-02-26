"""Pipeline execution handler for TUI backend."""

import copy
import os
import tempfile
import threading
import uuid
from pathlib import Path
from typing import TYPE_CHECKING, Any

import yaml

from ..validators import ConfigValidator
from .config import CONFIG_WORKSPACE_FILES, RUNTIME_CONFIG_KEY_MAP, ConfigHandler
from .files import FilesHandler
from .validation import ValidationHandler

if TYPE_CHECKING:
    from ..server import JsonRpcServer


class PipelineJob:
    """Represents a running pipeline job."""

    # Map internal stage names to TUI stage names
    STAGE_NAME_MAP = {
        "mol_prep": "mol_prep",
        "descriptors": "descriptors",
        "struct_filters": "struct_filters",
        "synthesis": "synthesis",
        "docking": "docking",
        "docking_filters": "docking_filters",
        "final_descriptors": "descriptors",
    }

    # Map TUI stage names to the config keys whose "run" flag they control
    TUI_STAGE_CONFIG_KEYS = {
        "mol_prep": ["config_mol_prep"],
        "descriptors": ["config_descriptors"],
        "struct_filters": ["config_structFilters"],
        "synthesis": ["config_synthesis"],
        "docking": ["config_docking"],
        "docking_filters": ["config_docking_filters"],
    }

    CONFIG_TYPE_TO_MAIN_KEY = {v: k for k, v in RUNTIME_CONFIG_KEY_MAP.items()}

    def __init__(
        self,
        job_id: str,
        stages: list[str],
        server: "JsonRpcServer",
        config_overrides: dict[str, Any] | None = None,
    ):
        self.job_id = job_id
        self.stages = stages
        self.server = server
        self.cancelled = False
        self.current_stage: str | None = None
        self.previous_stage: str | None = None
        self.progress: dict[str, int] = {s: 0 for s in stages}
        self._thread: threading.Thread | None = None
        self._temp_files: list[Path] = []
        self.config_overrides = self._validate_override_payload(config_overrides)

    @staticmethod
    def _validate_override_payload(
        config_overrides: dict[str, Any] | None,
    ) -> dict[str, dict[str, Any]]:
        if config_overrides is None:
            return {}
        if not isinstance(config_overrides, dict):
            raise ValueError("config_overrides must be a dictionary")

        allowed = set(CONFIG_WORKSPACE_FILES.keys())
        normalized: dict[str, dict[str, Any]] = {}
        for config_type, payload in config_overrides.items():
            if config_type not in allowed:
                raise ValueError(f"Unknown override config type: {config_type}")
            if not isinstance(payload, dict):
                raise ValueError(
                    f"Override payload for '{config_type}' must be a dictionary"
                )
            normalized[config_type] = payload
        return normalized

    def start(self):
        """Start the pipeline in a background thread."""
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def cancel(self):
        """Cancel the pipeline."""
        # TODO: Stages without progress events cannot be cancelled promptly because
        # cancellation is only checked in _progress_callback. Consider adding a
        # periodic cancellation check inside long-running stage loops.
        self.cancelled = True

    def _map_stage_name(self, internal_name: str) -> str:
        """Map internal stage name to TUI stage name."""
        return self.STAGE_NAME_MAP.get(internal_name, internal_name)

    def _notify(self, method: str, params: dict) -> None:
        """Send a notification via the server."""
        self.server.send_notification(method, params)

    def _log(self, level: str, message: str) -> None:
        """Send a log notification."""
        self._notify("log", {"level": level, "message": message})

    def _disable_unrequested_stages(self, config_dict: dict) -> None:
        """Disable pipeline stages not requested by the TUI user.

        For each sub-config whose TUI stage is not in self.stages, load the
        YAML, set ``run: false``, write to a temp file, and update config_dict
        to point to it so that the pipeline picks up the override.
        """
        # Collect all config keys that should remain enabled
        enabled_keys: set[str] = set()
        for tui_stage in self.stages:
            for key in self.TUI_STAGE_CONFIG_KEYS.get(tui_stage, []):
                enabled_keys.add(key)

        for tui_stage, cfg_keys in self.TUI_STAGE_CONFIG_KEYS.items():
            if tui_stage in self.stages:
                continue
            for cfg_key in cfg_keys:
                if cfg_key in enabled_keys:
                    # Another requested TUI stage shares this config key
                    continue
                cfg_path = config_dict.get(cfg_key)
                if not cfg_path:
                    continue
                try:
                    with open(cfg_path) as f:
                        sub_cfg = yaml.safe_load(f) or {}
                    sub_cfg["run"] = False
                    tmp = tempfile.NamedTemporaryFile(
                        mode="w", suffix=".yml", delete=False
                    )
                    os.chmod(tmp.name, 0o600)  # Restrict to owner-only read/write
                    yaml.safe_dump(sub_cfg, tmp)
                    tmp.close()
                    self._temp_files.append(Path(tmp.name))
                    config_dict[cfg_key] = tmp.name
                except Exception:
                    pass  # Best-effort; pipeline will use original config

    def _apply_stage_config_overrides(self, config_dict: dict[str, Any]) -> None:
        """Apply per-run stage config overrides via temporary YAML files."""
        for config_type, override_payload in self.config_overrides.items():
            if config_type == "main":
                continue

            config_key = self.CONFIG_TYPE_TO_MAIN_KEY.get(config_type)
            if not config_key:
                continue

            validation = ConfigValidator.validate(config_type, override_payload)
            if not validation["valid"]:
                raise ValueError(
                    f"Invalid override config for '{config_type}': "
                    + "; ".join(validation["errors"])
                )

            tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".yml", delete=False)
            os.chmod(tmp.name, 0o600)  # Restrict to owner-only read/write
            yaml.safe_dump(override_payload, tmp, sort_keys=False)
            tmp.close()

            self._temp_files.append(Path(tmp.name))
            config_dict[config_key] = tmp.name

    def _emit_stage_transition(self, tui_stage: str) -> None:
        """Emit stage_complete/stage_start notifications on stage change."""
        if self.previous_stage and self.previous_stage != tui_stage:
            self._notify("stage_complete", {"stage": self.previous_stage})
            self._log("info", f"Stage completed: {self.previous_stage}")

        if self.current_stage != tui_stage:
            self._notify("stage_start", {"stage": tui_stage})
            self._log("info", f"Starting stage: {tui_stage}")

        self.previous_stage = tui_stage
        self.current_stage = tui_stage

    def _progress_callback(self, event: dict) -> None:
        """Pipeline progress callback; raises InterruptedError on cancel."""
        if self.cancelled:
            raise InterruptedError("Pipeline cancelled")

        event_type = event.get("type")
        stage = str(event.get("stage", ""))
        tui_stage = self._map_stage_name(stage)
        self._emit_stage_transition(tui_stage)

        if event_type == "stage_progress":
            current = int(event.get("current", 0) or 0)
            total = int(event.get("total", 0) or 0)
            progress_pct = int((current / total) * 100) if total > 0 else 0
            progress_pct = max(0, min(100, progress_pct))
            message = event.get("message") or f"Running: {tui_stage}"
            self.progress.setdefault(tui_stage, 0)
            self.progress[tui_stage] = progress_pct
            self._notify(
                "progress",
                {
                    "stage": tui_stage,
                    "current": progress_pct,
                    "total": 100,
                    "message": message,
                },
            )
            return

        if event_type == "stage_start":
            self.progress[tui_stage] = 0
            message = event.get("message") or f"Starting: {tui_stage}"
            self._notify(
                "progress",
                {
                    "stage": tui_stage,
                    "current": 0,
                    "total": 100,
                    "message": message,
                },
            )
            return

        if event_type == "stage_complete":
            self.progress[tui_stage] = 100
            message = event.get("message") or f"Completed: {tui_stage}"
            self._notify(
                "progress",
                {
                    "stage": tui_stage,
                    "current": 100,
                    "total": 100,
                    "message": message,
                },
            )
            return

    def _run(self):
        """Execute the pipeline."""
        try:
            from hedgehog.configs.logger import load_config, logger
            from hedgehog.pipeline import calculate_metrics
            from hedgehog.utils.data_prep import prepare_input_data
            from hedgehog.utils.mol_index import assign_mol_idx

            config_handler = getattr(self.server, "config_handler", None)
            if not isinstance(config_handler, ConfigHandler):
                config_handler = ConfigHandler(self.server)

            main_override = self.config_overrides.get("main")
            if main_override is not None:
                config_dict = copy.deepcopy(main_override)
            else:
                config_path = config_handler.get_config_path("main")
                config_dict = load_config(str(config_path)) or {}

            main_validation = ConfigValidator.validate("main", config_dict)
            if not main_validation["valid"]:
                raise ValueError(
                    "Invalid main config override: "
                    + "; ".join(main_validation["errors"])
                )

            config_dict.update(config_handler.get_runtime_config_overrides())
            self._apply_stage_config_overrides(config_dict)
            self._disable_unrequested_stages(config_dict)

            self._log("info", "Preparing input data...")
            data = prepare_input_data(config_dict, logger)

            if "mol_idx" not in data.columns or data["mol_idx"].isna().all():
                folder_to_save = Path(config_dict.get("folder_to_save", "results"))
                data = assign_mol_idx(data, run_base=folder_to_save, logger=logger)

            self._log("info", f"Loaded {len(data)} molecules")
            self._log("info", "Starting pipeline execution...")

            success = calculate_metrics(data, config_dict, self._progress_callback)

            if self.current_stage and not self.cancelled:
                self._notify("stage_complete", {"stage": self.current_stage})
                self._log("info", f"Stage completed: {self.current_stage}")

            if not self.cancelled:
                self._notify(
                    "complete",
                    {
                        "job_id": self.job_id,
                        "success": success,
                        "results": {"molecules_processed": len(data)},
                    },
                )

        except InterruptedError:
            self._log("warn", "Pipeline was cancelled")
        except Exception as e:
            self._notify("error", {"message": str(e)})
            self._log("error", f"Pipeline error: {e}")
        finally:
            for tmp_path in self._temp_files:
                tmp_path.unlink(missing_ok=True)
            self._temp_files.clear()


_MAX_COMPLETED_JOBS = 50


class PipelineHandler:
    """Handler for pipeline-related RPC methods."""

    TUI_STAGE_TO_CONFIG_TYPE = {
        "mol_prep": "mol_prep",
        "descriptors": "descriptors",
        "struct_filters": "filters",
        "synthesis": "synthesis",
        "docking": "docking",
        "docking_filters": "docking_filters",
    }

    def __init__(self, server: "JsonRpcServer"):
        self.server = server
        self.jobs: dict[str, PipelineJob] = {}
        self._lock = threading.Lock()

    def _cleanup_completed_jobs(self) -> None:
        """Remove oldest completed jobs when exceeding the limit."""
        with self._lock:
            completed = [
                jid
                for jid, job in self.jobs.items()
                if not job._thread or not job._thread.is_alive()
            ]
            excess = len(completed) - _MAX_COMPLETED_JOBS
            if excess > 0:
                for jid in completed[:excess]:
                    del self.jobs[jid]

    def start_pipeline(
        self,
        stages: list[str],
        config_overrides: dict[str, Any] | None = None,
    ) -> str:
        """Start a new pipeline run."""
        self._cleanup_completed_jobs()
        job_id = str(uuid.uuid4())[:8]
        job = PipelineJob(
            job_id, stages, self.server, config_overrides=config_overrides
        )
        with self._lock:
            self.jobs[job_id] = job
        job.start()
        return job_id

    @staticmethod
    def _make_check(
        code: str,
        level: str,
        message: str,
        stage: str | None = None,
        field: str | None = None,
    ) -> dict[str, Any]:
        check = {
            "code": code,
            "level": level,
            "message": message,
        }
        if stage:
            check["stage"] = stage
        if field:
            check["field"] = field
        return check

    @staticmethod
    def _estimate_runtime(stages: list[str], molecule_count: int | None) -> str:
        """Estimate runtime bucket from selected stages and input size."""
        if molecule_count is None:
            return "unknown"

        score = 0
        if molecule_count >= 10000:
            score += 3
        elif molecule_count >= 3000:
            score += 2
        elif molecule_count >= 1000:
            score += 1

        if "synthesis" in stages:
            score += 1
        if "docking" in stages:
            score += 2
        if "docking_filters" in stages:
            score += 1

        if score <= 1:
            return "short"
        if score <= 3:
            return "medium"
        return "long"

    def preflight_pipeline(
        self,
        stages: list[str],
        config_overrides: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        """Run preflight checks for selected pipeline stages."""
        selected_stages = [s for s in stages if s in self.TUI_STAGE_TO_CONFIG_TYPE]
        if not selected_stages:
            raise ValueError("Select at least one valid stage")
        resolved_overrides = PipelineJob._validate_override_payload(config_overrides)

        config_handler = getattr(self.server, "config_handler", None)
        if not isinstance(config_handler, ConfigHandler):
            config_handler = ConfigHandler(self.server)

        validation_handler = getattr(self.server, "validation_handler", None)
        if not isinstance(validation_handler, ValidationHandler):
            validation_handler = ValidationHandler(self.server)

        files_handler = getattr(self.server, "files_handler", None)
        if not isinstance(files_handler, FilesHandler):
            files_handler = FilesHandler(self.server)

        result: dict[str, Any] = {
            "valid": True,
            "molecule_count": None,
            "estimated_runtime": "unknown",
            "checks": [],
            "stage_reports": [],
        }

        def add_global(
            code: str, level: str, message: str, field: str | None = None
        ) -> None:
            result["checks"].append(self._make_check(code, level, message, field=field))

        if "main" in resolved_overrides:
            main_config = copy.deepcopy(resolved_overrides["main"])
        else:
            try:
                main_config = config_handler.load_config("main")
            except Exception as exc:
                add_global(
                    "MAIN_CONFIG_LOAD_FAILED",
                    "error",
                    f"Failed to load main config: {exc}",
                )
                result["valid"] = False
                return result

        main_validation = ConfigValidator.validate("main", main_config)
        for msg in main_validation.get("errors", []):
            add_global("MAIN_CONFIG_INVALID", "error", msg)
        for msg in main_validation.get("warnings", []):
            add_global("MAIN_CONFIG_WARNING", "warning", msg)

        input_path = main_config.get("generated_mols_path")
        if input_path:
            try:
                input_validation = validation_handler.validate_input_file(input_path)
                for msg in input_validation.get("errors", []):
                    add_global(
                        "MAIN_INPUT_INVALID", "error", msg, field="generated_mols_path"
                    )
                for msg in input_validation.get("warnings", []):
                    add_global(
                        "MAIN_INPUT_WARNING",
                        "warning",
                        msg,
                        field="generated_mols_path",
                    )
            except Exception as exc:
                add_global(
                    "MAIN_INPUT_CHECK_FAILED",
                    "error",
                    f"Failed to validate input file: {exc}",
                    field="generated_mols_path",
                )

            try:
                count_result = files_handler.count_molecules(input_path)
                if "error" in count_result:
                    add_global(
                        "MAIN_COUNT_WARNING",
                        "warning",
                        f"Could not count molecules: {count_result['error']}",
                        field="generated_mols_path",
                    )
                else:
                    result["molecule_count"] = int(count_result.get("count", 0))
            except Exception as exc:
                add_global(
                    "MAIN_COUNT_FAILED",
                    "warning",
                    f"Could not count molecules: {exc}",
                    field="generated_mols_path",
                )

        output_path = main_config.get("folder_to_save")
        if output_path:
            try:
                output_validation = validation_handler.validate_output_directory(
                    output_path
                )
                for msg in output_validation.get("errors", []):
                    add_global(
                        "MAIN_OUTPUT_INVALID", "error", msg, field="folder_to_save"
                    )
                for msg in output_validation.get("warnings", []):
                    add_global(
                        "MAIN_OUTPUT_WARNING", "warning", msg, field="folder_to_save"
                    )
            except Exception as exc:
                add_global(
                    "MAIN_OUTPUT_CHECK_FAILED",
                    "error",
                    f"Failed to validate output directory: {exc}",
                    field="folder_to_save",
                )

        for stage in selected_stages:
            stage_checks: list[dict[str, Any]] = []
            config_type = self.TUI_STAGE_TO_CONFIG_TYPE[stage]

            def add_stage(
                code: str,
                level: str,
                message: str,
                field: str | None = None,
                _stage_checks: list[dict[str, Any]] = stage_checks,
                _stage_name: str = stage,
            ) -> None:
                _stage_checks.append(
                    self._make_check(
                        code, level, message, stage=_stage_name, field=field
                    )
                )

            if config_type in resolved_overrides:
                stage_config = copy.deepcopy(resolved_overrides[config_type])
            else:
                try:
                    stage_config = config_handler.load_config(config_type)
                except Exception as exc:
                    add_stage(
                        "STAGE_CONFIG_LOAD_FAILED",
                        "error",
                        f"Failed to load config: {exc}",
                    )
                    result["stage_reports"].append(
                        {"stage": stage, "status": "error", "checks": stage_checks}
                    )
                    continue

            stage_validation = ConfigValidator.validate(config_type, stage_config)
            for msg in stage_validation.get("errors", []):
                add_stage("STAGE_CONFIG_INVALID", "error", msg)
            for msg in stage_validation.get("warnings", []):
                add_stage("STAGE_CONFIG_WARNING", "warning", msg)

            if stage == "docking" and stage_config.get("run", False):
                receptor = stage_config.get("receptor_pdb")
                if not receptor:
                    add_stage(
                        "DOCKING_RECEPTOR_REQUIRED",
                        "error",
                        "receptor_pdb is required when docking is enabled",
                        field="receptor_pdb",
                    )
                else:
                    try:
                        receptor_validation = validation_handler.validate_receptor_pdb(
                            receptor
                        )
                        for msg in receptor_validation.get("errors", []):
                            add_stage(
                                "DOCKING_RECEPTOR_INVALID",
                                "error",
                                msg,
                                field="receptor_pdb",
                            )
                        for msg in receptor_validation.get("warnings", []):
                            add_stage(
                                "DOCKING_RECEPTOR_WARNING",
                                "warning",
                                msg,
                                field="receptor_pdb",
                            )
                    except Exception as exc:
                        add_stage(
                            "DOCKING_RECEPTOR_CHECK_FAILED",
                            "error",
                            f"Failed to validate receptor: {exc}",
                            field="receptor_pdb",
                        )

            if stage == "docking_filters" and stage_config.get("run", False):
                run_after_docking = bool(stage_config.get("run_after_docking", True))
                input_sdf = stage_config.get("input_sdf")
                receptor_pdb = stage_config.get("receptor_pdb")

                if not run_after_docking and not input_sdf:
                    add_stage(
                        "DOCKING_FILTERS_INPUT_REQUIRED",
                        "error",
                        "input_sdf is required when run_after_docking is false",
                        field="input_sdf",
                    )

                if input_sdf:
                    try:
                        input_validation = validation_handler.validate_input_file(
                            str(input_sdf)
                        )
                        for msg in input_validation.get("errors", []):
                            add_stage(
                                "DOCKING_FILTERS_INPUT_INVALID",
                                "error",
                                msg,
                                field="input_sdf",
                            )
                        for msg in input_validation.get("warnings", []):
                            add_stage(
                                "DOCKING_FILTERS_INPUT_WARNING",
                                "warning",
                                msg,
                                field="input_sdf",
                            )
                    except Exception as exc:
                        add_stage(
                            "DOCKING_FILTERS_INPUT_CHECK_FAILED",
                            "error",
                            f"Failed to validate input_sdf: {exc}",
                            field="input_sdf",
                        )

                if receptor_pdb:
                    try:
                        receptor_validation = validation_handler.validate_receptor_pdb(
                            str(receptor_pdb)
                        )
                        for msg in receptor_validation.get("errors", []):
                            add_stage(
                                "DOCKING_FILTERS_RECEPTOR_INVALID",
                                "error",
                                msg,
                                field="receptor_pdb",
                            )
                        for msg in receptor_validation.get("warnings", []):
                            add_stage(
                                "DOCKING_FILTERS_RECEPTOR_WARNING",
                                "warning",
                                msg,
                                field="receptor_pdb",
                            )
                    except Exception as exc:
                        add_stage(
                            "DOCKING_FILTERS_RECEPTOR_CHECK_FAILED",
                            "error",
                            f"Failed to validate receptor_pdb: {exc}",
                            field="receptor_pdb",
                        )

            if any(check["level"] == "error" for check in stage_checks):
                status = "error"
            elif any(check["level"] == "warning" for check in stage_checks):
                status = "warning"
            else:
                status = "ok"

            result["stage_reports"].append(
                {"stage": stage, "status": status, "checks": stage_checks}
            )

        all_stage_checks = [
            check
            for report in result["stage_reports"]
            for check in report.get("checks", [])
        ]
        all_checks = [*result["checks"], *all_stage_checks]

        result["valid"] = not any(check["level"] == "error" for check in all_checks)
        result["estimated_runtime"] = self._estimate_runtime(
            selected_stages, result["molecule_count"]
        )
        return result

    def get_progress(self, job_id: str) -> dict[str, Any]:
        """Get progress of a running pipeline."""
        with self._lock:
            if job_id not in self.jobs:
                raise ValueError(f"Job not found: {job_id}")
            job = self.jobs[job_id]
        return {
            "job_id": job_id,
            "stages": job.stages,
            "current_stage": job.current_stage,
            "progress": job.progress,
            "cancelled": job.cancelled,
        }

    def cancel_pipeline(self, job_id: str) -> bool:
        """Cancel a running pipeline."""
        with self._lock:
            if job_id not in self.jobs:
                raise ValueError(f"Job not found: {job_id}")
            job = self.jobs[job_id]
        job.cancel()
        return True
