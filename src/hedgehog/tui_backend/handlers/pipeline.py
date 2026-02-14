"""Pipeline execution handler for TUI backend."""

import tempfile
import threading
import uuid
from pathlib import Path
from typing import TYPE_CHECKING, Any

import yaml

if TYPE_CHECKING:
    from ..server import JsonRpcServer


def _find_project_root() -> Path:
    """Find project root by looking for pyproject.toml."""
    current = Path.cwd()
    for parent in [current, *current.parents]:
        if (parent / "pyproject.toml").exists():
            return parent
    return current


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
        "descriptors": ["config_descriptors"],
        "struct_filters": ["config_structFilters"],
        "synthesis": ["config_synthesis"],
        "docking": ["config_docking"],
        "docking_filters": ["config_docking_filters"],
    }

    def __init__(self, job_id: str, stages: list[str], server: "JsonRpcServer"):
        self.job_id = job_id
        self.stages = stages
        self.server = server
        self.cancelled = False
        self.current_stage: str | None = None
        self.previous_stage: str | None = None
        self.progress: dict[str, int] = {s: 0 for s in stages}
        self._thread: threading.Thread | None = None

    def start(self):
        """Start the pipeline in a background thread."""
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def cancel(self):
        """Cancel the pipeline."""
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
                    yaml.safe_dump(sub_cfg, tmp)
                    tmp.close()
                    config_dict[cfg_key] = tmp.name
                except Exception:
                    pass  # Best-effort; pipeline will use original config

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

            config_path = _find_project_root() / "src/hedgehog/configs/config.yml"
            config_dict = load_config(str(config_path))
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


_MAX_COMPLETED_JOBS = 50


class PipelineHandler:
    """Handler for pipeline-related RPC methods."""

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

    def start_pipeline(self, stages: list[str]) -> str:
        """Start a new pipeline run."""
        self._cleanup_completed_jobs()
        job_id = str(uuid.uuid4())[:8]
        job = PipelineJob(job_id, stages, self.server)
        with self._lock:
            self.jobs[job_id] = job
        job.start()
        return job_id

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
