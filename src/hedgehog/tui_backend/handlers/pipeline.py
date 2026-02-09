"""Pipeline execution handler for TUI backend."""

import threading
import uuid
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from ..server import JsonRpcServer


def _find_project_root() -> Path:
    """Find project root by looking for pyproject.toml."""
    current = Path.cwd()
    for parent in [current] + list(current.parents):
        if (parent / "pyproject.toml").exists():
            return parent
    return current


class PipelineJob:
    """Represents a running pipeline job."""

    # Map internal stage names to TUI stage names
    STAGE_NAME_MAP = {
        "struct_ini_filters": "struct_filters",
        "descriptors": "descriptors",
        "struct_filters": "struct_filters",
        "synthesis": "synthesis",
        "docking": "docking",
        "docking_filters": "docking_filters",
        "final_descriptors": "descriptors",
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

    def _run(self):
        """Execute the pipeline."""
        try:
            # Import hedgehog modules
            from hedgehog.configs.logger import load_config, logger
            from hedgehog.pipeline import calculate_metrics
            from hedgehog.utils.data_prep import prepare_input_data
            from hedgehog.utils.mol_index import assign_mol_idx

            # Load config
            config_path = _find_project_root() / "src/hedgehog/configs/config.yml"
            config_dict = load_config(str(config_path))

            # Set up progress callback for real pipeline
            def progress_callback(stage: str, current: int, total: int):
                if self.cancelled:
                    raise InterruptedError("Pipeline cancelled")

                # Map stage name
                tui_stage = self._map_stage_name(stage)

                # If stage changed, send stage_complete for previous and stage_start for new
                if self.previous_stage and self.previous_stage != tui_stage:
                    prev_tui_stage = self._map_stage_name(self.previous_stage)
                    self.server.send_notification(
                        "stage_complete",
                        {
                            "stage": prev_tui_stage,
                        },
                    )
                    self.server.send_notification(
                        "log",
                        {
                            "level": "info",
                            "message": f"Stage completed: {prev_tui_stage}",
                        },
                    )

                if self.current_stage != tui_stage:
                    self.server.send_notification(
                        "stage_start",
                        {
                            "stage": tui_stage,
                        },
                    )
                    self.server.send_notification(
                        "log",
                        {
                            "level": "info",
                            "message": f"Starting stage: {tui_stage}",
                        },
                    )

                self.previous_stage = stage
                self.current_stage = tui_stage

                # Send progress update - use percentage within stage context
                # Since we don't have fine-grained progress, use current/total as rough indicator
                progress_pct = int((current / total) * 100) if total > 0 else 0
                self.server.send_notification(
                    "progress",
                    {
                        "stage": tui_stage,
                        "current": progress_pct,
                        "total": 100,
                        "message": f"Running stage {current}/{total}: {tui_stage}",
                    },
                )

            # Prepare data
            self.server.send_notification(
                "log",
                {
                    "level": "info",
                    "message": "Preparing input data...",
                },
            )

            data = prepare_input_data(config_dict, logger)

            if "mol_idx" not in data.columns or data["mol_idx"].isna().all():
                folder_to_save = Path(config_dict.get("folder_to_save", "results"))
                data = assign_mol_idx(data, run_base=folder_to_save, logger=logger)

            self.server.send_notification(
                "log",
                {
                    "level": "info",
                    "message": f"Loaded {len(data)} molecules",
                },
            )

            # Run actual pipeline with progress callback
            self.server.send_notification(
                "log",
                {
                    "level": "info",
                    "message": "Starting pipeline execution...",
                },
            )

            success = calculate_metrics(data, config_dict, progress_callback)

            # Mark last stage as complete
            if self.current_stage and not self.cancelled:
                self.server.send_notification(
                    "stage_complete",
                    {
                        "stage": self.current_stage,
                    },
                )
                self.server.send_notification(
                    "log",
                    {
                        "level": "info",
                        "message": f"Stage completed: {self.current_stage}",
                    },
                )

            if not self.cancelled:
                self.server.send_notification(
                    "complete",
                    {
                        "job_id": self.job_id,
                        "success": success,
                        "results": {
                            "molecules_processed": len(data),
                        },
                    },
                )

        except InterruptedError:
            self.server.send_notification(
                "log",
                {
                    "level": "warn",
                    "message": "Pipeline was cancelled",
                },
            )
        except Exception as e:
            self.server.send_notification(
                "error",
                {
                    "message": str(e),
                },
            )
            self.server.send_notification(
                "log",
                {
                    "level": "error",
                    "message": f"Pipeline error: {e}",
                },
            )


_MAX_COMPLETED_JOBS = 50


class PipelineHandler:
    """Handler for pipeline-related RPC methods."""

    def __init__(self, server: "JsonRpcServer"):
        self.server = server
        self.jobs: dict[str, PipelineJob] = {}

    def _cleanup_completed_jobs(self) -> None:
        """Remove oldest completed jobs when exceeding the limit."""
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
        self.jobs[job_id] = job
        job.start()
        return job_id

    def get_progress(self, job_id: str) -> dict[str, Any]:
        """Get progress of a running pipeline."""
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
        if job_id not in self.jobs:
            raise ValueError(f"Job not found: {job_id}")

        self.jobs[job_id].cancel()
        return True
