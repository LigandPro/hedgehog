"""Job history handler for TUI backend."""

import json
import os
import threading
from copy import deepcopy
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..server import JsonRpcServer

HISTORY_SEED_ENV = "HEDGEHOG_TUI_HISTORY_SEED"
SEEDED_HISTORY = [
    {
        "id": "audit-run-001",
        "name": "Audit Run: Full Pipeline",
        "startTime": "2025-01-15T10:00:00",
        "endTime": "2025-01-15T10:11:32",
        "status": "completed",
        "config": {
            "inputPath": "src/hedgehog/configs/examples/moses_1000.csv",
            "outputPath": "results/audit_run_001",
            "stages": ["descriptors", "struct_filters", "synthesis"],
        },
        "results": {
            "inputMolecules": 1000,
            "outputMolecules": 742,
            "durationSeconds": 692,
        },
        "error": None,
    },
    {
        "id": "audit-run-002",
        "name": "Audit Run: Docking",
        "startTime": "2025-01-14T14:25:00",
        "endTime": "2025-01-14T14:31:08",
        "status": "error",
        "config": {
            "inputPath": "src/hedgehog/configs/examples/moses_1000.csv",
            "outputPath": "results/audit_run_002",
            "stages": ["docking"],
        },
        "results": None,
        "error": "Docking failed: receptor file validation error",
    },
    {
        "id": "audit-run-003",
        "name": "Audit Run: Retrosynthesis",
        "startTime": "2025-01-13T09:40:00",
        "endTime": None,
        "status": "running",
        "config": {
            "inputPath": "src/hedgehog/configs/examples/target_mols.csv",
            "outputPath": "results/audit_run_003",
            "stages": ["retrosynthesis"],
        },
        "results": None,
        "error": None,
    },
]


class HistoryHandler:
    """Handler for job history RPC methods."""

    def __init__(self, server: "JsonRpcServer"):
        self.server = server
        self.history_file = Path.home() / ".hedgehog" / "job_history.json"
        self._lock = threading.Lock()
        self._ensure_history_dir()

    def _ensure_history_dir(self):
        """Ensure the history directory exists."""
        self.history_file.parent.mkdir(parents=True, exist_ok=True)
        if not self.history_file.exists():
            self.history_file.write_text("[]")

    def _load_history(self) -> list[dict]:
        """Load history from file."""
        try:
            return json.loads(self.history_file.read_text())
        except (json.JSONDecodeError, FileNotFoundError):
            return []

    def _save_history(self, history: list[dict]):
        """Save history to file."""
        self.history_file.write_text(json.dumps(history, indent=2, default=str))

    def get_job_history(self, limit: int = 50) -> list[dict]:
        """Get job history."""
        with self._lock:
            history = self._load_history()
            if not history and self._history_seed_enabled():
                history = deepcopy(SEEDED_HISTORY)
            return history[:limit]

    def add_job(
        self,
        job_id: str,
        name: str | None = None,
        input_path: str = "",
        output_path: str = "",
        stages: list[str] | None = None,
    ) -> dict:
        """Add a new job to history."""
        with self._lock:
            history = self._load_history()

            job = {
                "id": job_id,
                "name": name or f"Job {job_id}",
                "startTime": datetime.now().isoformat(),
                "endTime": None,
                "status": "running",
                "config": {
                    "inputPath": input_path,
                    "outputPath": output_path,
                    "stages": stages or [],
                },
                "results": None,
                "error": None,
            }

            history.insert(0, job)
            history = history[:100]  # Keep last 100 jobs
            self._save_history(history)

            return job

    def update_job(
        self,
        job_id: str,
        status: str | None = None,
        results: dict | None = None,
        error: str | None = None,
    ) -> dict | None:
        """Update a job in history."""
        with self._lock:
            history = self._load_history()

            for job in history:
                if job["id"] == job_id:
                    if status is not None:
                        job["status"] = status
                        if status in ("completed", "cancelled", "error"):
                            job["endTime"] = datetime.now().isoformat()
                    if results is not None:
                        job["results"] = results
                    if error is not None:
                        job["error"] = error

                    self._save_history(history)
                    return job

            return None

    def delete_job(self, job_id: str) -> bool:
        """Delete a job from history."""
        with self._lock:
            history = self._load_history()
            original_len = len(history)

            history = [j for j in history if j["id"] != job_id]

            if len(history) < original_len:
                self._save_history(history)
                return True

            return False

    def _history_seed_enabled(self) -> bool:
        """Return whether deterministic seeded history mode is enabled."""
        value = os.getenv(HISTORY_SEED_ENV, "")
        return value.strip().lower() in {"1", "true", "yes", "on"}
