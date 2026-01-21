"""Pipeline execution handler for TUI backend."""
import threading
import time
import uuid
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from ..server import JsonRpcServer


class PipelineJob:
    """Represents a running pipeline job."""

    def __init__(self, job_id: str, stages: list[str], server: 'JsonRpcServer'):
        self.job_id = job_id
        self.stages = stages
        self.server = server
        self.cancelled = False
        self.current_stage: str | None = None
        self.progress: dict[str, int] = {s: 0 for s in stages}
        self._thread: threading.Thread | None = None

    def start(self):
        """Start the pipeline in a background thread."""
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def cancel(self):
        """Cancel the pipeline."""
        self.cancelled = True

    def _run(self):
        """Execute the pipeline."""
        try:
            # Import hedgehog modules
            from hedgehog.configs.logger import load_config
            from hedgehog.pipeline import calculate_metrics
            from hedgehog.utils.data_prep import prepare_input_data
            from hedgehog.utils.mol_index import assign_mol_idx

            # Load config
            config_path = Path('src/hedgehog/configs/config.yml')
            config_dict = load_config(str(config_path))
            
            # Set up progress callback
            def progress_callback(stage: str, current: int, total: int, message: str = ''):
                if self.cancelled:
                    raise InterruptedError('Pipeline cancelled')
                self.server.send_notification('progress', {
                    'stage': stage,
                    'current': current,
                    'total': total,
                    'message': message,
                })

            # Prepare data
            self.server.send_notification('log', {
                'level': 'info',
                'message': 'Preparing input data...',
            })
            
            data = prepare_input_data(config_dict, None)
            
            if 'mol_idx' not in data.columns or data['mol_idx'].isna().all():
                folder_to_save = Path(config_dict.get('folder_to_save', 'results'))
                data = assign_mol_idx(data, run_base=folder_to_save, logger=None)

            self.server.send_notification('log', {
                'level': 'info',
                'message': f'Loaded {len(data)} molecules',
            })

            # Run pipeline
            for stage in self.stages:
                if self.cancelled:
                    break
                    
                self.current_stage = stage
                self.server.send_notification('stage_start', {'stage': stage})
                self.server.send_notification('log', {
                    'level': 'info',
                    'message': f'Starting stage: {stage}',
                })
                
                # Simulate progress for now
                # In production, this would call the actual stage
                for i in range(10):
                    if self.cancelled:
                        break
                    time.sleep(0.5)
                    progress_callback(stage, i + 1, 10, f'Processing batch {i + 1}/10')
                
                if not self.cancelled:
                    self.server.send_notification('stage_complete', {'stage': stage})
                    self.progress[stage] = 100

            if not self.cancelled:
                self.server.send_notification('complete', {
                    'job_id': self.job_id,
                    'results': {
                        'stages_completed': self.stages,
                        'molecules_processed': len(data),
                    },
                })
            
        except InterruptedError:
            self.server.send_notification('log', {
                'level': 'warn',
                'message': 'Pipeline was cancelled',
            })
        except Exception as e:
            self.server.send_notification('error', {
                'message': str(e),
            })
            self.server.send_notification('log', {
                'level': 'error',
                'message': f'Pipeline error: {e}',
            })


class PipelineHandler:
    """Handler for pipeline-related RPC methods."""

    def __init__(self, server: 'JsonRpcServer'):
        self.server = server
        self.jobs: dict[str, PipelineJob] = {}

    def start_pipeline(self, stages: list[str]) -> str:
        """Start a new pipeline run."""
        job_id = str(uuid.uuid4())[:8]
        job = PipelineJob(job_id, stages, self.server)
        self.jobs[job_id] = job
        job.start()
        return job_id

    def get_progress(self, job_id: str) -> dict[str, Any]:
        """Get progress of a running pipeline."""
        if job_id not in self.jobs:
            raise ValueError(f'Job not found: {job_id}')
        
        job = self.jobs[job_id]
        return {
            'job_id': job_id,
            'stages': job.stages,
            'current_stage': job.current_stage,
            'progress': job.progress,
            'cancelled': job.cancelled,
        }

    def cancel_pipeline(self, job_id: str) -> bool:
        """Cancel a running pipeline."""
        if job_id not in self.jobs:
            raise ValueError(f'Job not found: {job_id}')
        
        self.jobs[job_id].cancel()
        return True
