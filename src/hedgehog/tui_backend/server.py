"""JSON-RPC server for TUI backend communication."""

import json
import sys
import threading
from typing import Any

from .handlers.config import ConfigHandler
from .handlers.files import FilesHandler
from .handlers.history import HistoryHandler
from .handlers.pipeline import PipelineHandler
from .handlers.validation import ValidationHandler


class JsonRpcServer:
    """JSON-RPC 2.0 server over stdio."""

    def __init__(self):
        self.handlers: dict[str, Any] = {}
        self.request_id = 0
        self._running = False
        self._lock = threading.Lock()

        # Register handlers
        self.config_handler = ConfigHandler(self)
        self.files_handler = FilesHandler(self)
        self.history_handler = HistoryHandler(self)
        self.pipeline_handler = PipelineHandler(self)
        self.validation_handler = ValidationHandler(self)

        self._register_methods()

    def _register_methods(self):
        """Register all RPC methods."""
        # Config methods
        self.handlers["load_config"] = self.config_handler.load_config
        self.handlers["save_config"] = self.config_handler.save_config
        self.handlers["validate_config"] = self.config_handler.validate_config

        # File methods
        self.handlers["list_files"] = self.files_handler.list_files
        self.handlers["list_directory"] = self.files_handler.list_directory
        self.handlers["count_molecules"] = self.files_handler.count_molecules

        # Pipeline methods
        self.handlers["start_pipeline"] = self.pipeline_handler.start_pipeline
        self.handlers["get_progress"] = self.pipeline_handler.get_progress
        self.handlers["cancel_pipeline"] = self.pipeline_handler.cancel_pipeline

        # Validation methods
        self.handlers["validate_input_file"] = (
            self.validation_handler.validate_input_file
        )
        self.handlers["validate_receptor_pdb"] = (
            self.validation_handler.validate_receptor_pdb
        )
        self.handlers["validate_output_directory"] = (
            self.validation_handler.validate_output_directory
        )
        self.handlers["validate_config_data"] = self.validation_handler.validate_config

        # History methods
        self.handlers["get_job_history"] = self.history_handler.get_job_history
        self.handlers["add_job"] = self.history_handler.add_job
        self.handlers["update_job"] = self.history_handler.update_job
        self.handlers["delete_job"] = self.history_handler.delete_job

    def send_response(
        self, request_id: int, result: Any = None, error: dict | None = None
    ):
        """Send JSON-RPC response."""
        response = {
            "jsonrpc": "2.0",
            "id": request_id,
        }
        if error:
            response["error"] = error
        else:
            response["result"] = result

        self._send_message(response)

    def send_notification(self, method: str, params: dict):
        """Send JSON-RPC notification (no id)."""
        notification = {
            "jsonrpc": "2.0",
            "method": method,
            "params": params,
        }
        self._send_message(notification)

    def _send_message(self, message: dict):
        """Send a message to stdout."""
        with self._lock:
            line = json.dumps(message) + "\n"
            sys.stdout.write(line)
            sys.stdout.flush()

    def handle_request(self, request: dict) -> None:
        """Handle a JSON-RPC request."""
        request_id = request.get("id")
        method = request.get("method", "")
        params = request.get("params", {})

        if method not in self.handlers:
            self.send_response(
                request_id,
                error={"code": -32601, "message": f"Method not found: {method}"},
            )
            return

        try:
            handler = self.handlers[method]
            result = handler(**params) if params else handler()
            self.send_response(request_id, result=result)
        except Exception as e:
            self.send_response(request_id, error={"code": -32000, "message": str(e)})

    def run(self):
        """Run the server, reading from stdin."""
        self._running = True

        # Signal ready
        sys.stderr.write("HEDGEHOG_TUI_READY\n")
        sys.stderr.flush()

        for line in sys.stdin:
            if not self._running:
                break

            line = line.strip()
            if not line:
                continue

            try:
                request = json.loads(line)
                # Handle request in separate thread for non-blocking
                threading.Thread(
                    target=self.handle_request, args=(request,), daemon=True
                ).start()
            except json.JSONDecodeError as e:
                # Send parse error
                self.send_response(
                    None, error={"code": -32700, "message": f"Parse error: {e}"}
                )

    def stop(self):
        """Stop the server."""
        self._running = False


def main() -> int:
    """Main entry point."""
    server = JsonRpcServer()
    try:
        server.run()
        return 0
    except KeyboardInterrupt:
        return 0
    except Exception as e:
        sys.stderr.write(f"Server error: {e}\n")
        return 1


if __name__ == "__main__":
    sys.exit(main())
