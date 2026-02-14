from __future__ import annotations

import logging
import re
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml
from rich.console import Console
from rich.logging import RichHandler

# Pattern to strip Rich markup tags like [#B29EEE], [bold], [/bold], etc.
# Excludes standard log levels: [INFO], [WARNING], [ERROR], [DEBUG], [CRITICAL]
_RICH_MARKUP_RE = re.compile(
    r"\[/?(?!INFO\]|WARNING\]|ERROR\]|DEBUG\]|CRITICAL\])[^\]]+\]"
)

CONFIG_PATH = Path(__file__).resolve().parent / "config.yml"


class LoggerSingleton:
    """Singleton responsible for configuring project logging."""

    _instance: LoggerSingleton | None = None
    _logger: logging.Logger | None = None
    _console: Console | None = None
    _log_directory: Path | None = None
    _file_handler_added: bool = False
    _file_handler: logging.Handler | None = None
    _lock = threading.Lock()

    def __new__(cls):
        """Create or reuse singleton instance (thread-safe)."""
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = super().__new__(cls)
                    cls._console = Console()
        return cls._instance

    def get_logger(self, name="run"):
        """Return shared logger instance configured with Rich handlers."""
        if self._logger is None:
            self._logger = self._setup_logger(name)
        return self._logger

    @property
    def console(self) -> Console:
        """Return the shared Rich Console used by the logging handler."""
        return self._console

    def configure_log_directory(self, folder_to_save: Path) -> None:
        """Configure the log directory and add file handler.

        This method should be called after CLI overrides are applied to ensure
        logs are written to the correct directory.

        Args:
            folder_to_save: Directory where log files should be written.
        """
        new_dir = Path(folder_to_save).resolve()

        # If we're reusing the same Python process for multiple runs, allow
        # reconfiguration so each run writes logs into its own folder.
        if self._log_directory is not None and new_dir != self._log_directory.resolve():
            if self._logger is not None and self._file_handler is not None:
                try:
                    self._logger.removeHandler(self._file_handler)
                    self._file_handler.close()
                except Exception:
                    pass
            self._file_handler_added = False
            self._file_handler = None

        self._log_directory = new_dir
        self._log_directory.mkdir(parents=True, exist_ok=True)

        # If the logger is already initialized, attach the file handler now.
        # If not, _setup_logger() will attach it lazily once the logger is created.
        if self._logger is not None:
            self._ensure_file_handler()

    def _setup_logger(self, name="run"):
        """Configure console handler for Rich logging output.

        File handler is added later via configure_log_directory() to ensure
        CLI-overridden folder_to_save is respected.
        """
        logger = logging.getLogger(name)
        logger.setLevel(logging.INFO)
        logger.handlers = []

        # Rich console handler - beautiful colored output
        console_handler = RichHandler(
            rich_tracebacks=True,
            tracebacks_show_locals=True,
            markup=True,
            show_time=True,
            show_level=True,
            show_path=False,
            console=self._console,
        )
        console_handler.setLevel(logging.INFO)
        logger.addHandler(console_handler)

        # If a log directory was configured before the logger existed, attach the
        # file handler now so each run gets a persistent log file.
        self._logger = logger
        self._ensure_file_handler()

        return logger

    def _ensure_file_handler(self) -> None:
        """Attach a file handler if a log directory is configured."""
        if self._file_handler_added:
            return
        if self._logger is None or self._log_directory is None:
            return

        current_time = datetime.now(tz=timezone.utc).strftime("%Y%m%d_%H%M%S")
        log_file = self._log_directory / f"run_{current_time}.log"
        # Use a plain FileHandler for clean, unwrapped logs on disk.
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setFormatter(_PlainTextFormatter())
        file_handler.setLevel(logging.INFO)
        self._logger.addHandler(file_handler)
        self._file_handler = file_handler
        self._file_handler_added = True


class _PlainTextFormatter(logging.Formatter):
    """Formatter that strips Rich markup tags for plain text log files."""

    def __init__(self):
        super().__init__(
            fmt="%(asctime)s [%(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )

    def format(self, record: logging.LogRecord) -> str:
        # Format first, then strip Rich markup from the final result
        result = super().format(record)
        return _RICH_MARKUP_RE.sub("", result)


def setup_logger(name="run"):
    """Initialise logger singleton and return configured logger."""
    return LoggerSingleton().get_logger(name)


def get_logger():
    """Return the shared logger instance."""
    return LoggerSingleton().get_logger()


class LazyLogger:
    """Proxy that lazily resolves attributes on first use."""

    def __getattr__(self, name):
        """Resolve attribute lookups against the underlying logger."""
        logger_instance = get_logger()
        return getattr(logger_instance, name)


logger = LazyLogger()


def load_config(config_path: str | Path = CONFIG_PATH) -> dict[str, Any]:
    """Load YAML configuration file.

    Parameters:
        config_path: Path to the YAML configuration file.

    Returns:
        dict[str, Any]: Parsed configuration dictionary.

    Raises:
        FileNotFoundError: If configuration file doesn't exist.
        yaml.YAMLError: If YAML parsing fails.
    """
    config_path = Path(config_path)
    root_logger = logging.getLogger(__name__)
    try:
        with config_path.open(encoding="utf-8") as file:
            return yaml.safe_load(file)
    except FileNotFoundError:
        root_logger.exception("Configuration file not found: %s", config_path)
        raise
    except yaml.YAMLError:
        root_logger.exception("Error parsing YAML configuration for %s", config_path)
        raise
