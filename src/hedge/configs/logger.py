"""Logger utilities with Rich formatting."""

from __future__ import annotations

import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Self

import yaml
from rich.console import Console
from rich.logging import RichHandler

CONFIG_PATH = Path(__file__).resolve().parent / "config.yml"


class LoggerSingleton:
    """Singleton responsible for configuring project logging."""

    _instance: LoggerSingleton | None = None
    _logger: logging.Logger | None = None
    _console: Console | None = None

    def __new__(cls):
        """Create or reuse singleton instance."""
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._console = Console()
        return cls._instance

    def get_logger(self, name="run"):
        """Return shared logger instance configured with Rich handlers."""
        if self._logger is None:
            self._logger = self._setup_logger(name)
        return self._logger

    def _setup_logger(self, name="run"):
        """Configure console and file handlers for Rich logging output."""
        config = load_config(CONFIG_PATH)
        folder_to_save = Path(config["folder_to_save"])
        logs_dir = folder_to_save
        logs_dir.mkdir(parents=True, exist_ok=True)

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

        # Rich file handler - plain text (no colors/markup) using native rich solution
        current_time = datetime.now(tz=timezone.utc).strftime("%Y%m%d_%H%M%S")
        log_file = logs_dir / f"{name}_{current_time}.log"
        file_console = Console(
            file=log_file.open("w", encoding="utf-8"),
            no_color=True,
            width=120,
        )
        file_handler = RichHandler(
            rich_tracebacks=True,
            tracebacks_show_locals=True,
            markup=True,
            show_time=True,
            show_level=True,
            show_path=False,
            console=file_console,
        )
        file_handler.setLevel(logging.INFO)
        logger.addHandler(file_handler)

        return logger


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

    Parameters
    ----------
    config_path:
        Path to the YAML configuration file.

    Returns
    -------
    dict[str, Any]
        Parsed configuration dictionary.

    Raises
    ------
    FileNotFoundError
        If configuration file doesn't exist.
    yaml.YAMLError
        If YAML parsing fails.
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
