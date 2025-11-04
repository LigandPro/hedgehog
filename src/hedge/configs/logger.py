import logging
import os
from datetime import datetime

import yaml
from rich.console import Console
from rich.logging import RichHandler


class LoggerSingleton:
    _instance = None
    _logger = None
    _console = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(LoggerSingleton, cls).__new__(cls)
            cls._console = Console()
        return cls._instance

    def get_logger(self, name="run"):
        if self._logger is None:
            self._logger = self._setup_logger(name)
        return self._logger

    def _setup_logger(self, name="run"):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        config_path = os.path.join(script_dir, "config.yml")
        config = yaml.safe_load(open(config_path))
        folder_to_save = config["folder_to_save"]
        logs_dir = f"{folder_to_save}/"
        os.makedirs(logs_dir, exist_ok=True)

        logger = logging.getLogger(name)
        logger.setLevel(logging.INFO)
        logger.handlers = []

        # Rich console handler - beautiful colored output
        console_handler = RichHandler(rich_tracebacks=True,
                                      tracebacks_show_locals=True,
                                      markup=True,
                                      show_time=True,
                                      show_level=True,
                                      show_path=False,
                                      console=self._console
                                     )
        console_handler.setLevel(logging.INFO)
        logger.addHandler(console_handler)

        # Rich file handler - plain text (no colors/markup) using native rich solution
        current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = os.path.join(logs_dir, f"{name}_{current_time}.log")
        file_console = Console(file=open(log_file, "w", encoding="utf-8"), no_color=True, width=120)
        file_handler = RichHandler(rich_tracebacks=True,
                                   tracebacks_show_locals=True,
                                   markup=True,
                                   show_time=True,
                                   show_level=True,
                                   show_path=False,
                                   console=file_console
                                  )
        file_handler.setLevel(logging.INFO)
        logger.addHandler(file_handler)

        return logger


def setup_logger(name="run"):
    return LoggerSingleton().get_logger(name)


def get_logger():
    return LoggerSingleton().get_logger()


# Lazy logger that initializes only when accessed
class LazyLogger:
    def __getattr__(self, name):
        logger = get_logger()
        return getattr(logger, name)


logger = LazyLogger()


def load_config(config_path="./src/hedge/configs/config.yml"):
    """Load YAML configuration file.
    
    Args:
        config_path: Path to the YAML configuration file
        
    Returns
    -------
        Dictionary containing configuration data
        
    Raises
    ------
        FileNotFoundError: If configuration file doesn't exist
        yaml.YAMLError: If YAML parsing fails
    """
    try:
        with open(config_path) as file:
            return yaml.safe_load(file)
    except FileNotFoundError:
        logger.error(f"Configuration file not found: {config_path}")
        raise
    except yaml.YAMLError as e:
        logger.error(f"Error parsing YAML configuration: {e}")
        raise

