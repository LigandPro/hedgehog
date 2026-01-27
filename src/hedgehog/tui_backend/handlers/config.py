"""Config handling for TUI backend."""

from pathlib import Path
from typing import TYPE_CHECKING, Any

import yaml

from ..validators import ConfigValidator

if TYPE_CHECKING:
    from ..server import JsonRpcServer

# Default config paths relative to project root
CONFIG_PATHS = {
    "main": "src/hedgehog/configs/config.yml",
    "descriptors": "src/hedgehog/configs/config_descriptors.yml",
    "filters": "src/hedgehog/configs/config_structFilters.yml",
    "synthesis": "src/hedgehog/configs/config_synthesis.yml",
    "docking": "src/hedgehog/configs/config_docking.yml",
}


class ConfigHandler:
    """Handler for config-related RPC methods."""

    def __init__(self, server: "JsonRpcServer"):
        self.server = server
        self._project_root = self._find_project_root()

    def _find_project_root(self) -> Path:
        """Find project root by looking for pyproject.toml."""
        current = Path.cwd()
        for parent in [current] + list(current.parents):
            if (parent / "pyproject.toml").exists():
                return parent
        return current

    def _get_config_path(self, config_type: str) -> Path:
        """Get the path to a config file."""
        if config_type not in CONFIG_PATHS:
            raise ValueError(f"Unknown config type: {config_type}")
        return self._project_root / CONFIG_PATHS[config_type]

    def load_config(self, config_type: str) -> dict[str, Any]:
        """Load a config file and return its contents."""
        config_path = self._get_config_path(config_type)

        if not config_path.exists():
            raise FileNotFoundError(f"Config file not found: {config_path}")

        with open(config_path) as f:
            return yaml.safe_load(f) or {}

    def save_config(self, config_type: str, data: dict[str, Any]) -> bool:
        """Save data to a config file."""
        config_path = self._get_config_path(config_type)

        # Backup existing config
        if config_path.exists():
            backup_path = config_path.with_suffix(".yml.bak")
            backup_path.write_text(config_path.read_text())

        with open(config_path, "w") as f:
            yaml.dump(
                data, f, default_flow_style=False, sort_keys=False, allow_unicode=True
            )

        return True

    def validate_config(self, config_type: str, data: dict[str, Any]) -> dict[str, Any]:
        """Validate config data using centralized ConfigValidator."""
        return ConfigValidator.validate(config_type, data)
