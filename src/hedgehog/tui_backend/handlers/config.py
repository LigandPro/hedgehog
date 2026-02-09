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
    "retrosynthesis": "modules/retrosynthesis/aizynthfinder/public/config.yml",
    "docking": "src/hedgehog/configs/config_docking.yml",
    "docking_filters": "src/hedgehog/configs/config_docking_filters.yml",
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
            data = yaml.safe_load(f) or {}

        # Special handling for retrosynthesis config - flatten nested structure
        if config_type == "retrosynthesis":
            data = self._flatten_retrosynthesis_config(data)

        return data

    def _flatten_retrosynthesis_config(self, data: dict[str, Any]) -> dict[str, Any]:
        """Flatten AiZynthFinder nested config to TUI flat structure."""
        flat = {}

        # Expansion models
        expansion = data.get("expansion", {})
        uspto = expansion.get("uspto", [])
        if len(uspto) >= 2:
            flat["expansion_uspto_model"] = uspto[0]
            flat["expansion_uspto_templates"] = uspto[1]
        ringbreaker = expansion.get("ringbreaker", [])
        if len(ringbreaker) >= 2:
            flat["expansion_ringbreaker_model"] = ringbreaker[0]
            flat["expansion_ringbreaker_templates"] = ringbreaker[1]

        # Filter models
        filter_data = data.get("filter", {})
        flat["filter_uspto"] = filter_data.get("uspto", "")

        # Stock databases
        stock = data.get("stock", {})
        flat["stock_zinc"] = stock.get("zinc", "")

        return flat

    def _unflatten_retrosynthesis_config(self, data: dict[str, Any]) -> dict[str, Any]:
        """Convert TUI flat structure back to AiZynthFinder nested config."""
        nested = {
            "expansion": {},
            "filter": {},
            "stock": {},
        }

        # Expansion models
        if data.get("expansion_uspto_model") or data.get("expansion_uspto_templates"):
            nested["expansion"]["uspto"] = [
                data.get("expansion_uspto_model", ""),
                data.get("expansion_uspto_templates", ""),
            ]
        if data.get("expansion_ringbreaker_model") or data.get(
            "expansion_ringbreaker_templates"
        ):
            nested["expansion"]["ringbreaker"] = [
                data.get("expansion_ringbreaker_model", ""),
                data.get("expansion_ringbreaker_templates", ""),
            ]

        # Filter models
        if data.get("filter_uspto"):
            nested["filter"]["uspto"] = data["filter_uspto"]

        # Stock databases
        if data.get("stock_zinc"):
            nested["stock"]["zinc"] = data["stock_zinc"]

        return nested

    def save_config(self, config_type: str, data: dict[str, Any]) -> bool:
        """Save data to a config file."""
        config_path = self._get_config_path(config_type)

        # Backup existing config
        if config_path.exists():
            backup_path = config_path.with_suffix(".yml.bak")
            backup_path.write_text(config_path.read_text())

        # Special handling for retrosynthesis config - unflatten to nested structure
        if config_type == "retrosynthesis":
            data = self._unflatten_retrosynthesis_config(data)

        with open(config_path, "w") as f:
            yaml.dump(
                data, f, default_flow_style=False, sort_keys=False, allow_unicode=True
            )

        return True

    def validate_config(self, config_type: str, data: dict[str, Any]) -> dict[str, Any]:
        """Validate config data using centralized ConfigValidator."""
        return ConfigValidator.validate(config_type, data)
