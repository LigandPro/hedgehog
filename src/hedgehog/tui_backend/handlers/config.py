"""Config handling for TUI backend."""

from pathlib import Path
from typing import TYPE_CHECKING, Any

import yaml

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
        """Validate config data."""
        errors = []

        if config_type == "main":
            if not data.get("generated_mols_path"):
                errors.append("generated_mols_path is required")
            if not data.get("folder_to_save"):
                errors.append("folder_to_save is required")
            if data.get("n_jobs", 0) < 1:
                errors.append("n_jobs must be at least 1")
            if data.get("sample_size", 0) < 1:
                errors.append("sample_size must be at least 1")

        elif config_type == "descriptors":
            borders = data.get("borders", {})
            for key in ["molWt_min", "molWt_max", "logP_min", "logP_max"]:
                if key in borders and not isinstance(borders[key], (int, float)):
                    errors.append(f"{key} must be a number")

        elif config_type == "synthesis":
            if data.get("sa_score_min", 0) < 0:
                errors.append("sa_score_min must be non-negative")
            if data.get("ra_score_min", 0) < 0 or data.get("ra_score_min", 0) > 1:
                errors.append("ra_score_min must be between 0 and 1")

        elif config_type == "docking":
            if not data.get("receptor_pdb"):
                errors.append("receptor_pdb is required")
            tools = data.get("tools", "")
            if tools not in ["both", "smina", "gnina"]:
                errors.append("tools must be one of: both, smina, gnina")

        return {
            "valid": len(errors) == 0,
            "errors": errors,
        }
