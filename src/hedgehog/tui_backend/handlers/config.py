"""Config handling for TUI backend."""

import hashlib
import shutil
from pathlib import Path
from typing import TYPE_CHECKING, Any

import yaml

from ..validators import ConfigValidator

if TYPE_CHECKING:
    from ..server import JsonRpcServer

# Default source config paths relative to project root
CONFIG_SOURCE_PATHS = {
    "main": "src/hedgehog/configs/config.yml",
    "mol_prep": "src/hedgehog/configs/config_mol_prep.yml",
    "descriptors": "src/hedgehog/configs/config_descriptors.yml",
    "filters": "src/hedgehog/configs/config_structFilters.yml",
    "synthesis": "src/hedgehog/configs/config_synthesis.yml",
    "retrosynthesis": "modules/retrosynthesis/aizynthfinder/public/config.yml",
    "docking": "src/hedgehog/configs/config_docking.yml",
    "docking_filters": "src/hedgehog/configs/config_docking_filters.yml",
}

# Workspace-local filenames for editable TUI copies
CONFIG_WORKSPACE_FILES = {
    "main": "config.yml",
    "mol_prep": "config_mol_prep.yml",
    "descriptors": "config_descriptors.yml",
    "filters": "config_structFilters.yml",
    "synthesis": "config_synthesis.yml",
    "retrosynthesis": "config_retrosynthesis.yml",
    "docking": "config_docking.yml",
    "docking_filters": "config_docking_filters.yml",
}

# Mapping from master config keys to config types used by TUI
RUNTIME_CONFIG_KEY_MAP = {
    "config_mol_prep": "mol_prep",
    "config_descriptors": "descriptors",
    "config_structFilters": "filters",
    "config_synthesis": "synthesis",
    "config_docking": "docking",
    "config_docking_filters": "docking_filters",
    "config_retrosynthesis": "retrosynthesis",
}


class ConfigHandler:
    """Handler for config-related RPC methods."""

    def __init__(self, server: "JsonRpcServer"):
        self.server = server
        self._project_root = self._find_project_root()
        self._user_config_root = self._build_user_config_root()

    def _find_project_root(self) -> Path:
        """Find project root by looking for pyproject.toml."""
        current = Path.cwd()
        for parent in [current, *current.parents]:
            if (parent / "pyproject.toml").exists():
                return parent
        return current

    def _build_user_config_root(self) -> Path:
        """Return per-project user config directory for editable TUI copies."""
        project = str(self._project_root.resolve())
        digest = hashlib.sha1(project.encode("utf-8")).hexdigest()[:12]
        workspace = "".join(
            ch if ch.isalnum() or ch in {"-", "_"} else "_"
            for ch in self._project_root.name
        ).strip("_")
        workspace = workspace or "workspace"
        return Path.home() / ".hedgehog" / "tui" / "configs" / f"{workspace}-{digest}"

    def _get_source_config_path(self, config_type: str) -> Path:
        """Get source template config path from repository checkout."""
        if config_type not in CONFIG_SOURCE_PATHS:
            raise ValueError(f"Unknown config type: {config_type}")
        return self._project_root / CONFIG_SOURCE_PATHS[config_type]

    def _get_user_config_path(self, config_type: str) -> Path:
        """Get editable user config path for a given type."""
        if config_type not in CONFIG_WORKSPACE_FILES:
            raise ValueError(f"Unknown config type: {config_type}")
        return self._user_config_root / CONFIG_WORKSPACE_FILES[config_type]

    def _ensure_user_config_exists(self, config_type: str) -> Path:
        """Create editable user config from source template on first access."""
        user_config_path = self._get_user_config_path(config_type)
        if user_config_path.exists():
            return user_config_path

        source_config_path = self._get_source_config_path(config_type)
        if not source_config_path.exists():
            raise FileNotFoundError(f"Config file not found: {source_config_path}")

        user_config_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_config_path, user_config_path)
        return user_config_path

    def get_config_path(self, config_type: str) -> Path:
        """Get editable config path used by TUI runtime."""
        return self._ensure_user_config_exists(config_type)

    def get_runtime_config_overrides(self) -> dict[str, str]:
        """Get config key overrides so pipeline uses TUI-managed config copies."""
        return {
            config_key: str(self.get_config_path(config_type))
            for config_key, config_type in RUNTIME_CONFIG_KEY_MAP.items()
        }

    def load_config(self, config_type: str) -> dict[str, Any]:
        """Load a config file and return its contents."""
        config_path = self.get_config_path(config_type)

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
        config_path = self.get_config_path(config_type)

        validation = ConfigValidator.validate(config_type, data)
        if not validation["valid"]:
            raise ValueError(f"Invalid config: {'; '.join(validation['errors'])}")

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
