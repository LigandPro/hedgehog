"""Centralized configuration validation for TUI backend."""

from typing import Any


class ConfigValidator:
    """Centralized validation logic for configuration types."""

    @staticmethod
    def validate(config_type: str, data: dict[str, Any]) -> dict[str, Any]:
        """Validate config data and return validation result.

        Args:
            config_type: Type of configuration ('main', 'descriptors', 'filters',
                'synthesis', 'docking')
            data: Configuration dictionary to validate

        Returns:
            Dictionary with 'valid' (bool), 'errors' (list), and optionally 'warnings' (list)
        """
        result = {
            "valid": True,
            "errors": [],
            "warnings": [],
        }

        validator_method = getattr(ConfigValidator, f"_validate_{config_type}", None)
        if validator_method:
            validator_method(data, result)
        else:
            # Unknown config type - no specific validation
            pass

        result["valid"] = len(result["errors"]) == 0
        return result

    @staticmethod
    def _validate_main(data: dict[str, Any], result: dict[str, Any]) -> None:
        """Validate main configuration."""
        required_fields = ["generated_mols_path", "folder_to_save"]
        for field in required_fields:
            if field not in data or not data[field]:
                result["errors"].append(f"Missing required field: {field}")

        if "n_jobs" in data:
            n_jobs = data["n_jobs"]
            if not isinstance(n_jobs, int) or n_jobs < 1:
                result["errors"].append("n_jobs must be a positive integer")

        if "sample_size" in data:
            sample_size = data["sample_size"]
            if sample_size is not None and (
                not isinstance(sample_size, int) or sample_size < 1
            ):
                result["errors"].append(
                    "sample_size must be a positive integer or null"
                )

    @staticmethod
    def _validate_descriptors(data: dict[str, Any], result: dict[str, Any]) -> None:
        """Validate descriptors configuration."""
        borders = data.get("borders", {})
        numeric_keys = ["molWt_min", "molWt_max", "logP_min", "logP_max"]
        for key in numeric_keys:
            if key in borders and not isinstance(borders[key], (int, float)):
                result["errors"].append(f"{key} must be a number")

    @staticmethod
    def _validate_filters(data: dict[str, Any], result: dict[str, Any]) -> None:
        """Validate structural filters configuration."""
        # No specific validation rules defined yet
        pass

    @staticmethod
    def _validate_synthesis(data: dict[str, Any], result: dict[str, Any]) -> None:
        """Validate synthesis configuration."""
        if "sa_score_min" in data:
            if data["sa_score_min"] < 0:
                result["errors"].append("sa_score_min must be non-negative")

        if "ra_score_min" in data:
            ra_min = data["ra_score_min"]
            if ra_min < 0 or ra_min > 1:
                result["errors"].append("ra_score_min must be between 0 and 1")

    @staticmethod
    def _validate_docking(data: dict[str, Any], result: dict[str, Any]) -> None:
        """Validate docking configuration."""
        if data.get("run", False):
            if not data.get("receptor_pdb"):
                result["errors"].append(
                    "receptor_pdb is required when docking is enabled"
                )

        tools = data.get("tools", "both")
        valid_tools = {"both", "smina", "gnina"}
        if tools not in valid_tools:
            result["errors"].append(f"tools must be one of: {', '.join(valid_tools)}")
