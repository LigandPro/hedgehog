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
            result["warnings"].append(
                f"No specific validator for config type: {config_type}"
            )

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
            if not isinstance(n_jobs, int) or (n_jobs != -1 and n_jobs < 1):
                result["errors"].append("n_jobs must be -1 or a positive integer")

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
    def _validate_mol_prep(data: dict[str, Any], result: dict[str, Any]) -> None:
        """Validate Mol Prep configuration."""
        n_jobs = data.get("n_jobs")
        if n_jobs is not None and (
            not isinstance(n_jobs, int) or (n_jobs != -1 and n_jobs < 1)
        ):
            result["errors"].append("n_jobs must be -1 or a positive integer")

    @staticmethod
    def _validate_filters(data: dict[str, Any], result: dict[str, Any]) -> None:
        """Validate structural filters configuration."""
        # No specific validation rules defined yet
        pass

    @staticmethod
    def _validate_synthesis(data: dict[str, Any], result: dict[str, Any]) -> None:
        """Validate synthesis configuration."""
        if "sa_score_min" in data:
            sa_min = data["sa_score_min"]
            if not isinstance(sa_min, (int, float)) or sa_min < 1 or sa_min > 10:
                result["errors"].append("sa_score_min must be between 1 and 10")

        if "sa_score_max" in data:
            sa_max = data["sa_score_max"]
            if not isinstance(sa_max, (int, float)) or sa_max < 1 or sa_max > 10:
                result["errors"].append("sa_score_max must be between 1 and 10")

        if "ra_score_min" in data:
            ra_min = data["ra_score_min"]
            if not isinstance(ra_min, (int, float)) or ra_min < 0 or ra_min > 1:
                result["errors"].append("ra_score_min must be between 0 and 1")

    @staticmethod
    def _validate_retrosynthesis(data: dict[str, Any], result: dict[str, Any]) -> None:
        """Validate retrosynthesis configuration.

        Retrosynthesis config uses complex nested AiZynthFinder format;
        detailed validation is not yet implemented.
        """
        pass

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

    @staticmethod
    def _validate_docking_filters(data: dict[str, Any], result: dict[str, Any]) -> None:
        """Validate docking filters configuration."""
        if not data.get("run", False):
            return

        aggregation = data.get("aggregation", {}) or {}
        mode = aggregation.get("mode", "all")
        if mode not in {"all", "any"}:
            result["errors"].append("aggregation.mode must be 'all' or 'any'")

        pose_quality = data.get("pose_quality", {}) or {}
        max_clashes = pose_quality.get("max_clashes")
        if max_clashes is not None and (
            not isinstance(max_clashes, int) or max_clashes < 0
        ):
            result["errors"].append(
                "pose_quality.max_clashes must be a non-negative integer"
            )

        conformer_dev = data.get("conformer_deviation", {}) or {}
        max_rmsd = conformer_dev.get("max_rmsd_to_conformer")
        if max_rmsd is not None and (
            not isinstance(max_rmsd, (int, float)) or max_rmsd < 0
        ):
            result["errors"].append(
                "conformer_deviation.max_rmsd_to_conformer must be a non-negative number"
            )

        shepherd_cfg = data.get("shepherd_score", {}) or {}
        shepherd_backend = shepherd_cfg.get("backend", "auto")
        valid_shepherd_backends = {"auto", "worker", "inprocess"}
        if shepherd_backend not in valid_shepherd_backends:
            result["errors"].append(
                "shepherd_score.backend must be one of: auto, worker, inprocess"
            )
