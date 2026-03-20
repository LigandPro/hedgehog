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
        borders = data.get("borders", {}) or {}
        if not isinstance(borders, dict):
            result["errors"].append("borders must be a mapping")
            return

        numeric_keys = ["molWt_min", "molWt_max", "logP_min", "logP_max"]
        for key in numeric_keys:
            if key in borders and not isinstance(borders[key], (int, float)):
                result["errors"].append(f"{key} must be a number")

        structural_constraints = data.get("structural_constraints")
        if structural_constraints is None:
            return

        if not isinstance(structural_constraints, dict):
            result["errors"].append("structural_constraints must be a mapping")
            return

        ConfigValidator._validate_integer_limit_map(
            structural_constraints=structural_constraints,
            map_key="type_limits",
            result=result,
        )
        ConfigValidator._validate_integer_limit_map(
            structural_constraints=structural_constraints,
            map_key="element_limits",
            result=result,
        )

        for scalar_key in [
            "max_n_or_o_atoms",
            "max_small_rings_3_4",
            "max_acyclic_chain_length",
        ]:
            if scalar_key not in structural_constraints:
                continue
            value = structural_constraints[scalar_key]
            if not ConfigValidator._is_non_negative_int(value):
                result["errors"].append(
                    f"structural_constraints.{scalar_key} must be a non-negative integer"
                )

    @staticmethod
    def _validate_integer_limit_map(
        structural_constraints: dict[str, Any],
        map_key: str,
        result: dict[str, Any],
    ) -> None:
        limit_map = structural_constraints.get(map_key)
        if limit_map is None:
            return

        if not isinstance(limit_map, dict):
            result["errors"].append(
                f"structural_constraints.{map_key} must be a mapping"
            )
            return

        for limit_key, limit_value in limit_map.items():
            if not isinstance(limit_key, str) or not limit_key.strip():
                result["errors"].append(
                    f"structural_constraints.{map_key} keys must be non-empty strings"
                )
                continue

            if not ConfigValidator._is_non_negative_int(limit_value):
                result["errors"].append(
                    f"structural_constraints.{map_key}.{limit_key} must be a non-negative integer"
                )

    @staticmethod
    def _is_non_negative_int(value: Any) -> bool:
        return isinstance(value, int) and not isinstance(value, bool) and value >= 0

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
        valid_tools = {"all", "smina", "gnina", "matcha"}
        legacy_valid_tools = {"both"}
        if isinstance(tools, str):
            tools_list = [
                chunk.strip().lower() for chunk in tools.split(",") if chunk.strip()
            ]
        elif isinstance(tools, (list, tuple)):
            tools_list = [
                str(chunk).strip().lower() for chunk in tools if str(chunk).strip()
            ]
        else:
            tools_list = [str(tools).strip().lower()]
        if not tools_list or any(
            tool not in valid_tools and tool not in legacy_valid_tools
            for tool in tools_list
        ):
            result["errors"].append(
                f"tools must contain only: {', '.join(sorted(valid_tools))}"
            )

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
