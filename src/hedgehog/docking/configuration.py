import os
import re
from pathlib import Path

from hedgehog._constants import TOOL_GNINA, TOOL_MATCHA, TOOL_SMINA


class DockingConfigError(ValueError):
    """Raised when a docking configuration cannot produce a valid run."""


_ENV_VAR_PATTERN = re.compile(r"\$(?:\{[^}]+\}|[A-Za-z_][A-Za-z0-9_]*)")
_TOP_LEVEL_KEYS = {
    "auto_run",
    "autobox_add",
    "autobox_ligand",
    "calculate_score_thresholds_from_targets",
    "center",
    "gnina_activate",
    "gnina_bin",
    "gnina_config",
    "gnina_container",
    "gnina_env_path",
    "gnina_ld_library_path",
    "gnina_ligands",
    "gnina_output_dir",
    "gnina_parallel_jobs",
    "gnina_parallel_jobs_max",
    "gnina_parallel_jobs_scale",
    "gnina_per_process_cpu",
    "matcha_config",
    "matcha_ligands",
    "per_molecule_docking",
    "prep_njobs",
    "prepare_ligands",
    "receptor_pdb",
    "run",
    "run_in_background",
    "score_thresholds",
    "size",
    "smina_bin",
    "smina_config",
    "smina_ligands",
    "smina_parallel_jobs",
    "smina_parallel_jobs_max",
    "smina_per_process_cpu",
    "tools",
}
_MATCHA_CONFIG_KEYS = {
    "autobox_ligand",
    "backend",
    "center",
    "center_x",
    "center_y",
    "center_z",
    "checkpoint_root",
    "checkpoint_run",
    "checkpoints",
    "checkout_dir",
    "config",
    "device",
    "docking_batch_limit",
    "gnina_batch_mode",
    "keep_workdir",
    "n_confs",
    "n_samples",
    "num_workers",
    "persistent_workers",
    "prefetch_factor",
    "repo_url",
    "run_name",
    "scorer",
    "scorer_minimize",
    "scorer_path",
    "target_name",
}


_COMMON_ENGINE_ARGUMENTS = {
    "autobox_add",
    "cpu",
    "energy_range",
    "exhaustiveness",
    "min_rmsd_filter",
    "num_modes",
    "seed",
}
_ENGINE_ARGUMENTS = {
    TOOL_SMINA: _COMMON_ENGINE_ARGUMENTS
    | {"atom_terms", "custom_scoring", "flex", "local_only", "minimize", "scoring"},
    TOOL_GNINA: _COMMON_ENGINE_ARGUMENTS
    | {
        "cnn",
        "cnn_model",
        "cnn_rotation",
        "cnn_scoring",
        "device",
        "flex",
        "local_only",
        "minimize",
        "scoring",
    },
}
_ENGINE_CONTROL_ARGUMENTS = {
    "activate",
    "autobox_ligand",
    "autobox_receptor_distance_warn",
    "bin",
    "center",
    "env_path",
    "ld_library_path",
    "no_gpu",
    "output_dir",
    "size",
}


def engine_arguments(tool: str) -> set[str]:
    """Return supported native INI arguments for one docking engine."""
    try:
        return _ENGINE_ARGUMENTS[tool]
    except KeyError as exc:
        raise ValueError(f"Unsupported docking config writer tool: {tool}") from exc


def validate_engine_config(tool: str, config: dict) -> None:
    """Reject engine keys that would otherwise be silently forwarded or ignored."""
    unknown = sorted(set(config) - engine_arguments(tool) - _ENGINE_CONTROL_ARGUMENTS)
    if unknown:
        raise ValueError(f"Unsupported {tool}_config options: {unknown}")


def _resolve_input_path(value: object, config_dir: Path, field: str) -> str | None:
    """Expand and resolve an input path relative to the docking config file."""
    if value in (None, ""):
        return None

    expanded = os.path.expandvars(os.path.expanduser(str(value)))
    unresolved = _ENV_VAR_PATTERN.search(expanded)
    if unresolved:
        raise DockingConfigError(
            f"{field} references an unset environment variable: {unresolved.group(0)}"
        )

    path = Path(expanded)
    if not path.is_absolute():
        path = config_dir / path
    return str(path.resolve())


def normalize_docking_config(
    cfg: dict, config_path: str | Path, tools: list[str] | None = None
) -> dict:
    """Return a copy with selected-tool input paths resolved from the config directory."""
    normalized = dict(cfg)
    selected_tools = set(tools or (TOOL_SMINA, TOOL_GNINA, TOOL_MATCHA))
    config_dir = Path(config_path).expanduser().resolve().parent

    normalized["receptor_pdb"] = _resolve_input_path(
        normalized.get("receptor_pdb"), config_dir, "receptor_pdb"
    )
    normalized["autobox_ligand"] = _resolve_input_path(
        normalized.get("autobox_ligand"), config_dir, "autobox_ligand"
    )

    for tool in (TOOL_SMINA, TOOL_GNINA, TOOL_MATCHA):
        key = f"{tool}_config"
        raw_tool_cfg = normalized.get(key)
        if raw_tool_cfg is not None and not isinstance(raw_tool_cfg, dict):
            raise DockingConfigError(f"{key} must be a mapping")
        tool_cfg = dict(raw_tool_cfg or {})
        if tool not in selected_tools:
            normalized[key] = tool_cfg
            continue
        if "autobox_ligand" in tool_cfg:
            tool_cfg["autobox_ligand"] = _resolve_input_path(
                tool_cfg.get("autobox_ligand"),
                config_dir,
                f"{key}.autobox_ligand",
            )
        normalized[key] = tool_cfg

    matcha_cfg = normalized["matcha_config"]
    if TOOL_MATCHA not in selected_tools:
        return normalized

    if "checkpoint_root" in matcha_cfg:
        matcha_cfg["checkpoint_root"] = _resolve_input_path(
            matcha_cfg.get("checkpoint_root"),
            config_dir,
            "matcha_config.checkpoint_root",
        )

    return normalized


def _has_vector(config: dict, name: str) -> bool:
    value = config.get(name)
    return isinstance(value, (list, tuple)) and len(value) == 3


def _validate_existing_file(value: object, field: str, errors: list[str]) -> None:
    if not value:
        errors.append(f"{field} is required")
        return
    path = Path(str(value))
    if not path.is_file():
        errors.append(f"{field} does not exist or is not a file: {path}")


def validate_docking_config(cfg: dict, tools: list[str]) -> None:
    """Validate only the engines selected for this run, then fail once with context."""
    errors: list[str] = []
    unknown_top_level = sorted(set(cfg) - _TOP_LEVEL_KEYS)
    if unknown_top_level:
        errors.append(
            "unsupported top-level option(s): " + ", ".join(unknown_top_level)
        )

    if TOOL_MATCHA in tools:
        matcha_value = cfg.get("matcha_config") or {}
        if isinstance(matcha_value, dict):
            unknown_matcha = sorted(set(matcha_value) - _MATCHA_CONFIG_KEYS)
            if unknown_matcha:
                errors.append(
                    "unsupported matcha_config option(s): " + ", ".join(unknown_matcha)
                )

    _validate_existing_file(cfg.get("receptor_pdb"), "receptor_pdb", errors)

    shared_autobox = cfg.get("autobox_ligand")
    shared_center = _has_vector(cfg, "center")
    shared_size = _has_vector(cfg, "size")

    for tool in tools:
        tool_cfg = cfg.get(f"{tool}_config") or {}
        if not isinstance(tool_cfg, dict):
            errors.append(f"{tool}_config must be a mapping")
            continue
        if tool in (TOOL_SMINA, TOOL_GNINA):
            try:
                validate_engine_config(tool, tool_cfg)
            except ValueError as exc:
                errors.append(str(exc))

        autobox = tool_cfg.get("autobox_ligand") or shared_autobox
        has_center = _has_vector(tool_cfg, "center") or shared_center
        has_size = _has_vector(tool_cfg, "size") or shared_size
        explicit_xyz = all(tool_cfg.get(f"center_{axis}") is not None for axis in "xyz")

        if autobox:
            _validate_existing_file(autobox, f"{tool} autobox_ligand", errors)
        elif tool in (TOOL_SMINA, TOOL_GNINA) and not (has_center and has_size):
            errors.append(
                f"{tool} requires autobox_ligand or complete center and size vectors"
            )
        elif tool == TOOL_MATCHA and not (has_center or explicit_xyz):
            errors.append(f"{tool} requires autobox_ligand or a complete center vector")

    matcha_value = cfg.get("matcha_config") or {}
    matcha_cfg = matcha_value if isinstance(matcha_value, dict) else {}
    backend = str(matcha_cfg.get("backend") or "matcha_cli").strip().lower()
    if TOOL_MATCHA in tools and backend == "docking":
        for field in ("checkpoint_root", "checkpoint_run"):
            if not matcha_cfg.get(field):
                errors.append(f"matcha_config.{field} is required for backend=docking")
        checkpoint_root = matcha_cfg.get("checkpoint_root")
        if checkpoint_root and not Path(str(checkpoint_root)).is_dir():
            errors.append(
                "matcha_config.checkpoint_root does not exist or is not a directory: "
                f"{checkpoint_root}"
            )
        training_config = (
            Path(str(checkpoint_root))
            / str(matcha_cfg.get("checkpoint_run"))
            / "config.yaml"
            if checkpoint_root and matcha_cfg.get("checkpoint_run")
            else None
        )
        if training_config and not training_config.is_file():
            errors.append(
                "derived Matcha training config does not exist or is not a file: "
                f"{training_config}"
            )

    if errors:
        details = "\n- ".join(errors)
        raise DockingConfigError(f"Invalid docking configuration:\n- {details}")
