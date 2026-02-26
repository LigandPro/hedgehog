from pathlib import Path

from hedgehog._constants import TOOL_GNINA
from hedgehog.configs.logger import logger
from hedgehog.docking.paths import _resolve_autobox_path, _resolve_path


def _skip_keys_for_tool(tool_name: str) -> set[str]:
    """Return the set of config keys to skip for a given docking tool."""
    keys = {"bin", "center", "size"}
    if tool_name == TOOL_GNINA:
        keys.update({"env_path", "ld_library_path", "activate", "output_dir", "no_gpu"})
    return keys


def _build_config_lines(
    receptor: str,
    ligand: str,
    output: str,
    tool_config: dict,
    cfg: dict,
    skip_keys: set[str],
    autobox_path=None,
) -> list[str]:
    """Build docking config lines shared between batch and per-molecule modes."""
    lines = [
        f"receptor = {receptor}",
        f"ligand = {ligand}",
        f"out = {output}",
    ]

    center = tool_config.get("center") or cfg.get("center")
    if center and isinstance(center, (list, tuple)) and len(center) >= 3:
        lines.extend(
            [
                f"center_x = {center[0]}",
                f"center_y = {center[1]}",
                f"center_z = {center[2]}",
            ]
        )

    size = tool_config.get("size") or cfg.get("size")
    if size and isinstance(size, (list, tuple)) and len(size) >= 3:
        lines.extend(
            [f"size_x = {size[0]}", f"size_y = {size[1]}", f"size_z = {size[2]}"]
        )

    if autobox_path:
        lines.append(f"autobox_ligand = {autobox_path}")

    for key, value in tool_config.items():
        if value is None or key in skip_keys or key == "autobox_ligand":
            continue
        if isinstance(value, (list, tuple)):
            lines.append(f"{key} = [{', '.join(str(v) for v in value)}]")
        elif isinstance(value, bool):
            lines.append(f"{key} = {str(value).lower()}")
        else:
            lines.append(f"{key} = {value}")

    return lines


def _resolve_tool_autobox(tool_config: dict) -> str | None:
    """Resolve the autobox_ligand path from tool config, returning absolute path or None."""
    autobox_ligand = tool_config.get("autobox_ligand")
    if not autobox_ligand:
        return None
    project_root = Path(__file__).parent.parent.parent.parent
    autobox_path = _resolve_autobox_path(autobox_ligand, project_root)
    if autobox_path:
        return str(autobox_path)
    logger.warning(
        "Could not resolve autobox_ligand path: %s. "
        "Using as-is (may fail if relative path is incorrect).",
        autobox_ligand,
    )
    return None


def _create_docking_config_file(
    cfg, ligands_dir, receptor, ligands_path, output_sdf, config_path, tool_name
):
    """Create docking tool config file from configuration arguments."""
    config_path.parent.mkdir(parents=True, exist_ok=True)
    Path(output_sdf).parent.mkdir(parents=True, exist_ok=True)

    tool_config = dict(cfg.get(f"{tool_name}_config", {}) or {})
    tool_config["num_modes"] = 1

    receptor = _resolve_path(receptor, ligands_dir)
    ligands_path = _resolve_path(ligands_path, ligands_dir)
    output_sdf = _resolve_path(output_sdf, ligands_dir)

    autobox_path = _resolve_tool_autobox(tool_config)
    if autobox_path:
        tool_config["autobox_ligand"] = autobox_path

    skip_keys = _skip_keys_for_tool(tool_name)
    lines = _build_config_lines(
        receptor, ligands_path, output_sdf, tool_config, cfg, skip_keys, autobox_path
    )

    with open(config_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    logger.debug("Created %s config file: %s", tool_name.upper(), config_path)
    return config_path


def _create_per_molecule_configs(
    cfg,
    ligands_dir,
    receptor,
    molecule_files,
    tool_name,
    cpu_override: int | None = None,
):
    """Create per-molecule docking config files."""
    configs_dir = ligands_dir / "_workdir" / "configs"
    results_dir = ligands_dir / "_workdir" / tool_name / "results"
    configs_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    tool_config = dict(cfg.get(f"{tool_name}_config", {}) or {})
    if cpu_override is not None:
        tool_config["cpu"] = int(cpu_override)
    tool_config["num_modes"] = 1

    receptor_abs = _resolve_path(receptor, ligands_dir)
    autobox_path = _resolve_tool_autobox(tool_config)
    skip_keys = _skip_keys_for_tool(tool_name)
    config_entries = []

    for mol_file in molecule_files:
        mol_id = mol_file.stem
        config_path = configs_dir / f"{tool_name}_{mol_id}.ini"
        output_sdf = results_dir / f"{mol_id}_out.sdf"

        lines = _build_config_lines(
            receptor_abs,
            str(mol_file.resolve()),
            str(output_sdf.resolve()),
            tool_config,
            cfg,
            skip_keys,
            autobox_path,
        )

        with open(config_path, "w") as f:
            f.write("\n".join(lines) + "\n")

        config_entries.append((mol_id, config_path, output_sdf))

    logger.info(
        "Created %d per-molecule config files in %s", len(config_entries), configs_dir
    )
    return config_entries


def _create_smina_config_file(
    cfg, ligands_dir, receptor, ligands_path, config_path, output_sdf
):
    """Create SMINA config file from configuration arguments."""
    return _create_docking_config_file(
        cfg, ligands_dir, receptor, ligands_path, output_sdf, config_path, "smina"
    )


def _create_gnina_config_file(
    cfg, ligands_dir, receptor, ligands_path, output_sdf, config_path
):
    """Create GNINA config file from configuration arguments."""
    return _create_docking_config_file(
        cfg, ligands_dir, receptor, ligands_path, output_sdf, config_path, "gnina"
    )


def _parse_bool_config(value, default: bool = False) -> bool:
    """Parse a boolean from config values that may be bool, int, float, or str."""
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "on"}
    return default
