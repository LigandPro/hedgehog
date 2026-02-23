from pathlib import Path

from hedgehog.configs.logger import logger
from hedgehog.docking.binary import _resolve_autobox_path, _resolve_path


def _create_docking_config_file(
    cfg, ligands_dir, receptor, ligands_path, output_sdf, config_path, tool_name
):
    """Create docking tool config file from configuration arguments.

    Shared implementation for SMINA and GNINA config file generation.

    Args:
        cfg: Configuration dictionary
        ligands_dir: Directory for ligands
        receptor: Path to receptor file
        ligands_path: Path to ligands file
        output_sdf: Path for output SDF
        config_path: Path for config file
        tool_name: 'smina' or 'gnina'

    Returns:
        Path to created config file
    """
    config_path.parent.mkdir(parents=True, exist_ok=True)
    Path(output_sdf).parent.mkdir(parents=True, exist_ok=True)

    tool_config = dict(cfg.get(f"{tool_name}_config", {}) or {})
    tool_config["num_modes"] = 1

    receptor = _resolve_path(receptor, ligands_dir)
    ligands_path = _resolve_path(ligands_path, ligands_dir)
    output_sdf = _resolve_path(output_sdf, ligands_dir)

    lines = [
        f"receptor = {receptor}",
        f"ligand = {ligands_path}",
        f"out = {output_sdf}",
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

    autobox_ligand = tool_config.get("autobox_ligand") or cfg.get("autobox_ligand")
    if autobox_ligand:
        project_root = Path(__file__).parent.parent.parent.parent
        autobox_path = _resolve_autobox_path(autobox_ligand, project_root)
        if autobox_path:
            tool_config["autobox_ligand"] = str(autobox_path)
        else:
            logger.warning(
                "Could not resolve autobox_ligand path: %s. "
                "Using as-is (may fail if relative path is incorrect).",
                autobox_ligand,
            )

    skip_keys = {"bin", "center", "size"}
    if tool_name == "gnina":
        skip_keys.update(
            {"env_path", "ld_library_path", "activate", "output_dir", "no_gpu"}
        )

    for key, value in tool_config.items():
        if value is None or key in skip_keys:
            continue
        if isinstance(value, (list, tuple)):
            lines.append(f"{key} = [{', '.join(str(v) for v in value)}]")
        elif isinstance(value, bool):
            lines.append(f"{key} = {str(value).lower()}")
        else:
            lines.append(f"{key} = {value}")

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
    """Create per-molecule docking config files.

    Args:
        cfg: Configuration dictionary
        ligands_dir: Base directory for docking
        receptor: Path to receptor file
        molecule_files: List of paths to individual molecule SDF files
        tool_name: 'smina' or 'gnina'

    Returns:
        List of tuples (mol_id, config_path, output_path)
    """
    configs_dir = ligands_dir / "_workdir" / "configs"
    results_dir = ligands_dir / "_workdir" / tool_name / "results"
    configs_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    tool_config = dict(cfg.get(f"{tool_name}_config", {}) or {})
    if cpu_override is not None:
        tool_config["cpu"] = int(cpu_override)
    tool_config["num_modes"] = 1

    # Resolve common settings once
    receptor_abs = _resolve_path(receptor, ligands_dir)

    center = tool_config.get("center") or cfg.get("center")
    size = tool_config.get("size") or cfg.get("size")

    autobox_ligand = tool_config.get("autobox_ligand") or cfg.get("autobox_ligand")
    autobox_path = None
    if autobox_ligand:
        project_root = Path(__file__).parent.parent.parent.parent
        autobox_path = _resolve_autobox_path(autobox_ligand, project_root)

    skip_keys = {"bin", "center", "size"}
    if tool_name == "gnina":
        skip_keys.update(
            {"env_path", "ld_library_path", "activate", "output_dir", "no_gpu"}
        )

    config_entries = []

    for mol_file in molecule_files:
        mol_id = mol_file.stem
        config_path = configs_dir / f"{tool_name}_{mol_id}.ini"
        output_sdf = results_dir / f"{mol_id}_out.sdf"

        lines = [
            f"receptor = {receptor_abs}",
            f"ligand = {mol_file.resolve()}",
            f"out = {output_sdf.resolve()}",
        ]

        if center and isinstance(center, (list, tuple)) and len(center) >= 3:
            lines.extend(
                [
                    f"center_x = {center[0]}",
                    f"center_y = {center[1]}",
                    f"center_z = {center[2]}",
                ]
            )

        if size and isinstance(size, (list, tuple)) and len(size) >= 3:
            lines.extend(
                [
                    f"size_x = {size[0]}",
                    f"size_y = {size[1]}",
                    f"size_z = {size[2]}",
                ]
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
