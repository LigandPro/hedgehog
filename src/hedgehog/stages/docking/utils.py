import json
import os
import shlex
import shutil
import subprocess
import time
import uuid
from datetime import datetime
from pathlib import Path

import pandas as pd

from hedgehog.configs.logger import load_config, logger
from hedgehog.utils.datamol_import import import_datamol_quietly
from hedgehog.utils.input_paths import find_latest_input_source

dm = import_datamol_quietly()


def _validate_optional_tool_path(tool_path, tool_label):
    """Validate optional external tool path and return usable value or None."""
    if not tool_path:
        return None

    path = Path(str(tool_path))
    if path.exists():
        if not path.is_file() or not os.access(path, os.X_OK):
            logger.warning(
                "%s is not executable: %s. Falling back to built-in behavior.",
                tool_label,
                tool_path,
            )
            return None
        return str(path)

    resolved = shutil.which(str(tool_path))
    if resolved:
        return resolved

    logger.warning(
        "%s not found: %s. Falling back to built-in behavior.",
        tool_label,
        tool_path,
    )
    return None


def _is_real_binary(path: str) -> bool:
    """Check if a file is a real compiled binary (ELF), not a script wrapper."""
    try:
        with open(path, "rb") as f:
            header = f.read(4)
        return header == b"\x7fELF"
    except Exception:
        return False


def _resolve_docking_binary(config_path: str, tool_name: str) -> str:
    """Resolve a docking binary path from config or PATH.

    Args:
        config_path: Path from config (absolute path or bare tool name).
        tool_name: Tool name for PATH lookup (e.g. 'smina', 'gnina').

    Returns:
        Resolved absolute path to the binary.

    Raises:
        FileNotFoundError: If the binary cannot be found.
    """
    if os.path.isabs(config_path) and os.path.isfile(config_path):
        return config_path

    found = shutil.which(tool_name)
    if found and _is_real_binary(found):
        return found
    elif found:
        logger.debug(
            "%s found at %s but is a script wrapper, not a real binary — skipping",
            tool_name,
            found,
        )

    if tool_name == "gnina":
        from hedgehog.setup import ensure_gnina

        try:
            return ensure_gnina()
        except RuntimeError as exc:
            raise FileNotFoundError(str(exc)) from exc

    raise FileNotFoundError(
        f"Docking binary '{tool_name}' not found. "
        f"Provide absolute path in config or ensure it's on PATH."
    )


def _find_latest_input_source(base_folder):
    """Find the most recent input source file for docking.

    Supports both new hierarchical structure and legacy flat structure.
    """
    path = find_latest_input_source(base_folder)
    if path:
        logger.debug("Using docking input: %s", path)
    return path


def _prepare_ligands_dataframe(df, output_csv):
    """Prepare ligands CSV from input DataFrame with SMILES validation."""
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    skipped_smiles = []

    for _, row in df.iterrows():
        smi = str(row["smiles"])
        try:
            mol = dm.to_mol(smi)
            if mol is None:
                skipped_smiles.append(smi)
                continue
        except Exception:
            skipped_smiles.append(smi)
            continue

        model_name = str(row["model_name"])
        mol_idx = str(row["mol_idx"])

        rows.append(
            {
                "smiles": smi,
                "name": mol_idx,
                "model_name": model_name,
                "mol_idx": mol_idx,
            }
        )

    output_df = pd.DataFrame(rows, columns=["smiles", "name", "model_name", "mol_idx"])
    output_df.to_csv(output_csv, index=False)

    if skipped_smiles:
        skip_path = output_csv.parent / "skipped_smiles.txt"
        with open(skip_path, "w") as f:
            for smi in skipped_smiles:
                f.write(f"{smi}\n")
        logger.warning(
            "Some SMILES could not be parsed for docking: %d/%d. See %s",
            len(skipped_smiles),
            len(df),
            skip_path,
        )
    return {
        "csv_path": str(output_csv),
        "total": len(df),
        "written": len(rows),
        "skipped": len(skipped_smiles),
    }


def _resolve_path(path, base_dir):
    """Resolve path to absolute, using base_dir if relative."""
    path_obj = Path(path)
    if path_obj.is_absolute():
        return str(path_obj.resolve())
    return str((base_dir / path).resolve())


def _resolve_autobox_path(autobox_ligand, project_root):
    """Resolve autobox_ligand path to absolute."""
    path = Path(autobox_ligand)
    if path.is_absolute():
        return path if path.exists() else None

    candidate = (project_root / autobox_ligand).resolve()
    if candidate.exists():
        return candidate

    if "data/" in autobox_ligand:
        data_path = autobox_ligand[autobox_ligand.find("data/") :]
        candidate = (project_root / data_path).resolve()
        if candidate.exists():
            return candidate

    return None


def _pdb_atom_coordinates(pdb_path: Path) -> list[tuple[float, float, float]]:
    """Extract atom coordinates from a PDB file (ATOM/HETATM records only)."""
    coords: list[tuple[float, float, float]] = []
    try:
        for line in pdb_path.read_text(encoding="utf-8", errors="ignore").splitlines():
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except Exception:
                continue
            coords.append((x, y, z))
    except OSError:
        return []
    return coords


def _sdf_center(sdf_path: Path) -> tuple[float, float, float] | None:
    """Compute center of the first conformer in an SDF file (mean of atom positions)."""
    try:
        from rdkit import Chem
    except Exception:
        return None

    suppl = Chem.SDMolSupplier(str(sdf_path))
    mol = next((m for m in suppl if m is not None), None)
    if mol is None or mol.GetNumConformers() == 0:
        return None
    conf = mol.GetConformer()
    n = mol.GetNumAtoms()
    sx = sy = sz = 0.0
    for i in range(n):
        p = conf.GetAtomPosition(i)
        sx += float(p.x)
        sy += float(p.y)
        sz += float(p.z)
    return (sx / n, sy / n, sz / n)


def _min_distance_to_point(
    coords: list[tuple[float, float, float]], point: tuple[float, float, float]
) -> float | None:
    """Compute the minimum Euclidean distance from coords to a point."""
    if not coords:
        return None
    px, py, pz = point
    best = None
    for x, y, z in coords:
        dx = x - px
        dy = y - py
        dz = z - pz
        d2 = dx * dx + dy * dy + dz * dz
        if best is None or d2 < best:
            best = d2
    return (best or 0.0) ** 0.5


def _warn_if_autobox_far_from_receptor(cfg: dict, tool_name: str) -> None:
    """Warn if autobox reference ligand seems far away from the receptor coordinates.

    This typically indicates a mismatched coordinate frame between receptor PDB
    and autobox ligand (e.g., different reference, different prepared structure),
    resulting in docking running in the wrong location.
    """
    try:
        receptor_raw = cfg.get("receptor_pdb")
        if not receptor_raw:
            return

        receptor_path = _resolve_receptor_path(receptor_raw)
        if receptor_path is None:
            return

        tool_cfg = cfg.get(f"{tool_name}_config", {}) or {}
        autobox_ligand = tool_cfg.get("autobox_ligand") or cfg.get("autobox_ligand")
        if not autobox_ligand:
            return

        project_root = Path(__file__).parent.parent.parent.parent.parent
        autobox_path = _resolve_autobox_path(str(autobox_ligand), project_root)
        if autobox_path is None:
            return

        center = _sdf_center(Path(autobox_path))
        if center is None:
            return

        protein_coords = _pdb_atom_coordinates(Path(receptor_path))
        min_dist = _min_distance_to_point(protein_coords, center)
        if min_dist is None:
            return

        warn_threshold = float(tool_cfg.get("autobox_receptor_distance_warn", 10.0))
        if min_dist > warn_threshold:
            logger.warning(
                "%s: Autobox reference ligand appears far from receptor (min dist %.2f A). "
                "This may indicate a wrong/mismatched autobox_ligand or receptor coordinate frame.",
                tool_name.upper(),
                min_dist,
            )
    except Exception:
        # Never fail docking due to a warning-only heuristic.
        return


def _count_box_warnings(log_path: Path) -> dict[str, int]:
    """Count GNINA/Vina 'outside box' warnings in a docking log."""
    import re

    patterns = ("outside box", "not within box")
    counts = {"lines": 0, "unique_molecules": 0}
    try:
        text = log_path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return counts

    mol_ids: set[str] = set()
    for line in text.splitlines():
        low = line.lower()
        if not any(p in low for p in patterns):
            continue
        counts["lines"] += 1
        m = re.match(r"^([^|]+?)\s*\|", line)
        if m:
            mol_id = m.group(1)
            if mol_id:
                mol_ids.add(mol_id.strip())
    counts["unique_molecules"] = len(mol_ids)
    return counts


def _gnina_zero_affinity_count(output_sdf: Path) -> tuple[int, int]:
    """Count how many poses have minimizedAffinity == 0.0 in GNINA SDF output."""
    try:
        from rdkit import Chem
    except Exception:
        return (0, 0)

    total = 0
    zero = 0
    suppl = Chem.SDMolSupplier(str(output_sdf))
    for mol in suppl:
        if mol is None:
            continue
        total += 1
        if mol.HasProp("minimizedAffinity"):
            try:
                if float(mol.GetProp("minimizedAffinity")) == 0.0:
                    zero += 1
            except Exception:
                continue
    return (zero, total)


def _emit_post_docking_warnings(
    tool_name: str,
    log_path: Path | None,
    output_sdf: Path | None = None,
) -> None:
    """Emit warnings about common docking failure modes based on tool logs/outputs."""
    try:
        if log_path and log_path.exists():
            counts = _count_box_warnings(log_path)
            if counts["lines"] > 0:
                logger.warning(
                    "%s log contains %d box warning lines across %d molecules. "
                    "This often indicates ligands started outside the search box or a misconfigured box.",
                    tool_name.upper(),
                    counts["lines"],
                    counts["unique_molecules"],
                )

        if tool_name.lower() == "gnina" and output_sdf and output_sdf.exists():
            zero, total = _gnina_zero_affinity_count(output_sdf)
            if total > 0 and zero > 0:
                logger.warning(
                    "GNINA output has minimizedAffinity == 0.0 for %d/%d poses. "
                    "If many, docking may have effectively failed (e.g., wrong box / no contacts).",
                    zero,
                    total,
                )
    except Exception:
        # Warning-only logic must never fail docking.
        return


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

    autobox_ligand = tool_config.get("autobox_ligand")
    if autobox_ligand:
        project_root = Path(__file__).parent.parent.parent.parent.parent
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

    autobox_ligand = tool_config.get("autobox_ligand")
    autobox_path = None
    if autobox_ligand:
        project_root = Path(__file__).parent.parent.parent.parent.parent
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


def _extract_pdb_output_from_cmd(cmd, ligands_dir):
    """Extract PDB output file paths from a protein preparation command.

    Args:
        cmd: Protein preparation command (list of args or string)
        ligands_dir: Directory for docking files

    Returns:
        Tuple of (output_filename, absolute_path) or (None, None) if not found
    """
    if not cmd:
        return None, None
    if isinstance(cmd, (list, tuple)):
        parts = [str(p) for p in cmd]
    else:
        # Best effort: handle quoted paths in legacy string commands.
        parts = shlex.split(str(cmd))
    for part in reversed(parts):
        if ".pdb" in part:
            if part.startswith("/"):
                return part.split("/")[-1], part
            return part, str(ligands_dir / part)
    return None, None


def _write_protein_prep_bash(f, protein_prep_cmd, ligands_dir, tool_name):
    """Write protein preparation bash code to script file.

    Shared implementation for SMINA and GNINA script generation.

    Args:
        f: File handle to write to
        protein_prep_cmd: Protein preparation command (list of args or string)
        ligands_dir: Directory for docking files
        tool_name: 'smina' or 'gnina' for log messages
    """
    f.write(f"cd {ligands_dir}\n")

    if isinstance(protein_prep_cmd, (list, tuple)):
        protein_prep_cmd_str = " ".join(shlex.quote(str(p)) for p in protein_prep_cmd)
    else:
        protein_prep_cmd_str = str(protein_prep_cmd)

    protein_output_file, protein_output_abs_path = _extract_pdb_output_from_cmd(
        protein_prep_cmd, ligands_dir
    )

    if protein_output_abs_path:
        f.write(f'mkdir -p "$(dirname "{protein_output_abs_path}")"\n')
        f.write(f'touch "{protein_output_abs_path}" 2>/dev/null || true\n')

    tool_label = f" for {tool_name.upper()}" if tool_name == "gnina" else ""
    f.write(f'echo "Running protein preparation{tool_label}..."\n')
    f.write(f"{protein_prep_cmd_str} || PREP_EXIT_CODE=$?\n")
    f.write('if [ ! -z "$PREP_EXIT_CODE" ]; then\n')
    f.write(
        '  echo "WARNING: Protein preparation command exited with code $PREP_EXIT_CODE"\n'
    )
    f.write("fi\n")

    if protein_output_abs_path:
        f.write("# If file exists in current dir but not at absolute path, move it\n")
        f.write(
            f'if [ -f "{protein_output_file}" ] && [ ! -f "{protein_output_abs_path}" ]; then\n'
        )
        f.write(f'  mv "{protein_output_file}" "{protein_output_abs_path}"\n')
        f.write("fi\n")
        f.write("# Check if prepared file exists and is valid (not empty)\n")
        f.write(
            f'if [ ! -f "{protein_output_abs_path}" ] && [ ! -f "{protein_output_file}" ]; then\n'
        )
        f.write(
            f'  echo "ERROR: Protein preparation failed - output file not found at {protein_output_abs_path}"\n'
        )
        f.write('  echo "Current directory: $(pwd)"\n')
        f.write('  echo "Listing files in current directory:"\n')
        f.write("  ls -la\n")
        f.write("  exit 1\n")
        f.write(
            f'elif [ -f "{protein_output_abs_path}" ] && [ ! -s "{protein_output_abs_path}" ]; then\n'
        )
        f.write(
            f'  echo "ERROR: Protein preparation failed - output file is empty: {protein_output_abs_path}"\n'
        )
        f.write("  exit 1\n")
        f.write(
            f'elif [ -f "{protein_output_file}" ] && [ ! -s "{protein_output_file}" ]; then\n'
        )
        f.write(
            f'  echo "ERROR: Protein preparation failed - output file is empty: {protein_output_file}"\n'
        )
        f.write("  exit 1\n")
        f.write("else\n")
        success_msg = "Protein preparation completed successfully"
        if tool_name == "gnina":
            success_msg += f": {protein_output_abs_path}"
        f.write(f'  echo "{success_msg}"\n')
        if tool_name == "gnina":
            f.write(
                f'  ls -lh "{protein_output_abs_path}" 2>/dev/null || ls -lh "{protein_output_file}" 2>/dev/null || true\n'
            )
        f.write("fi\n")


def _write_ligand_prep_bash(f, preparation_cmd, prepared_output_relative):
    """Write ligand preparation bash code to script file."""
    if isinstance(preparation_cmd, (list, tuple)):
        preparation_cmd_str = " ".join(
            shlex.quote(str(p)) for p in preparation_cmd if p is not None
        )
    else:
        preparation_cmd_str = str(preparation_cmd)
    f.write(f'mkdir -p "$(dirname "{prepared_output_relative}")"\n')
    f.write('echo "Running ligand preparation..."\n')
    f.write(f"{preparation_cmd_str}\n")
    f.write("PREP_EXIT_CODE=$?\n")
    f.write("if [ $PREP_EXIT_CODE -ne 0 ]; then\n")
    f.write(
        '  echo "WARNING: Ligand preparation command exited with code $PREP_EXIT_CODE"\n'
    )
    f.write("fi\n")
    _write_file_wait_check(
        f,
        prepared_output_relative,
        f"ERROR: Preparation tool failed - output file {prepared_output_relative} not found after waiting",
        "ligand preparation",
    )
    f.write("rm -f ligands_raw.smi || true\n")


def _write_receptor_check_bash(f, receptor):
    """Write receptor file check bash code for GNINA script."""
    f.write('echo "Checking receptor file before GNINA docking..."\n')
    f.write(
        "# Final check - wait a moment and verify file still exists and is readable\n"
    )
    f.write("sleep 1\n")
    f.write(f'if [ ! -f "{receptor}" ]; then\n')
    f.write(f'  echo "ERROR: Receptor file not found: {receptor}"\n')
    f.write('  echo "Current directory: $(pwd)"\n')
    f.write('  echo "Listing files in current directory:"\n')
    f.write("  ls -la\n")
    f.write(f'  if [ -d "$(dirname "{receptor}")" ]; then\n')
    f.write('    echo "Listing files in receptor directory:"\n')
    f.write(f'    ls -la "$(dirname "{receptor}")"\n')
    f.write("  fi\n")
    f.write("  exit 1\n")
    f.write(f'elif [ ! -s "{receptor}" ]; then\n')
    f.write(f'  echo "ERROR: Receptor file is empty: {receptor}"\n')
    f.write("  exit 1\n")
    f.write(f'elif [ ! -r "{receptor}" ]; then\n')
    f.write(f'  echo "ERROR: Receptor file is not readable: {receptor}"\n')
    f.write("  exit 1\n")
    f.write("else\n")
    f.write(f'  echo "Receptor file verified: {receptor}"\n')
    f.write(f'  ls -lh "{receptor}"\n')
    f.write("fi\n")


def _write_docking_command_bash(f, ligands_dir, docking_bin, config_file, tool_name):
    """Write docking command bash code to script file."""
    f.write(f"cd {ligands_dir}\n")
    try:
        config_rel = config_file.relative_to(ligands_dir)
    except ValueError:
        config_rel = config_file.name
    f.write(f'echo "Starting {tool_name.upper()} docking with config: {config_rel}"\n')
    f.write(f"{docking_bin} --config {config_rel}\n")
    f.write("if [ $? -eq 0 ]; then\n")
    f.write(f'  echo "{tool_name.upper()} docking completed successfully"\n')
    f.write("else\n")
    f.write(f'  echo "{tool_name.upper()} docking failed with exit code $?"\n')
    f.write("  exit 1\n")
    f.write("fi\n")


def _build_gnina_command_template(cfg: dict, gnina_bin: str, ligands_dir: Path) -> str:
    """Build a GNINA command template with config placeholder.

    The returned command must contain the ``__GNINA_CONFIG__`` placeholder that
    is replaced at script generation time.
    """

    def _enabled(value) -> bool:
        if isinstance(value, bool):
            return value
        if isinstance(value, (int, float)):
            return bool(value)
        if isinstance(value, str):
            return value.strip().lower() in {"1", "true", "yes", "on"}
        return bool(value)

    placeholder = "__GNINA_CONFIG__"
    gnina_cfg = cfg.get("gnina_config", {}) or {}
    no_gpu_enabled = _enabled(gnina_cfg.get("no_gpu"))
    no_gpu_flag = " --no_gpu" if no_gpu_enabled else ""
    host_cmd = f"{shlex.quote(str(gnina_bin))} --config {placeholder}{no_gpu_flag}"

    container_cfg = cfg.get("gnina_container", {}) or {}
    if not container_cfg.get("enabled", False):
        return host_cmd

    engine = str(container_cfg.get("engine", "docker")).strip().lower()
    if engine != "docker":
        logger.warning(
            "GNINA container engine '%s' is not supported, using host binary",
            engine,
        )
        return host_cmd

    image = container_cfg.get("image")
    if not image:
        logger.warning(
            "GNINA container mode enabled, but no image is configured. Using host binary."
        )
        return host_cmd

    container_bin = str(container_cfg.get("bin", "gnina")).strip() or "gnina"
    gpu_request = container_cfg.get("gpus", "all")
    mounts = container_cfg.get("mounts")
    if not mounts:
        mounts = ["/mnt:/mnt", "/home:/home", "/tmp:/tmp"]

    mount_entries = [str(entry).strip() for entry in mounts if str(entry).strip()]
    bind_ligands_dir = f"{ligands_dir}:{ligands_dir}"
    if bind_ligands_dir not in mount_entries:
        mount_entries.append(bind_ligands_dir)

    cmd_parts = [
        "docker run --rm",
        '--user "$(id -u):$(id -g)"',
        f"-w {shlex.quote(str(ligands_dir))}",
    ]
    if gpu_request not in (None, "", False):
        cmd_parts.append(f"--gpus {shlex.quote(str(gpu_request))}")
    if container_cfg.get("propagate_cuda_visible_devices", True):
        cmd_parts.append('-e CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-}"')
    for mount in mount_entries:
        cmd_parts.append(f"-v {shlex.quote(mount)}")
    cmd_parts.append(shlex.quote(str(image)))
    cmd_parts.append(shlex.quote(container_bin))
    cmd_parts.append(f"--config {placeholder}")
    if no_gpu_enabled:
        cmd_parts.append("--no_gpu")
    return " ".join(cmd_parts)


def _create_smina_script(
    ligands_dir,
    smina_bin,
    config_file,
    protein_prep_cmd,
    preparation_cmd=None,
    prepared_output_relative=None,
):
    """Create SMINA run script (legacy single-batch mode)."""
    script_path = ligands_dir / "_workdir" / "run_smina.sh"
    script_path.parent.mkdir(parents=True, exist_ok=True)
    with open(script_path, "w") as f:
        f.write("#!/usr/bin/env bash\n")
        f.write("set -eo pipefail\n")

        if protein_prep_cmd:
            _write_protein_prep_bash(f, protein_prep_cmd, ligands_dir, "smina")

        if preparation_cmd and prepared_output_relative:
            _write_ligand_prep_bash(f, preparation_cmd, prepared_output_relative)

        _write_docking_command_bash(f, ligands_dir, smina_bin, config_file, "smina")

    os.chmod(script_path, 0o755)
    return script_path


def _create_smina_per_molecule_script(
    ligands_dir,
    smina_bin,
    protein_prep_cmd,
):
    """Create SMINA run script with per-molecule processing and error handling.

    This script processes each molecule individually, allowing:
    - Individual molecule failures without stopping the entire run
    - Detailed per-molecule logging
    - Easy retry of failed molecules
    """
    script_path = ligands_dir / "_workdir" / "run_smina.sh"
    script_path.parent.mkdir(parents=True, exist_ok=True)
    configs_dir = ligands_dir / "_workdir" / "configs"
    logs_dir = ligands_dir / "_workdir" / "smina" / "logs"

    with open(script_path, "w") as f:
        f.write("#!/usr/bin/env bash\n")
        f.write("set -o pipefail\n\n")

        if protein_prep_cmd:
            _write_protein_prep_bash(f, protein_prep_cmd, ligands_dir, "smina")

        f.write(f"cd {ligands_dir}\n")
        f.write(f'mkdir -p "{logs_dir}"\n\n')

        f.write("# Per-molecule docking with error handling\n")
        f.write("FAILED=()\n")
        f.write("SUCCESS=0\n")
        f.write("TOTAL=0\n\n")

        f.write(f'for config in "{configs_dir}"/smina_*.ini; do\n')
        f.write('    [ -e "$config" ] || continue\n')
        f.write("    mol_id=$(basename \"$config\" .ini | sed 's/smina_//')\n")
        f.write("    TOTAL=$((TOTAL + 1))\n")
        f.write('    echo "[$TOTAL] Processing $mol_id..."\n\n')

        f.write(
            f'    if {smina_bin} --config "$config" 2>> "{logs_dir}/${{mol_id}}.log"; then\n'
        )
        f.write('        echo "  $mol_id: SUCCESS"\n')
        f.write("        SUCCESS=$((SUCCESS + 1))\n")
        f.write("    else\n")
        f.write("        EXIT_CODE=$?\n")
        f.write('        echo "  $mol_id: FAILED (exit $EXIT_CODE)"\n')
        f.write('        FAILED+=("$mol_id")\n')
        f.write("        if [ $EXIT_CODE -eq 143 ]; then\n")
        f.write(
            '            echo "    -> SIGTERM detected (timeout/memory limit/killed)"\n'
        )
        f.write("        fi\n")
        f.write("    fi\n")
        f.write("done\n\n")

        f.write("# Report summary\n")
        f.write('echo ""\n')
        f.write('echo "========================================"\n')
        f.write('echo "SMINA Docking Summary"\n')
        f.write('echo "========================================"\n')
        f.write('echo "Total molecules: $TOTAL"\n')
        f.write('echo "Successful: $SUCCESS"\n')
        f.write('echo "Failed: ${#FAILED[@]}"\n\n')

        f.write("if [ ${#FAILED[@]} -gt 0 ]; then\n")
        f.write('    echo ""\n')
        f.write('    echo "Failed molecules:"\n')
        f.write('    printf "  %s\\n" "${FAILED[@]}"\n')
        f.write(
            f'    printf "%s\\n" "${{FAILED[@]}}" > "{ligands_dir}/smina/failed_molecules.txt"\n'
        )
        f.write("fi\n\n")

        f.write("# Exit with success if at least one molecule succeeded\n")
        f.write("if [ $SUCCESS -gt 0 ]; then\n")
        f.write('    echo ""\n')
        f.write(
            '    echo "Docking completed with $SUCCESS/$TOTAL molecules successful"\n'
        )
        f.write("    exit 0\n")
        f.write("else\n")
        f.write('    echo ""\n')
        f.write('    echo "ERROR: All molecules failed docking"\n')
        f.write("    exit 1\n")
        f.write("fi\n")

    os.chmod(script_path, 0o755)
    return script_path


def _prepare_protein_for_docking(receptor_pdb, ligands_dir, protein_preparation_tool):
    """Prepare protein file for docking using external tool.

    Args:
        receptor_pdb: Path to the original receptor PDB file (must exist)
        ligands_dir: Directory where docking files are prepared
        protein_preparation_tool: Path to protein preparation tool

    Returns:
        Tuple of (prepared_receptor_path, preparation_cmd) or (original_path, None)
        Note: prepared_receptor_path is where the file WILL be after preparation, not where it currently is
    """
    protein_preparation_tool = _validate_optional_tool_path(
        protein_preparation_tool, "Protein preparation tool"
    )
    if not protein_preparation_tool:
        return receptor_pdb, None

    receptor_path = Path(receptor_pdb)
    if not receptor_path.exists():
        logger.warning(
            "Original receptor file not found: %s (resolved to: %s), skipping protein preprocessing",
            receptor_pdb,
            receptor_path,
        )
        return receptor_pdb, None

    prepared_output_path = ligands_dir / "_workdir" / "protein_prepared.pdb"
    prepared_output_path.parent.mkdir(parents=True, exist_ok=True)

    receptor_absolute = str(receptor_path.resolve())
    output_absolute = str(prepared_output_path.resolve())
    cmd_args = [protein_preparation_tool, receptor_absolute, output_absolute, "-WAIT"]

    logger.info(
        "Protein preprocessing will be performed in %s: %s",
        ligands_dir,
        " ".join(shlex.quote(str(p)) for p in cmd_args),
    )
    return str(prepared_output_path.resolve()), cmd_args


def _prepare_ligands_for_docking(
    ligands_csv, ligands_dir, ligand_preparation_tool, cfg, tool_name="docking"
):
    """Prepare ligands file for docking (shared for SMINA and GNINA).

    Args:
        ligands_csv: Path to input CSV file with ligands
        ligands_dir: Directory where docking files are prepared
        ligand_preparation_tool: Path to ligand preparation tool (optional)
        cfg: Configuration dict
        tool_name: Name of tool ('smina' or 'gnina') for output file naming

    Returns:
        Tuple of (ligands_path, preparation_cmd) where:
        - ligands_path: Absolute path to prepared SDF file
        - preparation_cmd: Command to run ligand preparation (None if already SDF or using RDKit)
    """
    ligands_arg = cfg.get(f"{tool_name}_ligands")
    if ligands_arg:
        ligands_val = str(ligands_arg)
    else:
        ligands_val = str(ligands_csv.relative_to(ligands_dir))

    sdf_extensions = (".sdf", ".sdf.gz", ".osd", ".mol2")
    needs_conversion = not ligands_val.lower().endswith(sdf_extensions)

    if not needs_conversion:
        ligands_path = Path(ligands_val)
        if not ligands_path.is_absolute():
            ligands_path = (ligands_dir / ligands_val).resolve()
        return str(ligands_path), None

    ligand_preparation_tool = _validate_optional_tool_path(
        ligand_preparation_tool, "Ligand preparation tool"
    )

    if ligand_preparation_tool:
        ligands_val_lower = ligands_val.lower()
        if ligands_val_lower.endswith(".csv"):
            input_format = "-icsv"
        elif ligands_val_lower.endswith((".smi", ".ismi", ".cmi")):
            input_format = "-ismi"
        else:
            input_format = "-icsv"

        prepared_output_path = (
            ligands_dir / "_workdir" / f"prepared_for_{tool_name}.sdf"
        )
        prepared_output_path.parent.mkdir(parents=True, exist_ok=True)
        prepared_output_abs = str(prepared_output_path.resolve())

        ligands_abs = str((ligands_dir / ligands_val).resolve())
        cmd_args = [
            ligand_preparation_tool,
            input_format,
            ligands_abs,
            "-osd",
            prepared_output_abs,
            "-WAIT",
        ]

        prep_njobs_cfg = cfg.get("prep_njobs", "auto")
        if str(prep_njobs_cfg).lower() == "auto":
            auto_njobs = (
                os.environ.get("SLURM_CPUS_PER_TASK")
                or os.environ.get("MOLSCORE_NJOBS")
                or 1
            )
            prep_njobs = _parse_positive_int(auto_njobs, 1)
        else:
            prep_njobs = _parse_positive_int(prep_njobs_cfg, 1)

        if prep_njobs > 1:
            cmd_args.append("-LOCAL")
            cmd_args.extend(["-HOST", f"localhost:{prep_njobs}"])
            cmd_args.extend(["-NJOBS", str(prep_njobs)])
        return str(prepared_output_path.resolve()), cmd_args
    else:
        ligands_path, _ = _convert_with_rdkit(ligands_csv, ligands_dir)
        return ligands_path, None


def _split_sdf_to_molecules(sdf_path: Path, molecules_dir: Path) -> list[Path]:
    """Split a multi-molecule SDF into individual molecule SDF files.

    Args:
        sdf_path: Path to input SDF file with multiple molecules
        molecules_dir: Directory to write individual molecule SDFs

    Returns:
        List of paths to individual molecule SDF files
    """
    try:
        from rdkit import Chem
    except ImportError as err:
        raise RuntimeError("RDKit not available for SDF splitting") from err

    molecules_dir.mkdir(parents=True, exist_ok=True)
    molecule_files = []

    suppl = Chem.SDMolSupplier(str(sdf_path))
    for mol in suppl:
        if mol is None:
            continue

        mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else None
        if not mol_name:
            mol_name = f"mol_{len(molecule_files):05d}"

        # Sanitize molecule name for filesystem and force uniqueness with index prefix
        base_name = "".join(c if c.isalnum() or c in "-_" else "_" for c in mol_name)
        if not base_name:
            base_name = "mol"
        base_name = base_name[:180]
        safe_name = f"{len(molecule_files):06d}_{base_name}"
        mol_path = molecules_dir / f"{safe_name}.sdf"

        writer = Chem.SDWriter(str(mol_path))
        writer.write(mol)
        writer.close()

        molecule_files.append(mol_path)

    logger.info(
        "Split SDF into %d individual molecule files in %s",
        len(molecule_files),
        molecules_dir,
    )
    return molecule_files


def _parse_positive_int(value, default: int) -> int:
    """Parse a positive integer with fallback."""
    try:
        parsed = int(value)
        if parsed > 0:
            return parsed
    except Exception:
        pass
    return default


def _resolve_gnina_parallel_jobs(cfg: dict, cpu_per_process: int) -> int:
    """Resolve per-molecule parallel process count for GNINA."""
    explicit = cfg.get("gnina_parallel_jobs")
    if explicit is not None:
        return _parse_positive_int(explicit, 1)

    cpus = os.environ.get("SLURM_CPUS_PER_TASK") or os.cpu_count() or cpu_per_process
    total_cpus = _parse_positive_int(cpus, cpu_per_process)

    scale_raw = cfg.get("gnina_parallel_jobs_scale", 1.0)
    try:
        scale = float(scale_raw)
    except Exception:
        scale = 1.0
    if scale <= 0:
        scale = 1.0

    jobs = int((total_cpus * scale) // max(1, cpu_per_process))
    jobs = max(1, jobs)

    jobs_max = cfg.get("gnina_parallel_jobs_max")
    if jobs_max is not None:
        jobs = min(jobs, _parse_positive_int(jobs_max, jobs))
    return jobs


def _materialize_prepared_ligands(
    prep_cmd,
    ligands_path: Path,
    ligands_dir: Path,
    tool_name: str,
) -> bool:
    """Run ligand preparation immediately to enable per-molecule docking."""
    if ligands_path.exists():
        return True
    if not prep_cmd:
        return ligands_path.exists()

    cmd = (
        [str(part) for part in prep_cmd]
        if isinstance(prep_cmd, (list, tuple))
        else shlex.split(str(prep_cmd))
    )
    logger.info(
        "%s per-molecule mode: preparing ligands before split",
        tool_name.upper(),
    )
    try:
        subprocess.run(
            cmd,
            check=True,
            cwd=str(ligands_dir),
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except subprocess.CalledProcessError as exc:
        logger.warning(
            "%s ligand preparation failed before per-molecule split (exit=%s), falling back to batch mode",
            tool_name.upper(),
            exc.returncode,
        )
        return False

    max_wait_sec = 300
    poll_sec = 2
    waited = 0
    while not ligands_path.exists() and waited < max_wait_sec:
        time.sleep(poll_sec)
        waited += poll_sec

    if not ligands_path.exists():
        logger.warning(
            "%s ligand preparation produced no output (%s), falling back to batch mode",
            tool_name.upper(),
            ligands_path,
        )
        return False
    return True


def _aggregate_docking_results(results_dir: Path, output_sdf: Path) -> int:
    """Aggregate per-molecule docking results into a single SDF file.

    Args:
        results_dir: Directory containing per-molecule result SDFs (*_out.sdf)
        output_sdf: Path to write aggregated output SDF

    Returns:
        Number of molecules successfully aggregated
    """
    try:
        from rdkit import Chem
    except ImportError as err:
        raise RuntimeError("RDKit not available for result aggregation") from err

    output_sdf.parent.mkdir(parents=True, exist_ok=True)

    result_files = sorted(results_dir.glob("*_out.sdf"))
    if not result_files:
        logger.warning("No result files found in %s", results_dir)
        return 0

    def _extract_pose_affinity(mol) -> float | None:
        for prop_name in ("minimizedAffinity", "affinity", "score"):
            if mol.HasProp(prop_name):
                try:
                    return float(mol.GetProp(prop_name))
                except Exception:
                    continue
        return None

    def _pick_best_pose(result_file: Path):
        best_mol = None
        best_affinity = float("inf")
        for mol in Chem.SDMolSupplier(str(result_file)):
            if mol is None:
                continue
            affinity = _extract_pose_affinity(mol)
            affinity_sort = affinity if affinity is not None else float("inf")
            if best_mol is None or affinity_sort < best_affinity:
                best_mol = mol
                best_affinity = affinity_sort
        return best_mol

    writer = Chem.SDWriter(str(output_sdf))
    count = 0

    for result_file in result_files:
        try:
            best_pose = _pick_best_pose(result_file)
            if best_pose is not None:
                writer.write(best_pose)
                count += 1
        except Exception as e:
            logger.warning("Failed to read result file %s: %s", result_file, e)
            continue

    writer.close()
    logger.info(
        "Aggregated %d best poses (single pose per molecule) from %d result files into %s",
        count,
        len(result_files),
        output_sdf,
    )
    return count


def _convert_with_rdkit(ligands_csv, ligands_dir):
    """Convert SMILES to SDF using RDKit as fallback.
    Assumes CSV contains smiles and name columns."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as err:
        raise RuntimeError("RDKit not available for ligand conversion") from err

    sdf_path = ligands_dir / "_workdir" / "ligands_prepared.sdf"
    sdf_path.parent.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(ligands_csv)
    smiles_series = df["smiles"]
    name_series = df["name"]

    writer = Chem.SDWriter(str(sdf_path))
    written_count = 0

    for smi, name in zip(
        smiles_series.astype(str), name_series.astype(str), strict=False
    ):
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            mol = Chem.AddHs(mol)
            try:
                AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                AllChem.UFFOptimizeMolecule(mol)
            except Exception:
                pass
            mol.SetProp("_Name", name)
            writer.write(mol)
            written_count += 1
        except Exception:
            continue

    writer.close()

    if written_count == 0:
        raise RuntimeError("RDKit conversion produced 0 molecules for GNINA SDF")

    logger.info("Converted %d molecules to SDF using RDKit", written_count)
    return str(sdf_path.resolve()), None


def _get_gnina_environment(cfg, base_folder):
    """Get GNINA activation command and LD_LIBRARY_PATH."""
    gnina_config = cfg.get("gnina_config", {})
    env_path = gnina_config.get("env_path") or cfg.get("gnina_env_path")

    gnina_activate = cfg.get("gnina_activate")
    if not gnina_activate and env_path:
        conda_sh = cfg.get("conda_sh")
        if not conda_sh:
            # Auto-detect conda.sh from common installation paths
            for candidate in ("miniforge", "miniconda3", "mambaforge", "anaconda3"):
                path = Path(os.path.expanduser(f"~/{candidate}/etc/profile.d/conda.sh"))
                if path.exists():
                    conda_sh = str(path)
                    break
            if not conda_sh:
                conda_sh = os.path.expanduser("~/miniconda3/etc/profile.d/conda.sh")
        gnina_activate = f"source {conda_sh} && conda activate {env_path}"

    ld_library_path = cfg.get("gnina_ld_library_path")
    if not ld_library_path and env_path:
        env_path_obj = Path(env_path)
        if env_path_obj.exists():
            torch_libs = list(env_path_obj.glob("lib/python*/site-packages/torch/lib"))
            if torch_libs and (torch_libs[0] / "libcudnn.so.9").exists():
                ld_library_path = str(torch_libs[0])
            else:
                lib_dirs = list(env_path_obj.glob("lib"))
                for lib_dir in lib_dirs:
                    if (lib_dir / "libcudnn.so.9").exists() or list(
                        lib_dir.glob("libcudnn.so*")
                    ):
                        ld_library_path = str(lib_dir)
                        break

    # Auto-detect from PyTorch or conda when no env_path is configured
    if not ld_library_path:
        ld_library_path = _auto_detect_cudnn_path()

    return gnina_activate, ld_library_path


def _auto_detect_cudnn_path() -> str | None:
    """Auto-detect LD_LIBRARY_PATH for GNINA from PyTorch or conda environments."""
    # 1. Check PyTorch bundled libraries
    try:
        import torch

        torch_lib = Path(torch.__file__).parent / "lib"
        if torch_lib.is_dir() and list(torch_lib.glob("libcudnn.so*")):
            return str(torch_lib)
    except ImportError:
        pass

    # 2. Check common conda installation paths
    for root_name in ("miniforge", "miniconda3", "mambaforge", "anaconda3"):
        torch_glob = (
            Path.home()
            / root_name
            / "lib"
            / "python*"
            / "site-packages"
            / "torch"
            / "lib"
        )
        import glob as _glob

        matches = _glob.glob(str(torch_glob))
        for match in matches:
            match_path = Path(match)
            if match_path.is_dir() and list(match_path.glob("libcudnn.so*")):
                return str(match_path)

    return None


def _get_gnina_output_directory(cfg, base_folder):
    """Get GNINA output directory path."""
    gnina_config = cfg.get("gnina_config", {})
    cfg_out_dir = gnina_config.get("output_dir") or cfg.get("gnina_output_dir")
    if cfg_out_dir:
        out_dir_candidate = Path(cfg_out_dir)
        return (
            out_dir_candidate
            if out_dir_candidate.is_absolute()
            else (base_folder / out_dir_candidate)
        )
    return base_folder / "stages" / "05_docking" / "gnina"


def _write_file_wait_check(f, output_file, error_msg, prep_type="preparation"):
    """Write bash code to wait for file and check if it exists.

    Generates a polling loop that waits up to 5 minutes for a file to appear,
    with detailed error reporting if the file is not found.

    Args:
        f: File handle to write to
        output_file: Path to the expected output file
        error_msg: Error message to display if file not found
        prep_type: Label for log messages (e.g., 'preparation', 'ligand preparation')
    """
    f.write(f"# Check for {prep_type} output\n")
    f.write(f'echo "Checking for {prep_type} output: {output_file}"\n')
    f.write("# Wait for file to appear (preparation may take time)\n")
    f.write("max_wait=300  # 5 minutes max wait\n")
    f.write("wait_interval=2  # Check every 2 seconds\n")
    f.write("waited=0\n")
    f.write(f'while [ ! -f "{output_file}" ] && [ $waited -lt $max_wait ]; do\n')
    f.write("  sleep $wait_interval\n")
    f.write("  waited=$((waited + wait_interval))\n")
    f.write("  if [ $((waited % 10)) -eq 0 ]; then\n")
    f.write('    echo -n "."\n')
    f.write("  fi\n")
    f.write("done\n")
    f.write('echo ""\n')
    f.write(f'if [ ! -f "{output_file}" ]; then\n')
    f.write(f'  echo "{error_msg}"\n')
    f.write('  echo "Current directory: $(pwd)"\n')
    f.write('  echo "Listing files in current directory:"\n')
    f.write("  ls -la\n")
    f.write(f'  if [ -d "$(dirname "{output_file}")" ]; then\n')
    f.write('    echo "Listing files in output directory:"\n')
    f.write(f'    ls -la "$(dirname "{output_file}")"\n')
    f.write("  fi\n")
    f.write("  exit 1\n")
    f.write("fi\n")
    f.write(f'echo "{prep_type} output file found: {output_file}"\n')


def _create_gnina_script(
    ligands_dir,
    gnina_command_template,
    config_file,
    activate_cmd,
    ld_library_path,
    preparation_cmd,
    prepared_output_relative,
    protein_preparation_cmd,
    receptor,
):
    """Create GNINA run script (legacy single-batch mode)."""
    script_path = ligands_dir / "_workdir" / "run_gnina.sh"
    script_path.parent.mkdir(parents=True, exist_ok=True)
    with open(script_path, "w") as f:
        f.write("#!/usr/bin/env bash\n")
        f.write("set -eo pipefail\n")

        if activate_cmd:
            f.write(f"{activate_cmd}\n")
        if ld_library_path:
            f.write(f'export LD_LIBRARY_PATH="{ld_library_path}:$LD_LIBRARY_PATH"\n')

        if protein_preparation_cmd:
            _write_protein_prep_bash(f, protein_preparation_cmd, ligands_dir, "gnina")

        if preparation_cmd and prepared_output_relative:
            _write_ligand_prep_bash(f, preparation_cmd, prepared_output_relative)

        if receptor:
            _write_receptor_check_bash(f, receptor)

        f.write(f"cd {ligands_dir}\n")
        try:
            config_rel = config_file.relative_to(ligands_dir)
        except ValueError:
            config_rel = config_file.name
        f.write(f'echo "Starting GNINA docking with config: {config_rel}"\n')
        gnina_cmd = gnina_command_template.replace(
            "__GNINA_CONFIG__", shlex.quote(str(config_rel))
        )
        f.write(f"{gnina_cmd}\n")
        f.write("if [ $? -eq 0 ]; then\n")
        f.write('  echo "GNINA docking completed successfully"\n')
        f.write("else\n")
        f.write('  echo "GNINA docking failed with exit code $?"\n')
        f.write("  exit 1\n")
        f.write("fi\n")

    os.chmod(script_path, 0o755)
    return script_path


def _create_gnina_per_molecule_script(
    ligands_dir,
    gnina_command_template,
    activate_cmd,
    ld_library_path,
    protein_preparation_cmd,
    receptor,
    parallel_jobs: int,
):
    """Create GNINA run script with per-molecule processing and error handling.

    This script processes each molecule individually, allowing:
    - Individual molecule failures without stopping the entire run
    - Detailed per-molecule logging
    - Easy retry of failed molecules
    """
    script_path = ligands_dir / "_workdir" / "run_gnina.sh"
    script_path.parent.mkdir(parents=True, exist_ok=True)
    configs_dir = ligands_dir / "_workdir" / "configs"
    logs_dir = ligands_dir / "_workdir" / "gnina" / "logs"

    with open(script_path, "w") as f:
        f.write("#!/usr/bin/env bash\n")
        f.write("set -o pipefail\n\n")

        if activate_cmd:
            f.write(f"{activate_cmd}\n")
        if ld_library_path:
            f.write(f'export LD_LIBRARY_PATH="{ld_library_path}:$LD_LIBRARY_PATH"\n')
        if activate_cmd or ld_library_path:
            f.write("\n")

        if protein_preparation_cmd:
            _write_protein_prep_bash(f, protein_preparation_cmd, ligands_dir, "gnina")

        if receptor:
            _write_receptor_check_bash(f, receptor)

        f.write(f"cd {ligands_dir}\n")
        f.write(f'mkdir -p "{logs_dir}"\n\n')

        f.write("# Per-molecule docking with bounded parallelism\n")
        f.write(f"MAX_JOBS={max(1, int(parallel_jobs))}\n")
        f.write('echo "Running GNINA per-molecule with MAX_JOBS=${MAX_JOBS}"\n')
        f.write(f'STATUS_DIR="{ligands_dir}/_workdir/gnina/status"\n')
        f.write('rm -rf "${STATUS_DIR}"\n')
        f.write('mkdir -p "${STATUS_DIR}"\n\n')

        per_mol_cmd = gnina_command_template.replace("__GNINA_CONFIG__", '"${config}"')
        f.write("TOTAL=0\n")
        f.write(f'for config in "{configs_dir}"/gnina_*.ini; do\n')
        f.write('  [ -e "${config}" ] || continue\n')
        f.write('  while [ "$(jobs -rp | wc -l)" -ge "${MAX_JOBS}" ]; do\n')
        f.write("    sleep 0.2\n")
        f.write("  done\n")
        f.write("  TOTAL=$((TOTAL + 1))\n")
        f.write("  (\n")
        f.write('    mol_id=$(basename "${config}" .ini | sed "s/gnina_//")\n')
        f.write('    echo "[${TOTAL}] Processing ${mol_id}..."\n')
        f.write(f'    if {per_mol_cmd} 2>> "{logs_dir}/${{mol_id}}.log"; then\n')
        f.write('      echo "${mol_id}" >> "${STATUS_DIR}/success.txt"\n')
        f.write("    else\n")
        f.write("      exit_code=$?\n")
        f.write('      echo "${mol_id}" >> "${STATUS_DIR}/failed.txt"\n')
        f.write(
            '      echo "${mol_id}:${exit_code}" >> "${STATUS_DIR}/failed_with_code.txt"\n'
        )
        f.write("      if [ ${exit_code} -eq 143 ]; then\n")
        f.write(
            '        echo "${mol_id}: SIGTERM detected (timeout/memory limit/killed)" >> "${STATUS_DIR}/signals.txt"\n'
        )
        f.write("      fi\n")
        f.write("    fi\n")
        f.write("  ) &\n")
        f.write("done\n")
        f.write("wait\n\n")

        f.write("SUCCESS=0\n")
        f.write("FAILED=0\n")
        f.write('if [ -f "${STATUS_DIR}/success.txt" ]; then\n')
        f.write('  SUCCESS=$(wc -l < "${STATUS_DIR}/success.txt")\n')
        f.write("fi\n")
        f.write('if [ -f "${STATUS_DIR}/failed.txt" ]; then\n')
        f.write('  FAILED=$(wc -l < "${STATUS_DIR}/failed.txt")\n')
        f.write("fi\n")

        f.write('echo ""\n')
        f.write('echo "========================================"\n')
        f.write('echo "GNINA Docking Summary"\n')
        f.write('echo "========================================"\n')
        f.write('echo "Total molecules: ${TOTAL}"\n')
        f.write('echo "Successful: ${SUCCESS}"\n')
        f.write('echo "Failed: ${FAILED}"\n')
        f.write('echo "Parallel jobs: ${MAX_JOBS}"\n\n')

        f.write('if [ "${FAILED}" -gt 0 ] && [ -f "${STATUS_DIR}/failed.txt" ]; then\n')
        f.write(
            f'  cp "${{STATUS_DIR}}/failed.txt" "{ligands_dir}/gnina/failed_molecules.txt"\n'
        )
        f.write('  echo "Failed molecules were saved to failed_molecules.txt"\n')
        f.write("fi\n\n")

        f.write('if [ "${SUCCESS}" -gt 0 ]; then\n')
        f.write(
            '  echo "Docking completed with ${SUCCESS}/${TOTAL} molecules successful"\n'
        )
        f.write("  exit 0\n")
        f.write("fi\n")
        f.write('echo "ERROR: All molecules failed docking"\n')
        f.write("exit 1\n")

    os.chmod(script_path, 0o755)
    return script_path


def _extract_prepared_output_from_cmd(prep_cmd):
    """Extract prepared output path from ligand preparation command.

    Args:
        prep_cmd: Ligand preparation command containing '-osd <output>'

    Returns:
        Output path string or None if not found
    """
    if not prep_cmd:
        return None

    if isinstance(prep_cmd, (list, tuple)):
        parts = [str(p) for p in prep_cmd]
    else:
        parts = shlex.split(str(prep_cmd))

    if "-osd" not in parts:
        return None

    idx = parts.index("-osd")
    if idx + 1 < len(parts) and parts[idx + 1]:
        return str(parts[idx + 1])
    return None


def _resolve_receptor_path(receptor_pdb, base_folder=None):
    """Resolve receptor path to absolute, checking multiple locations.

    Args:
        receptor_pdb: Original receptor path from config
        base_folder: Base folder for relative path resolution

    Returns:
        Resolved Path object or None if not found
    """
    if not receptor_pdb:
        return None

    receptor_path = Path(receptor_pdb)
    if receptor_path.is_absolute() and receptor_path.exists():
        return receptor_path

    project_root = Path(__file__).parent.parent.parent.parent.parent

    if not receptor_path.is_absolute():
        candidate = (project_root / receptor_pdb).resolve()
        if candidate.exists():
            return candidate

    candidate = Path(receptor_pdb).resolve()
    if candidate.exists():
        return candidate

    return None


def _get_receptor_and_prep_cmd(cfg, ligands_dir, protein_preparation_tool, tool_name):
    """Get receptor path and protein preparation command.

    Args:
        cfg: Configuration dictionary
        ligands_dir: Directory for docking files
        protein_preparation_tool: Path to protein preparation tool (or None)
        tool_name: 'smina' or 'gnina' for logging

    Returns:
        Tuple of (receptor_path, protein_prep_cmd) or (None, None) on error
    """
    original_receptor = cfg.get("receptor_pdb")
    if not original_receptor:
        logger.error("%s: receptor_pdb is missing in config", tool_name.upper())
        return None, None

    receptor_path = _resolve_receptor_path(original_receptor)
    if not receptor_path:
        logger.error(
            "%s: Receptor file not found: %s", tool_name.upper(), original_receptor
        )
        return None, None

    receptor = str(receptor_path)

    if protein_preparation_tool is None:
        if "protein_prepared.pdb" in receptor:
            logger.info("%s: Using prepared receptor: %s", tool_name.upper(), receptor)
        else:
            logger.info("%s: Using receptor: %s", tool_name.upper(), receptor)
        return receptor, None

    prepared_receptor, protein_prep_cmd = _prepare_protein_for_docking(
        receptor, ligands_dir, protein_preparation_tool
    )
    if prepared_receptor != receptor:
        cfg["receptor_pdb"] = prepared_receptor
        logger.info(
            "%s: Using prepared protein: %s", tool_name.upper(), prepared_receptor
        )
    return prepared_receptor, protein_prep_cmd


def _generate_job_id(tool="dock"):
    """Generate a unique job ID."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    unique_id = uuid.uuid4().hex[:8]
    return f"{tool}_{timestamp}_{unique_id}"


def _save_job_metadata(
    ligands_dir,
    source_file,
    num_ligands,
    receptor_pdb,
    tools_prepared,
    scripts_prepared,
    ligands_csv,
    ligands_stats,
    job_ids,
    overall_job_id,
):
    """Save job metadata to JSON file."""
    ligands_dir.mkdir(parents=True, exist_ok=True)
    metadata = {
        "job_id": overall_job_id,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "source_file": str(source_file),
        "num_ligands": num_ligands,
        "receptor_pdb": receptor_pdb,
        "tools_prepared": tools_prepared,
        "scripts": scripts_prepared,
        "ligands_csv": str(ligands_csv),
        "ligands_counts": ligands_stats,
        "jobs": {
            tool: {
                "name": tool,
                "job_id": job_id,
                "script": str(ligands_dir / "_workdir" / f"run_{tool}.sh"),
            }
            for tool, job_id in job_ids.items()
        },
    }
    meta_path = ligands_dir / "job_meta.json"
    with open(meta_path, "w") as f:
        json.dump(metadata, f, indent=2)


def _save_job_ids(ligands_dir, overall_job_id, job_ids):
    """Save job IDs to a simple text file."""
    ids_path = ligands_dir / "job_ids.txt"
    try:
        with open(ids_path, "w") as f:
            f.write(f"overall: {overall_job_id}\n")
            f.write(f"smina: {job_ids.get('smina', '')}\n")
            f.write(f"gnina: {job_ids.get('gnina', '')}\n")
    except Exception as e:
        logger.warning("Failed to write job_ids.txt: %s", e)


def _update_metadata_with_run_status(ligands_dir, run_status):
    """Update job metadata with run status."""
    meta_path = ligands_dir / "job_meta.json"
    try:
        metadata = {}
        if meta_path.exists():
            with open(meta_path) as f:
                metadata = json.load(f)
        metadata["run_status"] = run_status
        with open(meta_path, "w") as f:
            json.dump(metadata, f, indent=2)
    except Exception as e:
        logger.warning("Failed to update metadata with run status: %s", e)


def _count_lines(path: Path) -> int:
    """Count lines in a text file; returns 0 if missing/unreadable."""
    try:
        if not path.exists():
            return 0
        return sum(1 for _ in path.open("r", encoding="utf-8", errors="ignore"))
    except OSError:
        return 0


def _count_smina_done(ligands_dir: Path) -> int:
    results_dir = ligands_dir / "_workdir" / "smina" / "results"
    if results_dir.exists():
        done = 0
        for p in results_dir.glob("*_out.sdf"):
            try:
                if p.is_file() and p.stat().st_size > 0:
                    done += 1
            except OSError:
                continue
        return done
    output_sdf = ligands_dir / "smina" / "smina_out.sdf"
    return 1 if output_sdf.exists() and output_sdf.stat().st_size > 0 else 0


def _count_gnina_done(ligands_dir: Path, output_sdf: Path) -> int:
    status_dir = ligands_dir / "_workdir" / "gnina" / "status"
    if status_dir.exists():
        return _count_lines(status_dir / "success.txt") + _count_lines(
            status_dir / "failed.txt"
        )
    return 1 if output_sdf.exists() and output_sdf.stat().st_size > 0 else 0


def _run_docking_script(
    script_path: Path, working_dir, log_path, background, tick=None
):
    """Run a docking script and return execution status."""
    if not script_path.exists():
        logger.error("Script not found: %s", script_path)
        return {"status": "error", "log_path": str(log_path)}

    if background:
        with open(log_path, "ab") as logf:
            subprocess.Popen(
                ["./" + script_path.name],
                stdout=logf,
                stderr=logf,
                cwd=str(working_dir),
            )
        logger.info("Started %s in background. Log: %s", script_path.name, log_path)
        return {"status": "started_background", "log_path": str(log_path)}
    else:
        with open(log_path, "wb") as logf:
            proc = subprocess.Popen(
                ["./" + script_path.name],
                stdout=logf,
                stderr=logf,
                cwd=str(working_dir),
            )
            while proc.poll() is None:
                if tick:
                    tick()
                time.sleep(0.5)
            if tick:
                tick()

        returncode = proc.returncode or 0
        if returncode == 0:
            logger.info(
                "%s completed successfully. Log: %s", script_path.name, log_path
            )
            return {"status": "completed", "log_path": str(log_path)}

        if returncode == 143:
            logger.error(
                "%s terminated by SIGTERM (exit 143) - likely timeout, memory limit, or killed by system. "
                "Check system logs and consider reducing batch size or using per-molecule mode.",
                script_path.name,
            )
        elif returncode == 137:
            logger.error(
                "%s killed by SIGKILL (exit 137) - likely OOM killer. "
                "Consider reducing batch size or using per-molecule mode.",
                script_path.name,
            )
        else:
            logger.error(
                "%s failed with exit code %d. See log: %s",
                script_path.name,
                returncode,
                log_path,
            )
        return {"status": "failed", "log_path": str(log_path), "exit_code": returncode}


def _run_smina(ligands_dir, background, job_id, tick=None):
    """Run SMINA docking."""
    workdir = ligands_dir / "_workdir"
    script_path = workdir / "run_smina.sh"
    log_path = workdir / "smina_run.log"
    status = _run_docking_script(script_path, workdir, log_path, background, tick=tick)
    if job_id:
        status["job_id"] = job_id

    # Aggregate per-molecule results if they exist
    results_dir = workdir / "smina" / "results"
    output_sdf = ligands_dir / "smina" / "smina_out.sdf"

    if not background and results_dir.exists():
        try:
            count = _aggregate_docking_results(results_dir, output_sdf)
            status["aggregated_molecules"] = count
            logger.info("Aggregated %d SMINA docking results", count)
        except Exception as e:
            logger.warning("Failed to aggregate SMINA results: %s", e)

    smina_out_dir = ligands_dir / "smina_results"
    status["results_dir"] = str(smina_out_dir) if smina_out_dir.exists() else None
    return status


def _run_gnina(ligands_dir, output_sdf, background, job_id, tick=None):
    """Run GNINA docking."""
    workdir = ligands_dir / "_workdir"
    script_path = workdir / "run_gnina.sh"
    log_path = workdir / "gnina_run.log"
    status = _run_docking_script(script_path, workdir, log_path, background, tick=tick)
    if job_id:
        status["job_id"] = job_id

    # Aggregate per-molecule results if they exist
    results_dir = workdir / "gnina" / "results"

    if not background and results_dir.exists():
        try:
            count = _aggregate_docking_results(results_dir, output_sdf)
            status["aggregated_molecules"] = count
            logger.info("Aggregated %d GNINA docking results", count)
        except Exception as e:
            logger.warning("Failed to aggregate GNINA results: %s", e)

    status["output"] = str(output_sdf)
    status["log"] = str(log_path)
    return status


def _parse_tools_config(cfg):
    """Parse tools configuration into a list of tool names."""
    tools_cfg = cfg.get("tools", "both")

    if isinstance(tools_cfg, str):
        tools_list = (
            [t.strip().lower() for t in tools_cfg.split(",")]
            if "," in tools_cfg
            else [tools_cfg.strip().lower()]
        )
    elif isinstance(tools_cfg, (list, tuple)):
        tools_list = [str(t).strip().lower() for t in tools_cfg]
    else:
        tools_list = ["both"]

    if "both" in tools_list or not tools_list:
        return ["smina", "gnina"]

    return [t for t in tools_list if t in ["smina", "gnina"]]


def _prepare_receptor_if_needed(
    cfg, ligands_dir, protein_preparation_tool, base_folder=None
):
    """Resolve receptor path and update config. Actual preparation happens in script."""
    original_receptor = cfg.get("receptor_pdb")
    if not original_receptor:
        return

    receptor_path = Path(original_receptor)
    if not receptor_path.is_absolute():
        project_root = Path(__file__).parent.parent.parent.parent.parent
        receptor_path = (project_root / original_receptor).resolve()
        if not receptor_path.exists():
            receptor_path = Path(original_receptor).resolve()

    if not receptor_path.exists():
        logger.warning(
            "Receptor file not found: %s (resolved to: %s)",
            original_receptor,
            receptor_path,
        )
        return

    cfg["receptor_pdb"] = str(receptor_path)


def _setup_smina(
    cfg,
    ligands_dir,
    ligands_csv,
    protein_preparation_tool,
    base_folder=None,
    ligand_preparation_tool=None,
):
    """Setup SMINA docking configuration and script with per-molecule processing."""
    try:
        if protein_preparation_tool is None:
            _prepare_receptor_if_needed(cfg, ligands_dir, None, base_folder)

        receptor, protein_prep_cmd = _get_receptor_and_prep_cmd(
            cfg, ligands_dir, protein_preparation_tool, "smina"
        )
        if not receptor:
            return None

        smina_config = cfg.get("smina_config", {})
        smina_bin_cfg = smina_config.get("bin") or cfg.get("smina_bin", "smina")
        smina_bin = _resolve_docking_binary(smina_bin_cfg, "smina")

        ligands_path, prep_cmd = _prepare_ligands_for_docking(
            ligands_csv, ligands_dir, ligand_preparation_tool, cfg, tool_name="smina"
        )

        smina_output_dir = ligands_dir / "smina"
        smina_output_dir.mkdir(parents=True, exist_ok=True)
        output_sdf = smina_output_dir / "smina_out.sdf"

        # Check if per-molecule mode is enabled (default: True)
        per_molecule_mode = cfg.get("per_molecule_docking", True)

        if per_molecule_mode and Path(ligands_path).exists():
            # Split SDF into per-molecule files
            molecules_dir = ligands_dir / "_workdir" / "molecules"
            molecule_files = _split_sdf_to_molecules(Path(ligands_path), molecules_dir)

            if molecule_files:
                # Create per-molecule configs
                _create_per_molecule_configs(
                    cfg, ligands_dir, receptor, molecule_files, "smina"
                )

                # Create per-molecule script
                script_path = _create_smina_per_molecule_script(
                    ligands_dir,
                    smina_bin,
                    protein_prep_cmd,
                )
                logger.info(
                    "SMINA per-molecule configuration prepared for %d molecules",
                    len(molecule_files),
                )
                return script_path
            else:
                logger.warning("No molecules found in SDF, falling back to batch mode")

        # Fallback to legacy batch mode
        config_file = ligands_dir / "_workdir" / "smina_config.ini"
        _create_smina_config_file(
            cfg, ligands_dir, receptor, ligands_path, config_file, output_sdf
        )

        prepared_output_relative = _extract_prepared_output_from_cmd(prep_cmd)
        script_path = _create_smina_script(
            ligands_dir,
            smina_bin,
            config_file,
            protein_prep_cmd,
            prep_cmd,
            prepared_output_relative,
        )
        logger.info("SMINA batch configuration prepared")
        return script_path
    except Exception as e:
        logger.error("Failed to setup SMINA: %s", e)
        import traceback

        logger.debug(traceback.format_exc())
        return None


def _setup_gnina(
    cfg,
    base_folder,
    ligands_dir,
    ligands_csv,
    ligand_preparation_tool,
    protein_preparation_tool,
):
    """Setup GNINA docking configuration and script with per-molecule processing."""
    try:
        original_receptor = cfg.get("receptor_pdb")
        if not original_receptor:
            logger.error("GNINA: receptor_pdb is missing in config")
            return None

        if "protein_prepared.pdb" in original_receptor:
            logger.warning(
                "GNINA: Config has prepared path: %s, this should have been restored",
                original_receptor,
            )
            project_root = Path(__file__).parent.parent.parent.parent.parent
            possible_originals = [
                project_root / "data/test/7EW9_apo.pdb",
            ]
            original_receptor = None
            for possible in possible_originals:
                if possible.exists():
                    original_receptor = str(possible.resolve())
                    break
            if not original_receptor:
                logger.error("GNINA: Could not find original receptor file")
                return None
            cfg["receptor_pdb"] = original_receptor

        receptor, protein_prep_cmd = _get_receptor_and_prep_cmd(
            cfg, ligands_dir, protein_preparation_tool, "gnina"
        )
        if not receptor:
            return None

        ligands_path, prep_cmd = _prepare_ligands_for_docking(
            ligands_csv, ligands_dir, ligand_preparation_tool, cfg, tool_name="gnina"
        )

        gnina_dir = _get_gnina_output_directory(cfg, base_folder)
        gnina_dir.mkdir(parents=True, exist_ok=True)
        output_sdf = gnina_dir / "gnina_out.sdf"

        gnina_config = cfg.get("gnina_config", {})
        gnina_bin_cfg = gnina_config.get("bin") or cfg.get("gnina_bin", "gnina")
        gnina_bin = _resolve_docking_binary(gnina_bin_cfg, "gnina")
        gnina_command_template = _build_gnina_command_template(
            cfg, gnina_bin, ligands_dir
        )
        activate_cmd, ld_library_path = _get_gnina_environment(cfg, base_folder)

        # Check if per-molecule mode is enabled (default: True)
        per_molecule_mode = cfg.get("per_molecule_docking", True)
        gnina_cpu_default = gnina_config.get("cpu", 8)
        cpu_per_process = _parse_positive_int(
            cfg.get("gnina_per_process_cpu", gnina_cpu_default), 8
        )
        parallel_jobs = _resolve_gnina_parallel_jobs(cfg, cpu_per_process)

        if per_molecule_mode:
            ligands_path_obj = Path(ligands_path)
            if prep_cmd:
                ready = _materialize_prepared_ligands(
                    prep_cmd,
                    ligands_path_obj,
                    ligands_dir,
                    tool_name="gnina",
                )
            else:
                ready = ligands_path_obj.exists()

            if not ready:
                logger.warning(
                    "GNINA per-molecule mode requested but prepared ligands are unavailable, falling back to batch mode"
                )
            else:
                logger.info(
                    "GNINA per-molecule mode enabled: cpu_per_process=%d, parallel_jobs=%d",
                    cpu_per_process,
                    parallel_jobs,
                )

                # Split SDF into per-molecule files
                molecules_dir = ligands_dir / "_workdir" / "molecules"
                molecule_files = _split_sdf_to_molecules(
                    ligands_path_obj, molecules_dir
                )

                if molecule_files:
                    # Create per-molecule configs
                    _create_per_molecule_configs(
                        cfg,
                        ligands_dir,
                        receptor,
                        molecule_files,
                        "gnina",
                        cpu_override=cpu_per_process,
                    )

                    # Create per-molecule script
                    script_path = _create_gnina_per_molecule_script(
                        ligands_dir,
                        gnina_command_template,
                        activate_cmd,
                        ld_library_path,
                        protein_prep_cmd,
                        receptor,
                        parallel_jobs,
                    )
                    logger.info(
                        "GNINA per-molecule configuration prepared for %d molecules",
                        len(molecule_files),
                    )
                    return script_path
                else:
                    logger.warning(
                        "No molecules found in SDF, falling back to batch mode"
                    )

        # Fallback to legacy batch mode
        config_file = ligands_dir / "_workdir" / "gnina_config.ini"
        _create_gnina_config_file(
            cfg, ligands_dir, receptor, ligands_path, output_sdf, config_file
        )

        prepared_output_relative = _extract_prepared_output_from_cmd(prep_cmd)
        script_path = _create_gnina_script(
            ligands_dir,
            gnina_command_template,
            config_file,
            activate_cmd,
            ld_library_path,
            prep_cmd,
            prepared_output_relative,
            protein_prep_cmd,
            receptor,
        )
        logger.info("GNINA batch configuration prepared")
        return script_path
    except Exception as e:
        logger.error("Failed to setup GNINA: %s", e)
        return None


def run_docking(config, reporter=None):
    """Main docking orchestration function."""
    cfg = load_config(config["config_docking"])
    if not cfg.get("run", False):
        logger.info("Docking disabled in config")
        return False

    base_folder = Path(config["folder_to_save"]).resolve()
    source = _find_latest_input_source(base_folder)
    if source is None:
        logger.warning(
            "No pass*SMILES.csv or sampled_molecules.csv found for docking input"
        )
        return False

    try:
        df = pd.read_csv(source)
    except Exception as e:
        logger.error("Failed to read docking input %s: %s", source, e)
        return False

    ligands_dir = base_folder / "stages" / "05_docking"
    ligands_csv = ligands_dir / "ligands.csv"

    try:
        ligands_stats = _prepare_ligands_dataframe(df, ligands_csv)
    except ValueError as e:
        logger.error("Ligand preparation failed: %s", e)
        return False

    ligand_preparation_tool = _validate_optional_tool_path(
        config.get("ligand_preparation_tool"), "Ligand preparation tool"
    )
    protein_preparation_tool = _validate_optional_tool_path(
        config.get("protein_preparation_tool"), "Protein preparation tool"
    )
    tools_list = _parse_tools_config(cfg)
    logger.info("Docking tools configured: %s", tools_list)

    # Heuristic warnings about common configuration issues.
    if "gnina" in tools_list:
        _warn_if_autobox_far_from_receptor(cfg, "gnina")
    if "smina" in tools_list:
        _warn_if_autobox_far_from_receptor(cfg, "smina")

    prepared_receptor_path = None
    if protein_preparation_tool:
        _prepare_receptor_if_needed(
            cfg, ligands_dir, protein_preparation_tool, base_folder
        )
        original_receptor = cfg.get("receptor_pdb")
        if original_receptor:
            original_receptor_path = Path(original_receptor)
            if original_receptor_path.exists():
                prepared_receptor_path, prep_cmd = _prepare_protein_for_docking(
                    original_receptor, ligands_dir, protein_preparation_tool
                )
                if prepared_receptor_path != original_receptor and prep_cmd:
                    try:
                        import time

                        _PROTEIN_PREP_TIMEOUT = 600  # 10 minutes
                        result = subprocess.run(
                            prep_cmd,
                            shell=False,
                            capture_output=True,
                            text=True,
                            cwd=str(ligands_dir),
                            timeout=_PROTEIN_PREP_TIMEOUT,
                        )
                        if result.returncode != 0:
                            stderr = (result.stderr or "").strip()
                            stdout = (result.stdout or "").strip()
                            if (
                                result.returncode == 127
                                or "not found" in stderr.lower()
                                or "no such file" in stderr.lower()
                            ):
                                logger.warning(
                                    "Protein preparation tool is unavailable at runtime. "
                                    "Skipping receptor preprocessing and using original receptor: %s",
                                    original_receptor,
                                )
                                cfg["receptor_pdb"] = str(original_receptor_path)
                            else:
                                logger.error(
                                    "Protein preparation command failed: %s",
                                    stderr or stdout,
                                )
                                return False
                        else:
                            prepared_path = Path(prepared_receptor_path)
                            max_wait = 300
                            wait_interval = 2
                            waited = 0
                            while waited < max_wait:
                                if (
                                    prepared_path.exists()
                                    and prepared_path.stat().st_size > 0
                                ):
                                    logger.info(
                                        "Protein prepared successfully: %s",
                                        prepared_receptor_path,
                                    )
                                    cfg["receptor_pdb"] = prepared_receptor_path
                                    break
                                time.sleep(wait_interval)
                                waited += wait_interval
                                if waited % 60 == 0:
                                    logger.info(
                                        "Waiting for protein preparation... (%ds)",
                                        waited,
                                    )

                            if not prepared_path.exists():
                                logger.error(
                                    "Protein preparation failed - output file not found after %ds: %s",
                                    max_wait,
                                    prepared_receptor_path,
                                )
                                logger.error("Command output: %s", result.stdout)
                                logger.error("Command error: %s", result.stderr)
                                return False
                            elif prepared_path.stat().st_size == 0:
                                logger.error(
                                    "Protein preparation failed - output file is empty: %s",
                                    prepared_receptor_path,
                                )
                                return False
                    except FileNotFoundError:
                        logger.warning(
                            "Protein preparation tool is unavailable at runtime. "
                            "Skipping receptor preprocessing and using original receptor: %s",
                            original_receptor,
                        )
                        cfg["receptor_pdb"] = str(original_receptor_path)
                    except Exception as e:
                        logger.error("Failed to prepare protein: %s", e)
                        import traceback

                        logger.debug(traceback.format_exc())
                        return False

    scripts_prepared = []
    job_ids = {}
    overall_job_id = _generate_job_id("dock")

    if "smina" in tools_list:
        script = _setup_smina(
            cfg,
            ligands_dir,
            ligands_csv,
            None,
            base_folder,
            ligand_preparation_tool=ligand_preparation_tool,
        )
        if script:
            scripts_prepared.append(str(script))
            job_ids["smina"] = _generate_job_id("smina")

    if "gnina" in tools_list:
        try:
            script = _setup_gnina(
                cfg,
                base_folder,
                ligands_dir,
                ligands_csv,
                ligand_preparation_tool,
                None,
            )
            if script:
                scripts_prepared.append(str(script))
                job_ids["gnina"] = _generate_job_id("gnina")
        except Exception as e:
            logger.warning("GNINA setup failed, continuing without GNINA: %s", e)
            tools_list = [t for t in tools_list if t != "gnina"]

    if not scripts_prepared:
        logger.error("No docking tools were successfully configured")
        return False
    try:
        original_receptor = cfg.get("receptor_pdb")
        _save_job_metadata(
            ligands_dir,
            source,
            len(df),
            original_receptor,
            list(job_ids.keys()),
            scripts_prepared,
            ligands_csv,
            ligands_stats,
            job_ids,
            overall_job_id,
        )
        _save_job_ids(ligands_dir, overall_job_id, job_ids)
        logger.info("Docking job ID: %s", overall_job_id)
    except Exception as e:
        logger.warning("Failed to save metadata: %s", e)

    auto_run = cfg.get("auto_run", True)
    if auto_run:
        background = bool(cfg.get("run_in_background", False))
        run_status = {}

        selected_tools = [
            t for t in tools_list if t in ["smina", "gnina"] and t in job_ids
        ]
        progress_total = 0
        tool_totals: dict[str, int] = {}
        gnina_output_sdf = (
            _get_gnina_output_directory(cfg, base_folder) / "gnina_out.sdf"
        )
        if reporter is not None and not background and selected_tools:
            configs_dir = ligands_dir / "_workdir" / "configs"
            for tool in selected_tools:
                count = 0
                if configs_dir.exists():
                    count = len(list(configs_dir.glob(f"{tool}_*.ini")))
                tool_totals[tool] = count if count > 0 else 1
            progress_total = sum(tool_totals.values()) or 1

            def _count_done() -> int:
                done = 0
                if "smina" in tool_totals:
                    done += min(_count_smina_done(ligands_dir), tool_totals["smina"])
                if "gnina" in tool_totals:
                    done += min(
                        _count_gnina_done(ligands_dir, gnina_output_sdf),
                        tool_totals["gnina"],
                    )
                return done

            def _make_tick(tool_name: str):
                def _tick() -> None:
                    reporter.progress(
                        _count_done(), progress_total, message=f"Docking ({tool_name})"
                    )

                return _tick

            reporter.progress(0, progress_total, message="Docking")

        if "smina" in job_ids:
            logger.info("Running SMINA docking")
            tick = (
                _make_tick("smina")
                if reporter is not None and not background and progress_total
                else None
            )
            run_status["smina"] = _run_smina(
                ligands_dir, background, job_ids["smina"], tick=tick
            )
            if not background and run_status["smina"].get("status") == "completed":
                log_path = ligands_dir / "_workdir" / "smina_run.log"
                _emit_post_docking_warnings("smina", log_path)

        if "gnina" in job_ids:
            logger.info("Running GNINA docking")
            try:
                gnina_dir = _get_gnina_output_directory(cfg, base_folder)
                output_sdf = gnina_dir / "gnina_out.sdf"
                tick = (
                    _make_tick("gnina")
                    if reporter is not None and not background and progress_total
                    else None
                )
                run_status["gnina"] = _run_gnina(
                    ligands_dir, output_sdf, background, job_ids["gnina"], tick=tick
                )
                if not background and run_status["gnina"].get("status") == "completed":
                    log_path = ligands_dir / "_workdir" / "gnina_run.log"
                    _emit_post_docking_warnings(
                        "gnina", log_path, output_sdf=output_sdf
                    )
            except Exception as e:
                logger.error("GNINA execution failed: %s", e)
                run_status["gnina"] = {"status": "failed", "error": str(e)}

        try:
            _update_metadata_with_run_status(ligands_dir, run_status)
        except Exception as e:
            logger.warning("Failed to update metadata with run status: %s", e)

        if not background:
            selected_tools = [t for t in tools_list if t in ["smina", "gnina"]]
            completed_tools = [
                t
                for t in selected_tools
                if run_status.get(t, {}).get("status") == "completed"
            ]
            failed_tools = [
                t
                for t in selected_tools
                if run_status.get(t, {}).get("status") == "failed"
            ]

            if reporter is not None and progress_total:
                reporter.progress(
                    progress_total, progress_total, message="Docking complete"
                )

            if failed_tools:
                logger.error("Docking tools failed: %s", ", ".join(failed_tools))

            if len(completed_tools) == len(selected_tools):
                return True
            elif len(completed_tools) > 0:
                logger.warning(
                    "Only %d/%d docking tools completed successfully",
                    len(completed_tools),
                    len(selected_tools),
                )
                return False
            else:
                logger.error("All docking tools failed")
                return False

    # If docking wasn't executed automatically, still emit warnings based on any
    # existing logs/outputs in the run folder. This makes reruns/debugging cheaper.
    if not auto_run:
        try:
            if "smina" in tools_list:
                _emit_post_docking_warnings(
                    "smina", ligands_dir / "_workdir" / "smina_run.log"
                )
            if "gnina" in tools_list:
                gnina_dir = _get_gnina_output_directory(cfg, base_folder)
                _emit_post_docking_warnings(
                    "gnina",
                    ligands_dir / "_workdir" / "gnina_run.log",
                    output_sdf=gnina_dir / "gnina_out.sdf",
                )
        except Exception:
            pass

    return True
