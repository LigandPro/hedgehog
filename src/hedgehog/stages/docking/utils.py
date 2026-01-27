import json
import os
import subprocess
import uuid
from datetime import datetime
from pathlib import Path

import datamol as dm
import pandas as pd

from hedgehog.configs.logger import load_config, logger
from hedgehog.utils.input_paths import find_latest_input_source


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

    tool_config = cfg.get(f"{tool_name}_config", {})

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
        skip_keys.update({"env_path", "ld_library_path", "activate", "output_dir"})

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
        cmd: Protein preparation command string
        ligands_dir: Directory for docking files

    Returns:
        Tuple of (output_filename, absolute_path) or (None, None) if not found
    """
    parts = cmd.split()
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
        protein_prep_cmd: Protein preparation command string
        ligands_dir: Directory for docking files
        tool_name: 'smina' or 'gnina' for log messages
    """
    f.write(f"cd {ligands_dir}\n")

    protein_output_file, protein_output_abs_path = _extract_pdb_output_from_cmd(
        protein_prep_cmd, ligands_dir
    )

    if protein_output_abs_path:
        f.write(f'mkdir -p "$(dirname "{protein_output_abs_path}")"\n')
        f.write(f'touch "{protein_output_abs_path}" 2>/dev/null || true\n')

    tool_label = f" for {tool_name.upper()}" if tool_name == "gnina" else ""
    f.write(f'echo "Running protein preparation{tool_label}..."\n')
    f.write(f"{protein_prep_cmd} || PREP_EXIT_CODE=$?\n")
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
    f.write(f'mkdir -p "$(dirname "{prepared_output_relative}")"\n')
    f.write('echo "Running ligand preparation..."\n')
    f.write(f"{preparation_cmd}\n")
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
    f.write(
        f'echo "Starting {tool_name.upper()} docking with config: {config_file.name}"\n'
    )
    f.write(f"{docking_bin} --config {config_file.name}\n")
    f.write("if [ $? -eq 0 ]; then\n")
    f.write(f'  echo "{tool_name.upper()} docking completed successfully"\n')
    f.write("else\n")
    f.write(f'  echo "{tool_name.upper()} docking failed with exit code $?"\n')
    f.write("  exit 1\n")
    f.write("fi\n")


def _create_smina_script(
    ligands_dir,
    smina_bin,
    config_file,
    protein_prep_cmd,
    preparation_cmd=None,
    prepared_output_relative=None,
):
    """Create SMINA run script."""
    script_path = ligands_dir / "run_smina.sh"
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

    prepared_output_path = ligands_dir / "protein_prepared.pdb"

    receptor_absolute = str(receptor_path.resolve())
    output_absolute = str(prepared_output_path.resolve())
    cmd_line = f"cd {ligands_dir} && {protein_preparation_tool} {receptor_absolute} {output_absolute}"

    logger.info("Protein preprocessing will be performed: %s", cmd_line)
    return str(prepared_output_path.resolve()), cmd_line


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
    if ligand_preparation_tool:
        ligands_val_lower = ligands_val.lower()
        if ligands_val_lower.endswith(".csv"):
            input_format = "-icsv"
        elif ligands_val_lower.endswith((".smi", ".ismi", ".cmi")):
            input_format = "-ismi"
        else:
            input_format = "-icsv"

        prepared_output_path = ligands_dir / f"prepared_for_{tool_name}.sdf"
        prepared_output_path.parent.mkdir(parents=True, exist_ok=True)
        prepared_output_relative = str(prepared_output_path.relative_to(ligands_dir))

        cmd_line = f"{ligand_preparation_tool} {input_format} {ligands_val} -osd {prepared_output_relative}"
        return str(prepared_output_path.resolve()), cmd_line
    else:
        ligands_path, _ = _convert_with_rdkit(ligands_csv, ligands_dir)
        return ligands_path, None


def _convert_with_rdkit(ligands_csv, ligands_dir):
    """Convert SMILES to SDF using RDKit as fallback.
    Assumes CSV contains smiles and name columns."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError as err:
        raise RuntimeError("RDKit not available for ligand conversion") from err

    sdf_path = ligands_dir / "ligands_prepared.sdf"
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
        conda_sh = cfg.get(
            "conda_sh", os.path.expanduser("~/miniconda3/etc/profile.d/conda.sh")
        )
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

    return gnina_activate, ld_library_path


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
    gnina_bin,
    config_file,
    activate_cmd,
    ld_library_path,
    preparation_cmd,
    prepared_output_relative,
    protein_preparation_cmd,
    receptor,
):
    """Create GNINA run script."""
    script_path = ligands_dir / "run_gnina.sh"
    with open(script_path, "w") as f:
        f.write("#!/usr/bin/env bash\n")
        f.write("set -eo pipefail\n")

        if activate_cmd:
            f.write(f"{activate_cmd}\n")
            if ld_library_path:
                f.write(
                    f'export LD_LIBRARY_PATH="{ld_library_path}:$LD_LIBRARY_PATH"\n'
                )

        if protein_preparation_cmd:
            _write_protein_prep_bash(f, protein_preparation_cmd, ligands_dir, "gnina")

        if preparation_cmd and prepared_output_relative:
            _write_ligand_prep_bash(f, preparation_cmd, prepared_output_relative)

        if receptor:
            _write_receptor_check_bash(f, receptor)

        _write_docking_command_bash(f, ligands_dir, gnina_bin, config_file, "gnina")

    os.chmod(script_path, 0o755)
    return script_path


def _extract_prepared_output_from_cmd(prep_cmd):
    """Extract prepared output path from ligand preparation command.

    Args:
        prep_cmd: Ligand preparation command containing '-osd <output>'

    Returns:
        Output path string or None if not found
    """
    if not prep_cmd or "-osd" not in prep_cmd:
        return None
    parts = prep_cmd.split()
    idx = parts.index("-osd")
    if idx + 1 < len(parts):
        return parts[idx + 1]
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
                "script": str(ligands_dir / f"run_{tool}.sh"),
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
        logger.warning(f"Failed to update metadata with run status: {e}")


def _run_docking_script(script_path: Path, working_dir, log_path, background):
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
            subprocess.run(
                ["./" + script_path.name],
                check=True,
                stdout=logf,
                stderr=logf,
                cwd=str(working_dir),
            )
        logger.info("%s completed successfully. Log: %s", script_path.name, log_path)
        return {"status": "completed", "log_path": str(log_path)}


def _run_smina(ligands_dir, background, job_id):
    """Run SMINA docking."""
    script_path = ligands_dir / "run_smina.sh"
    log_path = ligands_dir / "smina_run.log"
    status = _run_docking_script(script_path, ligands_dir, log_path, background)
    if job_id:
        status["job_id"] = job_id
    smina_out_dir = ligands_dir / "smina_results"
    status["results_dir"] = str(smina_out_dir) if smina_out_dir.exists() else None
    return status


def _run_gnina(ligands_dir, output_sdf, background, job_id):
    """Run GNINA docking."""
    script_path = ligands_dir / "run_gnina.sh"
    log_path = ligands_dir / "gnina_run.log"
    status = _run_docking_script(script_path, ligands_dir, log_path, background)
    if job_id:
        status["job_id"] = job_id
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
    cfg, ligands_dir, ligands_csv, protein_preparation_tool, base_folder=None
):
    """Setup SMINA docking configuration and script."""
    try:
        if protein_preparation_tool is None:
            _prepare_receptor_if_needed(cfg, ligands_dir, None, base_folder)

        receptor, protein_prep_cmd = _get_receptor_and_prep_cmd(
            cfg, ligands_dir, protein_preparation_tool, "smina"
        )
        if not receptor:
            return None

        smina_config = cfg.get("smina_config", {})
        smina_bin = smina_config.get("bin") or cfg.get("smina_bin", "smina")

        ligand_preparation_tool = cfg.get("ligand_preparation_tool")
        ligands_path, prep_cmd = _prepare_ligands_for_docking(
            ligands_csv, ligands_dir, ligand_preparation_tool, cfg, tool_name="smina"
        )

        smina_output_dir = ligands_dir / "smina"
        smina_output_dir.mkdir(parents=True, exist_ok=True)
        output_sdf = smina_output_dir / "smina_out.sdf"

        config_file = ligands_dir / "smina_config.ini"
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
        logger.info("SMINA configuration prepared")
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
    """Setup GNINA docking configuration and script."""
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
        gnina_bin = gnina_config.get("bin") or cfg.get("gnina_bin", "gnina")
        activate_cmd, ld_library_path = _get_gnina_environment(cfg, base_folder)

        config_file = ligands_dir / "gnina_config.ini"
        _create_gnina_config_file(
            cfg, ligands_dir, receptor, ligands_path, output_sdf, config_file
        )

        prepared_output_relative = _extract_prepared_output_from_cmd(prep_cmd)
        script_path = _create_gnina_script(
            ligands_dir,
            gnina_bin,
            config_file,
            activate_cmd,
            ld_library_path,
            prep_cmd,
            prepared_output_relative,
            protein_prep_cmd,
            receptor,
        )
        logger.info("GNINA configuration prepared")
        return script_path
    except Exception as e:
        logger.error("Failed to setup GNINA: %s", e)
        return None


def run_docking(config):
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

    ligand_preparation_tool = config.get("ligand_preparation_tool")
    protein_preparation_tool = config.get("protein_preparation_tool")
    tools_list = _parse_tools_config(cfg)
    logger.info("Docking tools configured: %s", tools_list)

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

                        result = subprocess.run(
                            prep_cmd,
                            shell=True,
                            capture_output=True,
                            text=True,
                            cwd=str(ligands_dir),
                        )
                        if result.returncode != 0:
                            logger.error(
                                "Protein preparation command failed: %s", result.stderr
                            )
                            return False

                        prepared_path = Path(prepared_receptor_path)
                        max_wait = 60
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
                                break
                            time.sleep(wait_interval)
                            waited += wait_interval
                            if waited % 60 == 0:
                                logger.info(
                                    "Waiting for protein preparation... (%ds)", waited
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

                        cfg["receptor_pdb"] = prepared_receptor_path
                    except Exception as e:
                        logger.error("Failed to prepare protein: %s", e)
                        import traceback

                        logger.debug(traceback.format_exc())
                        return False

    scripts_prepared = []
    job_ids = {}
    overall_job_id = _generate_job_id("dock")

    if "smina" in tools_list:
        script = _setup_smina(cfg, ligands_dir, ligands_csv, None, base_folder)
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

        if "smina" in job_ids:
            logger.info("Running SMINA docking")
            run_status["smina"] = _run_smina(ligands_dir, background, job_ids["smina"])

        if "gnina" in job_ids:
            logger.info("Running GNINA docking")
            try:
                gnina_dir = _get_gnina_output_directory(cfg, base_folder)
                output_sdf = gnina_dir / "gnina_out.sdf"
                run_status["gnina"] = _run_gnina(
                    ligands_dir, output_sdf, background, job_ids["gnina"]
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

            if failed_tools:
                logger.error(f"Docking tools failed: {', '.join(failed_tools)}")

            if len(completed_tools) == len(selected_tools):
                return True
            elif len(completed_tools) > 0:
                logger.warning(
                    f"Only {len(completed_tools)}/{len(selected_tools)} docking tools completed successfully"
                )
                return False
            else:
                logger.error("All docking tools failed")
                return False

    return True
