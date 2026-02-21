import os
import shlex
from pathlib import Path

from hedgehog.configs.logger import logger


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
