import os
import shlex
import shutil
import subprocess
from pathlib import Path

from hedgehog._constants import TOOL_GNINA, TOOL_MATCHA, TOOL_SMINA
from hedgehog.configs.logger import logger
from hedgehog.docking.binaries import _resolve_docking_binary
from hedgehog.docking.config_writer import (
    _create_gnina_config_file,
    _create_per_molecule_configs,
    _create_smina_config_file,
    _parse_bool_config,
)
from hedgehog.docking.ligand_prep import (
    _extract_prepared_output_from_cmd,
    _materialize_prepared_ligands,
    _parse_positive_int,
    _prepare_ligands_for_docking,
    _split_sdf_to_molecules,
)
from hedgehog.docking.metadata import _generate_job_id
from hedgehog.docking.paths import _emit_post_docking_warnings, _resolve_autobox_path
from hedgehog.docking.receptor_prep import (
    _get_receptor_and_prep_cmd,
    _prepare_receptor_if_needed,
    _restore_gnina_receptor,
)
from hedgehog.setup import ensure_matcha_checkout


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

    tool_label = f" for {tool_name.upper()}" if tool_name == TOOL_GNINA else ""
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
        if tool_name == TOOL_GNINA:
            success_msg += f": {protein_output_abs_path}"
        f.write(f'  echo "{success_msg}"\n')
        if tool_name == TOOL_GNINA:
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


def _build_docker_gnina_command(
    container_cfg: dict,
    ligands_dir: Path,
    no_gpu_enabled: bool,
    placeholder: str,
) -> str:
    """Build a Docker-wrapped GNINA command template."""
    container_bin = str(container_cfg.get("bin", "gnina")).strip() or "gnina"
    gpu_request = container_cfg.get("gpus", "all")
    mounts = container_cfg.get("mounts") or ["/mnt:/mnt", "/home:/home", "/tmp:/tmp"]

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
    cmd_parts.append(shlex.quote(str(container_cfg["image"])))
    cmd_parts.append(shlex.quote(container_bin))
    cmd_parts.append(f"--config {placeholder}")
    if no_gpu_enabled:
        cmd_parts.append("--no_gpu")
    return " ".join(cmd_parts)


def _build_gnina_command_template(cfg: dict, gnina_bin: str, ligands_dir: Path) -> str:
    """Build a GNINA command template with config placeholder.

    The returned command must contain the ``__GNINA_CONFIG__`` placeholder that
    is replaced at script generation time.
    """
    placeholder = "__GNINA_CONFIG__"
    gnina_cfg = cfg.get("gnina_config", {}) or {}
    no_gpu_enabled = _parse_bool_config(gnina_cfg.get("no_gpu"))
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

    return _build_docker_gnina_command(
        container_cfg, ligands_dir, no_gpu_enabled, placeholder
    )


def _project_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _resolve_executable(bin_name: str) -> str | None:
    """Resolve an executable path from config or PATH."""
    candidate = Path(bin_name).expanduser()
    if candidate.exists() and candidate.is_file() and os.access(candidate, os.X_OK):
        return str(candidate)
    return shutil.which(bin_name)


def _resolve_matcha_path(path_value: str | None, *extra_roots: Path) -> str | None:
    """Resolve a Matcha-related input path to an existing absolute path."""
    if not path_value:
        return None
    path = Path(str(path_value)).expanduser()
    if path.is_absolute() and path.exists():
        return str(path.resolve())
    for root in (*extra_roots, _project_root(), Path.cwd()):
        candidate = (root / path).resolve()
        if candidate.exists():
            return str(candidate)
    return None


def _resolve_matcha_autobox(cfg: dict, matcha_cfg: dict) -> str | None:
    """Resolve Matcha autobox path with fallback to GNINA/SMINA autobox config."""
    raw_autobox = (
        matcha_cfg.get("autobox_ligand")
        or (cfg.get("gnina_config", {}) or {}).get("autobox_ligand")
        or (cfg.get("smina_config", {}) or {}).get("autobox_ligand")
    )
    if not raw_autobox:
        return None
    resolved = _resolve_matcha_path(str(raw_autobox))
    if resolved:
        return resolved
    fallback = _resolve_autobox_path(str(raw_autobox), _project_root())
    return str(fallback) if fallback else None


def _resolve_matcha_center(
    cfg: dict, matcha_cfg: dict
) -> tuple[float, float, float] | None:
    """Resolve Matcha box center from Matcha/GNINA/SMINA configs."""
    for config_dict in (
        matcha_cfg,
        cfg.get("gnina_config", {}) or {},
        cfg.get("smina_config", {}) or {},
    ):
        center = config_dict.get("center")
        if isinstance(center, (list, tuple)) and len(center) >= 3:
            return (float(center[0]), float(center[1]), float(center[2]))
        cx = config_dict.get("center_x")
        cy = config_dict.get("center_y")
        cz = config_dict.get("center_z")
        if cx is not None and cy is not None and cz is not None:
            return (float(cx), float(cy), float(cz))
    return None


def _build_matcha_command(
    cfg: dict,
    ligands_dir: Path,
    receptor: str,
    ligands_path: str,
) -> tuple[list[str], str]:
    """Build Matcha CLI command and return (command_parts, run_name)."""
    matcha_cfg = cfg.get("matcha_config", {}) or {}
    matcha_repo = ensure_matcha_checkout(
        _project_root(),
        checkout_dir=matcha_cfg.get("checkout_dir"),
        pr_number=matcha_cfg.get("pr_number"),
    )
    uv_bin = str(matcha_cfg.get("uv_bin") or "uv").strip() or "uv"
    resolved_uv = _resolve_executable(uv_bin)
    if resolved_uv is None:
        raise FileNotFoundError(
            f"Matcha launcher '{uv_bin}' was not found. Install uv or configure matcha_config.uv_bin."
        )

    run_name = str(matcha_cfg.get("run_name") or "matcha_run").strip() or "matcha_run"
    output_dir = ligands_dir / "matcha"
    command = [
        resolved_uv,
        "run",
        "--project",
        str(matcha_repo),
        "matcha",
        "--receptor",
        receptor,
        "--ligand-dir",
        ligands_path,
        "--out",
        str(output_dir.resolve()),
        "--run-name",
        run_name,
        "--overwrite",
        "--device",
        str(matcha_cfg.get("device") or "auto"),
        "--n-samples",
        str(int(matcha_cfg.get("n_samples", 20))),
        "--scorer",
        str(matcha_cfg.get("scorer") or "gnina"),
    ]

    n_confs = matcha_cfg.get("n_confs")
    if n_confs is not None:
        command.extend(["--n-confs", str(int(n_confs))])

    checkpoints = _resolve_matcha_path(matcha_cfg.get("checkpoints"), matcha_repo)
    if checkpoints:
        command.extend(["--checkpoints", checkpoints])

    matcha_config_path = _resolve_matcha_path(matcha_cfg.get("config"), matcha_repo)
    if matcha_config_path:
        command.extend(["--config", matcha_config_path])

    scorer_path = _resolve_matcha_path(matcha_cfg.get("scorer_path"), matcha_repo)
    if not scorer_path:
        gnina_bin = ((cfg.get("gnina_config", {}) or {}).get("bin") or "").strip()
        scorer_path = _resolve_executable(gnina_bin) if gnina_bin else None
    if scorer_path:
        command.extend(["--scorer-path", scorer_path])

    if _parse_bool_config(matcha_cfg.get("scorer_minimize", True), True):
        command.append("--scorer-minimize")
    else:
        command.append("--no-scorer-minimize")

    if _parse_bool_config(matcha_cfg.get("physical_only", False), False):
        command.append("--physical-only")
    else:
        command.append("--keep-all-poses")

    if _parse_bool_config(matcha_cfg.get("keep_workdir", False), False):
        command.append("--keep-workdir")

    gpus = matcha_cfg.get("gpus")
    if gpus:
        command.extend(["--gpus", str(gpus)])

    autobox_path = _resolve_matcha_autobox(cfg, matcha_cfg)
    center = _resolve_matcha_center(cfg, matcha_cfg)
    if autobox_path:
        command.extend(["--autobox-ligand", autobox_path])
    elif center is not None:
        command.extend(
            [
                "--center-x",
                str(center[0]),
                "--center-y",
                str(center[1]),
                "--center-z",
                str(center[2]),
            ]
        )

    return command, run_name


def _create_matcha_script(
    ligands_dir: Path,
    matcha_command: list[str],
    run_name: str,
    preparation_cmd=None,
    prepared_output_relative=None,
):
    """Create Matcha run script."""
    script_path = ligands_dir / "_workdir" / "run_matcha.sh"
    script_path.parent.mkdir(parents=True, exist_ok=True)
    with open(script_path, "w") as f:
        f.write("#!/usr/bin/env bash\n")
        f.write("set -eo pipefail\n")
        f.write(f"cd {ligands_dir}\n")

        if preparation_cmd and prepared_output_relative:
            _write_ligand_prep_bash(f, preparation_cmd, prepared_output_relative)

        f.write(f'echo "Starting MATCHA docking run: {shlex.quote(run_name)}"\n')
        f.write(" ".join(shlex.quote(part) for part in matcha_command) + "\n")
        f.write("if [ $? -eq 0 ]; then\n")
        f.write('  echo "MATCHA docking completed successfully"\n')
        f.write("else\n")
        f.write('  echo "MATCHA docking failed with exit code $?"\n')
        f.write("  exit 1\n")
        f.write("fi\n")

    os.chmod(script_path, 0o700)
    return script_path


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

    os.chmod(
        script_path, 0o700
    )  # Owner-only: generated script, no need for group/world access
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

    os.chmod(
        script_path, 0o700
    )  # Owner-only: generated script, no need for group/world access
    return script_path


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

    os.chmod(
        script_path, 0o700
    )  # Owner-only: generated script, no need for group/world access
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

    os.chmod(
        script_path, 0o700
    )  # Owner-only: generated script, no need for group/world access
    return script_path


def _detect_conda_sh() -> str:
    """Auto-detect conda.sh from common installation paths."""
    for candidate in ("miniforge", "miniconda3", "mambaforge", "anaconda3"):
        path = Path(os.path.expanduser(f"~/{candidate}/etc/profile.d/conda.sh"))
        if path.exists():
            return str(path)
    return os.path.expanduser("~/miniconda3/etc/profile.d/conda.sh")


def _resolve_gnina_activate(cfg, env_path):
    """Resolve GNINA conda activation command."""
    gnina_activate = cfg.get("gnina_activate")
    if gnina_activate or not env_path:
        return gnina_activate

    conda_sh = cfg.get("conda_sh") or _detect_conda_sh()
    return f"source {conda_sh} && conda activate {env_path}"


def _resolve_env_library_paths(env_path):
    """Collect library paths from a conda/virtualenv environment directory."""
    env_lib_paths: list[str] = []
    env_path_obj = Path(env_path)
    if not env_path_obj.exists():
        return env_lib_paths

    torch_libs = sorted(env_path_obj.glob("lib/python*/site-packages/torch/lib"))
    env_lib_paths.extend(str(path) for path in torch_libs if path.is_dir())

    nvidia_libs = sorted(env_path_obj.glob("lib/python*/site-packages/nvidia/*/lib"))
    env_lib_paths.extend(str(path) for path in nvidia_libs if path.is_dir())

    lib_dir = env_path_obj / "lib"
    if lib_dir.is_dir():
        env_lib_paths.append(str(lib_dir))

    return env_lib_paths


def _get_gnina_environment(cfg, base_folder):
    """Get GNINA activation command and LD_LIBRARY_PATH."""
    gnina_config = cfg.get("gnina_config", {})
    env_path = gnina_config.get("env_path") or cfg.get("gnina_env_path")

    gnina_activate = _resolve_gnina_activate(cfg, env_path)

    ld_library_path = cfg.get("gnina_ld_library_path")
    if not ld_library_path and env_path:
        ld_library_path = _join_existing_library_paths(
            _resolve_env_library_paths(env_path)
        )

    # Auto-detect from PyTorch or conda when no env_path is configured
    if not ld_library_path:
        ld_library_path = _auto_detect_cudnn_path()

    return gnina_activate, ld_library_path


def _join_existing_library_paths(paths: list[str]) -> str | None:
    """Join existing library directories into LD_LIBRARY_PATH string."""
    resolved: list[str] = []
    seen: set[str] = set()
    for raw in paths:
        path = Path(raw).expanduser()
        if not path.is_dir():
            continue
        normalized = str(path.resolve())
        if normalized in seen:
            continue
        seen.add(normalized)
        resolved.append(normalized)
    if not resolved:
        return None
    return os.pathsep.join(resolved)


def _auto_detect_cudnn_path() -> str | None:
    """Auto-detect LD_LIBRARY_PATH for GNINA from PyTorch or conda environments."""
    try:
        from hedgehog.setup._gnina import _collect_gnina_library_paths
    except Exception:
        pass
        return None

    return _join_existing_library_paths(_collect_gnina_library_paths())


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


def _count_visible_nvidia_gpus() -> int:
    """Return count of visible NVIDIA GPUs available to this process."""
    visible_devices = os.environ.get("CUDA_VISIBLE_DEVICES", "").strip()
    if visible_devices:
        tokens = [item.strip() for item in visible_devices.split(",") if item.strip()]
        if any(
            token in {"-1", "none", "void"} for token in (t.lower() for t in tokens)
        ):
            return 0
        return len(tokens)

    nvidia_smi = shutil.which("nvidia-smi")
    if not nvidia_smi:
        return 0
    try:
        # Security: using list args (shell=False) to prevent command injection
        result = subprocess.run(
            [nvidia_smi, "-L"],
            capture_output=True,
            text=True,
            timeout=5,
        )
    except (OSError, subprocess.TimeoutExpired):
        return 0
    if result.returncode != 0:
        return 0
    lines = [line for line in result.stdout.splitlines() if line.strip()]
    return len(lines)


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


def _resolve_gnina_parallelism(cfg, gnina_config):
    """Resolve GNINA CPU-per-process, GPU count, and parallel jobs."""
    gnina_cpu_default = gnina_config.get("cpu", 8)
    cpu_per_process = _parse_positive_int(
        cfg.get("gnina_per_process_cpu", gnina_cpu_default), 8
    )
    no_gpu_enabled = _parse_bool_config(gnina_config.get("no_gpu", False))
    gpu_count = 0 if no_gpu_enabled else _count_visible_nvidia_gpus()
    parallel_jobs = _resolve_gnina_parallel_jobs(cfg, cpu_per_process)
    return cpu_per_process, gpu_count, parallel_jobs


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
        smina_bin_cfg = smina_config.get("bin") or cfg.get("smina_bin", TOOL_SMINA)
        smina_bin = _resolve_docking_binary(smina_bin_cfg, TOOL_SMINA)

        ligands_path, prep_cmd = _prepare_ligands_for_docking(
            ligands_csv, ligands_dir, ligand_preparation_tool, cfg, tool_name=TOOL_SMINA
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


def _setup_gnina_per_molecule(
    cfg,
    ligands_dir,
    ligands_path,
    prep_cmd,
    receptor,
    protein_prep_cmd,
    gnina_command_template,
    activate_cmd,
    ld_library_path,
    cpu_per_process,
    parallel_jobs,
    gpu_count,
):
    """Setup GNINA per-molecule docking. Returns script path or None on fallback."""
    ligands_path_obj = Path(ligands_path)
    if prep_cmd:
        ready = _materialize_prepared_ligands(
            prep_cmd, ligands_path_obj, ligands_dir, tool_name=TOOL_GNINA
        )
    else:
        ready = ligands_path_obj.exists()

    if not ready:
        logger.warning(
            "GNINA per-molecule mode requested but prepared ligands are unavailable, falling back to batch mode"
        )
        return None

    logger.info(
        "GNINA per-molecule mode enabled: cpu_per_process=%d, parallel_jobs=%d, gpu_count=%d",
        cpu_per_process,
        parallel_jobs,
        gpu_count,
    )

    molecules_dir = ligands_dir / "_workdir" / "molecules"
    molecule_files = _split_sdf_to_molecules(ligands_path_obj, molecules_dir)
    if not molecule_files:
        logger.warning("No molecules found in SDF, falling back to batch mode")
        return None

    _create_per_molecule_configs(
        cfg,
        ligands_dir,
        receptor,
        molecule_files,
        "gnina",
        cpu_override=cpu_per_process,
    )
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


def _setup_gnina_batch_mode(
    cfg,
    ligands_dir,
    receptor,
    ligands_path,
    prep_cmd,
    output_sdf,
    gnina_command_template,
    activate_cmd,
    ld_library_path,
    protein_prep_cmd,
):
    """Setup GNINA in legacy batch mode. Returns script path."""
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
        original_receptor = _restore_gnina_receptor(cfg)
        if not original_receptor:
            return None

        receptor, protein_prep_cmd = _get_receptor_and_prep_cmd(
            cfg, ligands_dir, protein_preparation_tool, "gnina"
        )
        if not receptor:
            return None

        ligands_path, prep_cmd = _prepare_ligands_for_docking(
            ligands_csv, ligands_dir, ligand_preparation_tool, cfg, tool_name=TOOL_GNINA
        )

        gnina_dir = _get_gnina_output_directory(cfg, base_folder)
        gnina_dir.mkdir(parents=True, exist_ok=True)
        output_sdf = gnina_dir / "gnina_out.sdf"

        gnina_config = cfg.get("gnina_config", {})
        gnina_bin_cfg = gnina_config.get("bin") or cfg.get("gnina_bin", TOOL_GNINA)
        gnina_bin = _resolve_docking_binary(gnina_bin_cfg, TOOL_GNINA)
        gnina_command_template = _build_gnina_command_template(
            cfg, gnina_bin, ligands_dir
        )
        activate_cmd, ld_library_path = _get_gnina_environment(cfg, base_folder)

        cpu_per_process, gpu_count, parallel_jobs = _resolve_gnina_parallelism(
            cfg, gnina_config
        )

        if cfg.get("per_molecule_docking", True):
            script = _setup_gnina_per_molecule(
                cfg,
                ligands_dir,
                ligands_path,
                prep_cmd,
                receptor,
                protein_prep_cmd,
                gnina_command_template,
                activate_cmd,
                ld_library_path,
                cpu_per_process,
                parallel_jobs,
                gpu_count,
            )
            if script:
                return script

        return _setup_gnina_batch_mode(
            cfg,
            ligands_dir,
            receptor,
            ligands_path,
            prep_cmd,
            output_sdf,
            gnina_command_template,
            activate_cmd,
            ld_library_path,
            protein_prep_cmd,
        )
    except Exception as e:
        logger.error("Failed to setup GNINA: %s", e)
        return None


def _setup_docking_tools(
    cfg, tools_list, base_folder, ligands_dir, ligands_csv, ligand_preparation_tool
):
    """Configure docking tools and return (scripts_prepared, job_ids, tools_list).

    tools_list may be modified in-place if GNINA setup fails.
    """
    scripts_prepared = []
    job_ids = {}

    if TOOL_SMINA in tools_list:
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
            job_ids[TOOL_SMINA] = _generate_job_id(TOOL_SMINA)

    if TOOL_GNINA in tools_list:
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
                job_ids[TOOL_GNINA] = _generate_job_id(TOOL_GNINA)
        except Exception as e:
            logger.warning("GNINA setup failed, continuing without GNINA: %s", e)
            tools_list[:] = [t for t in tools_list if t != TOOL_GNINA]

    if TOOL_MATCHA in tools_list:
        try:
            receptor, _ = _get_receptor_and_prep_cmd(
                cfg, ligands_dir, None, TOOL_MATCHA
            )
            if receptor:
                ligands_path, prep_cmd = _prepare_ligands_for_docking(
                    ligands_csv,
                    ligands_dir,
                    ligand_preparation_tool,
                    cfg,
                    tool_name=TOOL_MATCHA,
                )
                matcha_command, run_name = _build_matcha_command(
                    cfg, ligands_dir, receptor, ligands_path
                )
                script = _create_matcha_script(
                    ligands_dir,
                    matcha_command,
                    run_name,
                    prep_cmd,
                    _extract_prepared_output_from_cmd(prep_cmd),
                )
                scripts_prepared.append(str(script))
                job_ids[TOOL_MATCHA] = _generate_job_id(TOOL_MATCHA)
        except Exception as e:
            logger.warning("Matcha setup failed, continuing without Matcha: %s", e)
            tools_list[:] = [t for t in tools_list if t != TOOL_MATCHA]

    return scripts_prepared, job_ids


def _emit_manual_mode_warnings(cfg, tools_list, ligands_dir, base_folder):
    """Emit post-docking warnings in manual (non-auto_run) mode."""
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
