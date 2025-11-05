import json
import os
import shlex
import subprocess
import time
import uuid
from collections.abc import Sequence
from datetime import datetime
from pathlib import Path
from typing import Any

import datamol as dm
import pandas as pd

from hedge.configs.logger import load_config, logger


def _find_latest_input_source(base_folder: Path) -> Path | None:
    """Find the most recent input source file for docking."""
    candidates = [
        base_folder / "Synthesis" / "passSynthesisSMILES.csv",
        base_folder / "StructFilters" / "passStructFiltersSMILES.csv",
        base_folder / "Descriptors" / "passDescriptorsSMILES.csv",
        base_folder / "sampledMols.csv",
    ]
    for path in candidates:
        if path.exists():
            return path
    return None


def _prepare_ligands_dataframe(df: pd.DataFrame, output_csv: Path) -> dict[str, Any]:
    """Prepare ligands CSV from input DataFrame with SMILES validation."""
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    names = df["mol_idx"].astype(str)
    rows = []
    skipped_smiles = []

    for idx, row in df.iterrows():
        smi = str(row["smiles"])
        try:
            mol = dm.to_mol(smi)
            if mol is None:
                skipped_smiles.append(smi)
                continue
        except (ValueError, TypeError, RuntimeError):
            skipped_smiles.append(smi)
            continue

        model_name = str(row["model_name"])
        mol_idx = str(row["mol_idx"])

        rows.append(
            {
                "smiles": smi,
                "name": str(names.iloc[idx]),
                "model_name": model_name,
                "mol_idx": mol_idx,
            }
        )

    output_df = pd.DataFrame(rows, columns=["smiles", "name", "model_name", "mol_idx"])
    output_df.to_csv(output_csv, index=False)

    if skipped_smiles:
        skip_path = output_csv.parent / "skipped_smiles.txt"
        with skip_path.open("w") as f:
            for smi in skipped_smiles:
                f.write(f"{smi}\n")
        logger.warning(
            f"Some SMILES could not be parsed for docking: {len(skipped_smiles)}/{len(df)}. "
            f"See {skip_path}"
        )
    return {
        "csv_path": str(output_csv),
        "total": len(df),
        "written": len(rows),
        "skipped": len(skipped_smiles),
    }


def _load_smina_ini_content(cfg: dict[str, Any]) -> list[str]:
    """Load SMINA INI configuration content."""
    inline_content = cfg.get("smina_ini_content")
    if inline_content and isinstance(inline_content, str):
        return inline_content.splitlines()

    ini_path = Path(cfg.get("smina_ini", ""))
    if ini_path.exists():
        try:
            with ini_path.open() as f:
                return f.read().splitlines()
        except (OSError, IOError) as e:
            logger.warning("Failed to read smina_ini file %s: %s", ini_path, e)

    return [
        "output-dir = smina_results",
        "screen-type = smina",
        'metadata-template = {"software": "smina"}',
    ]


def _update_ini_key_value(lines: list[str], key: str, value: str) -> list[str]:
    """Update or add a key-value pair in INI file lines."""
    updated_lines = []
    found = False
    for line in lines:
        if line.strip().startswith(f"{key} "):
            updated_lines.append(f"{key} = {value}")
            found = True
        else:
            updated_lines.append(line)
    if not found:
        updated_lines.append(f"{key} = {value}")
    return updated_lines


def _create_smina_config(
    cfg: dict[str, Any], ligands_csv: Path, output_ini: Path
) -> None:
    """Create SMINA INI configuration file."""
    lines = _load_smina_ini_content(cfg)

    if not any(line.strip().startswith("screen-type") for line in lines):
        lines.insert(0, "screen-type = smina")

    lines = _update_ini_key_value(lines, "input-files", f"[{ligands_csv.name}]")
    lines = _update_ini_key_value(lines, "smiles-col", "0")
    lines = _update_ini_key_value(lines, "name-col", "1")

    receptor = cfg.get("receptor_pdb")
    if receptor:
        lines = _update_ini_key_value(lines, "receptors", f"[{receptor}]")

    output_ini.parent.mkdir(parents=True, exist_ok=True)
    with output_ini.open("w") as f:
        f.write("\n".join(lines) + "\n")


def _create_smina_script(
    ligands_dir: Path,
    ini_file: Path,
    activate_cmd: str | None,
    protein_prep_cmd: str | None,
) -> Path:
    """Create SMINA run script."""
    script_path = ligands_dir / "run_smina.sh"
    with script_path.open("w") as f:
        f.write("#!/usr/bin/env bash\n")
        f.write("set -eo pipefail\n")
        if activate_cmd:
            f.write(f"{activate_cmd}\n")

        if protein_prep_cmd:
            f.write(f"cd {ligands_dir}\n")
            parts = protein_prep_cmd.split()
            protein_output_file = None
            protein_output_abs_path = None
            for part in reversed(parts):
                if ".pdb" in part:
                    if part.startswith("/"):
                        protein_output_file = part.split("/")[-1]
                        protein_output_abs_path = part
                    else:
                        protein_output_file = part
                        protein_output_abs_path = str(ligands_dir / part)
                    break

            if protein_output_abs_path:
                f.write(f'mkdir -p "$(dirname "{protein_output_abs_path}")"\n')
                f.write(f'touch "{protein_output_abs_path}" 2>/dev/null || true\n')

            f.write(
                'echo "Running protein preparation..."\n'
                f"{protein_prep_cmd} || PREP_EXIT_CODE=$?\n"
                'if [ ! -z "$PREP_EXIT_CODE" ]; then\n'
                '  echo "WARNING: Protein preparation command exited with code $PREP_EXIT_CODE"\n'
                "fi\n"
            )

            if protein_output_abs_path:
                f.write(
                    "# If file exists in current dir but not at absolute path, move it\n"
                    f'if [ -f "{protein_output_file}" ] && [ ! -f "{protein_output_abs_path}" ]; then\n'
                    f'  mv "{protein_output_file}" "{protein_output_abs_path}"\n'
                    "fi\n"
                    "# Check if prepared file exists and is valid (not empty)\n"
                    f'if [ ! -f "{protein_output_abs_path}" ] && [ ! -f "{protein_output_file}" ]; then\n'
                    f'  echo "ERROR: Protein preparation failed - output file not found at {protein_output_abs_path}"\n'
                    '  echo "Current directory: $(pwd)"\n'
                    '  echo "Listing files in current directory:"\n'
                    "  ls -la\n"
                    "  exit 1\n"
                    f'elif [ -f "{protein_output_abs_path}" ] && [ ! -s "{protein_output_abs_path}" ]; then\n'
                    f'  echo "ERROR: Protein preparation failed - output file is empty: {protein_output_abs_path}"\n'
                    "  exit 1\n"
                    f'elif [ -f "{protein_output_file}" ] && [ ! -s "{protein_output_file}" ]; then\n'
                    f'  echo "ERROR: Protein preparation failed - output file is empty: {protein_output_file}"\n'
                    "  exit 1\n"
                    "else\n"
                    '  echo "Protein preparation completed successfully"\n'
                    "fi\n"
                )

        f.write(f"cd {ligands_dir}\npyscreener --config {ini_file.name}\n")
    script_path.chmod(0o755)  # noqa: S103
    return script_path


def _prepare_protein_for_docking(
    receptor_pdb: str,
    ligands_dir: Path,
    protein_preparation_tool: str,
) -> tuple[str, str | None]:
    """Prepare protein file for docking using external tool.

    Args:
        receptor_pdb: Path to the original receptor PDB file (must exist)
        ligands_dir: Directory where docking files are prepared
        protein_preparation_tool: Path to protein preparation tool

    Returns
    -------
        Tuple of (prepared_receptor_path, preparation_cmd) or (original_path, None).
        Note: prepared_receptor_path is where the file WILL be after preparation,
        not where it currently is.
    """
    if not protein_preparation_tool:
        return receptor_pdb, None

    receptor_path = Path(receptor_pdb)
    if not receptor_path.exists():
        logger.warning(
            f"Original receptor file not found: {receptor_pdb} (resolved to: {receptor_path}), skipping protein preprocessing"
        )
        return receptor_pdb, None

    prepared_output_path = ligands_dir / "protein_prepared.pdb"

    receptor_absolute = str(receptor_path.resolve())
    output_absolute = str(prepared_output_path.resolve())
    cmd_line = f"cd {ligands_dir} && {protein_preparation_tool} {receptor_absolute} {output_absolute}"

    logger.info(f"Protein preprocessing will be performed: {cmd_line}")
    return str(prepared_output_path.resolve()), cmd_line


def _prepare_ligands_for_gnina(
    ligands_csv: Path,
    ligands_dir: Path,
    ligand_preparation_tool: str | None,
    cfg: dict[str, Any],
) -> tuple[str, str | None]:
    """Prepare ligands file for GNINA docking."""
    ligands_arg = cfg.get("gnina_ligands")
    if ligands_arg:
        ligands_val = str(ligands_arg)
    else:
        ligands_val = str(ligands_csv.relative_to(ligands_dir))

    sdf_extensions = (".sdf", ".sdf.gz", ".osd", ".mol2")
    needs_conversion = not ligands_val.lower().endswith(sdf_extensions)

    if not needs_conversion:
        return ligands_val, None

    if ligand_preparation_tool:
        ligands_val_lower = ligands_val.lower()
        if ligands_val_lower.endswith(".csv"):
            input_format = "-icsv"
        elif ligands_val_lower.endswith((".smi", ".ismi", ".cmi")):
            input_format = "-ismi"
        else:
            input_format = "-icsv"

        prepared_output_path = ligands_dir / "docking" / "prepared_for_gnina.sdf"
        prepared_output_path.parent.mkdir(parents=True, exist_ok=True)
        prepared_output_relative = str(prepared_output_path.relative_to(ligands_dir))

        cmd_line = f"{ligand_preparation_tool} {input_format} {ligands_val} -osd {prepared_output_relative}"
        return str(prepared_output_path.resolve()), cmd_line
    return _convert_with_rdkit(ligands_csv, ligands_dir)


def _convert_with_rdkit(ligands_csv: Path, ligands_dir: Path) -> tuple[str, None]:
    """Convert SMILES to SDF using RDKit as fallback.

    Assumes CSV contains smiles and name columns.
    """
    try:
        from rdkit import Chem  # noqa: PLC0415
        from rdkit.Chem import AllChem  # noqa: PLC0415
    except ImportError:
        msg = "RDKit not available for ligand conversion"
        raise RuntimeError(msg) from None

    sdf_path = ligands_dir / "ligands_prepared.sdf"
    df = pd.read_csv(ligands_csv)
    smiles_series = df["smiles"]
    name_series = df["name"]

    writer = Chem.SDWriter(str(sdf_path))
    written_count = 0

    for smi, name in zip(smiles_series.astype(str), name_series.astype(str)):
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            mol = Chem.AddHs(mol)
            try:
                AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                AllChem.UFFOptimizeMolecule(mol)
            except (ValueError, RuntimeError):
                pass
            mol.SetProp("_Name", name)
            writer.write(mol)
            written_count += 1
        except (ValueError, TypeError, RuntimeError):
            continue

    writer.close()

    if written_count == 0:
        raise RuntimeError("RDKit conversion produced 0 molecules for GNINA SDF")

    logger.info(f"Converted {written_count} molecules to SDF using RDKit")
    return str(sdf_path.resolve()), None


def _build_gnina_command(  # noqa: C901
    cfg: dict[str, Any],
    receptor: str,
    ligands_path: str,
    output_sdf: Path,
) -> list[str]:
    """Build GNINA command line arguments."""
    user_path = cfg.get("gnina_path")
    user_args = cfg.get("gnina_args")
    if user_path and user_args:
        formatted = user_args.format(
            receptor=receptor, ligands=ligands_path, out=str(output_sdf)
        )
        return [str(user_path)] + shlex.split(formatted)

    gnina_bin = cfg.get("gnina_bin", "gnina")
    cmd = [gnina_bin, "-r", receptor, "-l", ligands_path, "-o", str(output_sdf)]

    center = cfg.get("center")
    if isinstance(center, (list, tuple)) and len(center) >= 3:  # noqa: PLR2004
        cmd.extend(
            [
                "--center_x",
                str(center[0]),
                "--center_y",
                str(center[1]),
                "--center_z",
                str(center[2]),
            ]
        )

    size = cfg.get("size")
    if isinstance(size, (list, tuple)) and len(size) >= 3:  # noqa: PLR2004
        cmd.extend(
            [
                "--size_x",
                str(size[0]),
                "--size_y",
                str(size[1]),
                "--size_z",
                str(size[2]),
            ]
        )

    autobox_ligand = cfg.get("gnina_autobox_ligand")
    if autobox_ligand:
        autobox_ligand_path = Path(autobox_ligand)
        if not autobox_ligand_path.is_absolute():
            project_root = Path(__file__).parent.parent.parent.parent.parent
            autobox_ligand_path = (project_root / autobox_ligand).resolve()
            if not autobox_ligand_path.exists():
                autobox_ligand_path = Path(autobox_ligand).resolve()

        if not autobox_ligand_path.exists():
            raise FileNotFoundError(
                f"GNINA: autobox_ligand file not found: {autobox_ligand} (resolved to: {autobox_ligand_path}). "
                f"Please check that the file exists or remove gnina_autobox_ligand from config if not needed."
            )
        cmd.extend(["--autobox_ligand", str(autobox_ligand_path)])

    optional_params = {
        "exhaustiveness": "--exhaustiveness",
        "num_modes": "--num_modes",
        "gnina_autobox_add": "--autobox_add",
        "gnina_cpu": "--cpu",
        "gnina_seed": "--seed",
    }

    for cfg_key, arg_name in optional_params.items():
        value = cfg.get(cfg_key)
        if value is not None:
            cmd.extend([arg_name, str(value)])

    extra_args = cfg.get("gnina_extra_args")
    if extra_args:
        extra = str(extra_args).strip()
        if extra:
            cmd.extend(extra.split())

    return cmd


def _get_gnina_environment(
    cfg: dict[str, Any],
    base_folder: Path,  # noqa: ARG001
) -> tuple[str | None, str | None]:
    """Get GNINA activation command and LD_LIBRARY_PATH."""
    gnina_activate = cfg.get("gnina_activate")
    if not gnina_activate:
        env_path = cfg.get("gnina_env_path")
        if env_path:
            conda_sh = cfg.get(
                "conda_sh",
                os.path.expanduser("~/miniconda3/etc/profile.d/conda.sh"),  # noqa: PTH111
            )
            gnina_activate = f"source {conda_sh} && conda activate {env_path}"

    ld_library_path = cfg.get("gnina_ld_library_path")
    if not ld_library_path:
        env_path = cfg.get("gnina_env_path")
        if env_path:
            env_path_obj = Path(env_path)
            torch_libs = list(env_path_obj.glob("lib/python*/site-packages/torch/lib"))
            if torch_libs and (torch_libs[0] / "libcudnn.so.9").exists():
                ld_library_path = str(torch_libs[0])

    return gnina_activate, ld_library_path


def _get_gnina_output_directory(
    cfg: dict[str, Any],
    base_folder: Path,
) -> Path:
    """Get GNINA output directory path."""
    cfg_out_dir = cfg.get("gnina_output_dir")
    if cfg_out_dir:
        out_dir_candidate = Path(cfg_out_dir)
        return (
            out_dir_candidate
            if out_dir_candidate.is_absolute()
            else (base_folder / out_dir_candidate)
        )
    return base_folder / "Docking" / "gnina_results"


def _write_file_wait_check(
    f: Any,  # noqa: ANN401
    output_file: str,
    error_msg: str,
    prep_type: str = "preparation",
) -> None:
    """Write bash code to wait for file and check if it exists."""
    f.write(
        f"# Check for {prep_type} output (prep is synchronous, should exist immediately)\n"
        f'echo "Checking for {prep_type} output: {output_file}"\n'
        f'if [ -f "{output_file}" ]; then\n'
        f'  echo "{prep_type} output file found: {output_file}"\n'
        "else\n"
        f'  echo "Waiting for {prep_type} output file (max 10 seconds)..."\n'
        "  max_wait=10  # 10 seconds max wait (prep should be synchronous)\n"
        "  wait_interval=1  # Check every 1 second\n"
        "  waited=0\n"
        f'  while [ ! -f "{output_file}" ] && [ $waited -lt $max_wait ]; do\n'
        "    sleep $wait_interval\n"
        "    waited=$((waited + wait_interval))\n"
        '    echo -n "."\n'
        "  done\n"
        '  echo ""\n'
        "fi\n"
        f'if [ ! -f "{output_file}" ]; then\n'
        f'  echo "{error_msg}"\n'
        '  echo "Current directory: $(pwd)"\n'
        '  echo "Listing files in current directory:"\n'
        "  ls -la\n"
        f'  if [ -d "$(dirname "{output_file}")" ]; then\n'
        '    echo "Listing files in output directory:"\n'
        f'    ls -la "$(dirname "{output_file}")"\n'
        "  fi\n"
        "  exit 1\n"
        "fi\n"
        f'echo "{prep_type} output file found: {output_file}"\n'
    )


def _create_gnina_script(  # noqa: C901, PLR0912
    ligands_dir: Path,
    cmd: list[str],
    output_sdf: Path,
    activate_cmd: str | None,
    ld_library_path: str | None,
    preparation_cmd: str | None,
    prepared_output_relative: str | None,
    protein_preparation_cmd: str | None,
) -> Path:
    """Create GNINA run script."""
    script_path = ligands_dir / "run_gnina.sh"
    with script_path.open("w") as f:
        f.write("#!/usr/bin/env bash\nset -eo pipefail\n")

        if activate_cmd:
            f.write(f"{activate_cmd}\n")
            if ld_library_path:
                f.write(
                    f'export LD_LIBRARY_PATH="{ld_library_path}:$LD_LIBRARY_PATH"\n'
                )

        if protein_preparation_cmd:
            parts = protein_preparation_cmd.split()
            protein_output_file = None
            protein_output_abs_path = None
            for part in reversed(parts):
                if ".pdb" in part:
                    if part.startswith("/"):
                        protein_output_file = part.split("/")[-1]
                        protein_output_abs_path = part
                    else:
                        protein_output_file = part
                        protein_output_abs_path = str(ligands_dir / part)
                    break
            f.write(f"cd {ligands_dir}\n")
            if protein_output_abs_path:
                f.write(
                    f'mkdir -p "$(dirname "{protein_output_abs_path}")"\n'
                    f'touch "{protein_output_abs_path}" 2>/dev/null || true\n'
                )
            f.write(
                'echo "Running protein preparation for GNINA..."\n'
                f"{protein_preparation_cmd} || PREP_EXIT_CODE=$?\n"
                'if [ ! -z "$PREP_EXIT_CODE" ]; then\n'
                '  echo "WARNING: Protein preparation command exited with code $PREP_EXIT_CODE"\n'
                "fi\n"
            )
            if protein_output_file:
                f.write(
                    "# If file exists in current dir but not at absolute path, move it\n"
                    f'if [ -f "{protein_output_file}" ] && [ ! -f "{protein_output_abs_path}" ]; then\n'
                    f'  mv "{protein_output_file}" "{protein_output_abs_path}"\n'
                    "fi\n"
                    "# Verify prepared file exists and is not empty\n"
                    f'if [ ! -f "{protein_output_abs_path}" ] && [ ! -f "{protein_output_file}" ]; then\n'
                    f'  echo "ERROR: Protein preparation failed - output file not found at {protein_output_abs_path}"\n'
                    '  echo "Current directory: $(pwd)"\n'
                    '  echo "Listing files in current directory:"\n'
                    "  ls -la\n"
                    "  exit 1\n"
                    f'elif [ -f "{protein_output_abs_path}" ] && [ ! -s "{protein_output_abs_path}" ]; then\n'
                    f'  echo "ERROR: Protein preparation failed - output file is empty: {protein_output_abs_path}"\n'
                    "  exit 1\n"
                    f'elif [ -f "{protein_output_file}" ] && [ ! -s "{protein_output_file}" ]; then\n'
                    f'  echo "ERROR: Protein preparation failed - output file is empty: {protein_output_file}"\n'
                    "  exit 1\n"
                    "else\n"
                    f'  echo "Protein preparation completed successfully: {protein_output_abs_path}"\n'
                    f'  ls -lh "{protein_output_abs_path}" 2>/dev/null || ls -lh "{protein_output_file}" 2>/dev/null || true\n'
                    "fi\n"
                )
            else:
                f.write(
                    "if [ $? -ne 0 ]; then\n"
                    '  echo "ERROR: Protein preparation failed"\n'
                    "  exit 1\n"
                    "fi\n"
                )

        if preparation_cmd and prepared_output_relative:
            f.write(
                f'mkdir -p "$(dirname "{prepared_output_relative}")"\n'
                f"{preparation_cmd}\n"
            )
            _write_file_wait_check(
                f,
                prepared_output_relative,
                f"ERROR: Preparation tool failed - output file {prepared_output_relative} not found after waiting",
                "ligand preparation",
            )
            f.write("rm -f ligands_raw.smi || true\n")

        receptor_path_in_cmd = None
        for i, arg in enumerate(cmd):
            if arg == "-r" and i + 1 < len(cmd):
                receptor_path_in_cmd = cmd[i + 1]
                break

        if receptor_path_in_cmd:
            f.write(
                'echo "Checking receptor file before GNINA docking..."\n'
                "# Final check - wait a moment and verify file still exists and is readable\n"
                "sleep 1\n"
                f'if [ ! -f "{receptor_path_in_cmd}" ]; then\n'
                f'  echo "ERROR: Receptor file not found: {receptor_path_in_cmd}"\n'
                '  echo "Current directory: $(pwd)"\n'
                '  echo "Listing files in current directory:"\n'
                "  ls -la\n"
                f'  if [ -d "$(dirname "{receptor_path_in_cmd}")" ]; then\n'
                '    echo "Listing files in receptor directory:"\n'
                f'    ls -la "$(dirname "{receptor_path_in_cmd}")"\n'
                "  fi\n"
                "  exit 1\n"
                f'elif [ ! -s "{receptor_path_in_cmd}" ]; then\n'
                f'  echo "ERROR: Receptor file is empty: {receptor_path_in_cmd}"\n'
                "  exit 1\n"
                f'elif [ ! -r "{receptor_path_in_cmd}" ]; then\n'
                f'  echo "ERROR: Receptor file is not readable: {receptor_path_in_cmd}"\n'
                "  exit 1\n"
                "else\n"
                f'  echo "Receptor file verified: {receptor_path_in_cmd}"\n'
                f'  ls -lh "{receptor_path_in_cmd}"\n'
                "fi\n"
            )

        f.write(
            f'mkdir -p "$(dirname "{output_sdf!s}")"\n'
            f': > "{output_sdf!s}"\n'
            'echo "Starting GNINA docking..."\n'
            f"{' '.join(cmd)}\n"
            "if [ $? -eq 0 ]; then\n"
            '  echo "GNINA docking completed successfully"\n'
            "else\n"
            '  echo "GNINA docking failed with exit code $?"\n'
            "  exit 1\n"
            "fi\n"
        )

    script_path.chmod(0o755)  # noqa: S103
    return script_path


def _generate_job_id(tool: str = "dock") -> str:
    """Generate a unique job ID."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")  # noqa: DTZ005
    unique_id = uuid.uuid4().hex[:8]
    return f"{tool}_{timestamp}_{unique_id}"


def _save_job_metadata(
    ligands_dir: Path,
    source_file: Path,
    num_ligands: int,
    receptor_pdb: str,
    tools_prepared: list[str],
    scripts_prepared: list[str],
    ligands_csv: Path,
    ligands_stats: dict[str, Any],
    job_ids: dict[str, str],
    overall_job_id: str,
) -> None:
    """Save job metadata to JSON file."""
    ligands_dir.mkdir(parents=True, exist_ok=True)
    metadata = {
        "job_id": overall_job_id,
        "timestamp": datetime.now().isoformat(timespec="seconds"),  # noqa: DTZ005
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
    with meta_path.open("w") as f:
        json.dump(metadata, f, indent=2)


def _save_job_ids(
    ligands_dir: Path,
    overall_job_id: str,
    job_ids: dict[str, str],
) -> None:
    """Save job IDs to a simple text file."""
    ids_path = ligands_dir / "job_ids.txt"
    try:
        with ids_path.open("w") as f:
            f.write(
                f"overall: {overall_job_id}\n"
                f"smina: {job_ids.get('smina', '')}\n"
                f"gnina: {job_ids.get('gnina', '')}\n"
            )
    except (OSError, IOError) as e:
        logger.warning("Failed to write job_ids.txt: %s", e)


def _update_metadata_with_run_status(
    ligands_dir: Path,
    run_status: dict[str, Any],
) -> None:
    """Update job metadata with run status."""
    meta_path = ligands_dir / "job_meta.json"
    try:
        metadata = {}
        if meta_path.exists():
            with meta_path.open() as f:
                metadata = json.load(f)
        metadata["run_status"] = run_status
        with meta_path.open("w") as f:
            json.dump(metadata, f, indent=2)
    except (OSError, IOError, json.JSONDecodeError) as e:
        logger.warning("Failed to update metadata with run status: %s", e)


def _run_docking_script(
    script_path: Path,
    working_dir: Path,
    log_path: Path,
    background: bool,  # noqa: FBT001
) -> dict[str, Any]:
    """Run a docking script and return execution status."""
    if not script_path.exists():
        logger.error("Script not found: %s", script_path)
        return {"status": "error", "log_path": str(log_path)}

    if background:
        with log_path.open("ab") as logf:
            subprocess.Popen(  # noqa: S603
                ["./" + script_path.name],
                stdout=logf,
                stderr=logf,
                cwd=str(working_dir),
            )
        logger.info("Started %s in background. Log: %s", script_path.name, log_path)
        return {"status": "started_background", "log_path": str(log_path)}
    with log_path.open("wb") as logf:
        subprocess.run(  # noqa: S603
            ["./" + script_path.name],
            check=True,
            stdout=logf,
            stderr=logf,
            cwd=str(working_dir),
        )
    logger.info("%s completed successfully. Log: %s", script_path.name, log_path)
    return {"status": "completed", "log_path": str(log_path)}


def _run_smina(
    ligands_dir: Path,
    background: bool,  # noqa: FBT001
    job_id: str,
) -> dict[str, Any]:
    """Run SMINA docking."""
    script_path = ligands_dir / "run_smina.sh"
    log_path = ligands_dir / "smina_run.log"
    status = _run_docking_script(script_path, ligands_dir, log_path, background)
    if job_id:
        status["job_id"] = job_id
    smina_out_dir = ligands_dir / "smina_results"
    status["results_dir"] = str(smina_out_dir) if smina_out_dir.exists() else None
    return status


def _run_gnina(
    ligands_dir: Path,
    output_sdf: Path,
    background: bool,  # noqa: FBT001
    job_id: str,
) -> dict[str, Any]:
    """Run GNINA docking."""
    script_path = ligands_dir / "run_gnina.sh"
    log_path = ligands_dir / "gnina_run.log"
    status = _run_docking_script(script_path, ligands_dir, log_path, background)
    if job_id:
        status["job_id"] = job_id
    status["output"] = str(output_sdf)
    status["log"] = str(log_path)
    return status


def _parse_tools_config(cfg: dict[str, Any]) -> list[str]:
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
    cfg: dict[str, Any],
    ligands_dir: Path,  # noqa: ARG001
    protein_preparation_tool: str | None,  # noqa: ARG001
    base_folder: Path | None = None,  # noqa: ARG001
) -> None:
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
            f"Receptor file not found: {original_receptor} (resolved to: {receptor_path})"
        )
        return

    cfg["receptor_pdb"] = str(receptor_path)


def _setup_smina(
    cfg: dict[str, Any],
    ligands_dir: Path,
    ligands_csv: Path,
    protein_preparation_tool: str | None,
    base_folder: Path | None = None,
) -> Path | None:
    """Set up SMINA docking configuration and script."""
    try:
        if protein_preparation_tool is None:
            _prepare_receptor_if_needed(cfg, ligands_dir, None, base_folder)

        receptor = cfg.get("receptor_pdb")
        if not receptor:
            logger.error("SMINA: receptor_pdb is missing in config")
            return None

        receptor_path = Path(receptor)
        if not receptor_path.exists():
            logger.error(f"SMINA: Receptor file not found: {receptor}")
            return None

        protein_prep_cmd = None
        if "protein_prepared.pdb" in str(receptor):
            logger.info(f"SMINA: Using prepared receptor: {receptor}")
        else:
            logger.info(f"SMINA: Using receptor: {receptor}")

        logger.debug(f"SMINA: Final protein_prep_cmd for script: {protein_prep_cmd}")

        ini_file = ligands_dir / "smina_auto.ini"
        _create_smina_config(cfg, ligands_csv, ini_file)
        activate_cmd = cfg.get("pyscreener_activate")
        script_path = _create_smina_script(
            ligands_dir, ini_file, activate_cmd, protein_prep_cmd
        )
    except (OSError, ValueError, RuntimeError) as e:
        logger.exception("Failed to setup SMINA: %s", e)
        return None
    else:
        logger.info("SMINA configuration prepared")
        return script_path


def _setup_gnina(  # noqa: C901, PLR0912, PLR0915
    cfg: dict[str, Any],
    base_folder: Path,
    ligands_dir: Path,
    ligands_csv: Path,
    ligand_preparation_tool: str | None,
    protein_preparation_tool: str | None,
) -> Path | None:
    """Set up GNINA docking configuration and script."""
    try:
        original_receptor = cfg.get("receptor_pdb")
        if not original_receptor:
            logger.error("GNINA: receptor_pdb is missing in config")
            return None

        if "protein_prepared.pdb" in original_receptor:
            logger.warning(
                f"GNINA: Config has prepared path: {original_receptor}, this should have been restored"
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

        receptor_path = Path(original_receptor)
        if not receptor_path.is_absolute():
            project_root = Path(__file__).parent.parent.parent.parent.parent
            receptor_path = (project_root / original_receptor).resolve()
            if not receptor_path.exists():
                receptor_path = Path(original_receptor).resolve()

        if not receptor_path.exists():
            logger.error(
                f"GNINA: receptor_pdb file not found: {original_receptor} (resolved to: {receptor_path})"
            )
            return None

        original_receptor = str(receptor_path)
        if protein_preparation_tool is None:
            receptor = cfg.get("receptor_pdb", original_receptor)
            protein_prep_cmd = None
            if "protein_prepared.pdb" in str(receptor):
                logger.info(f"GNINA: Using prepared receptor: {receptor}")
            else:
                logger.info(f"GNINA: Using receptor: {receptor}")
        else:
            receptor, protein_prep_cmd = _prepare_protein_for_docking(
                original_receptor, ligands_dir, protein_preparation_tool
            )
            if receptor != original_receptor:
                cfg["receptor_pdb"] = receptor
                logger.info(f"GNINA: Using prepared protein: {receptor}")
            else:
                cfg["receptor_pdb"] = original_receptor

        ligands_path, prep_cmd = _prepare_ligands_for_gnina(
            ligands_csv, ligands_dir, ligand_preparation_tool, cfg
        )
        gnina_dir = _get_gnina_output_directory(cfg, base_folder)
        gnina_dir.mkdir(parents=True, exist_ok=True)
        output_sdf = gnina_dir / "gnina_out.sdf"

        cmd = _build_gnina_command(cfg, receptor, ligands_path, output_sdf)
        activate_cmd, ld_library_path = _get_gnina_environment(cfg, base_folder)

        prepared_output_relative = None
        if prep_cmd and "-osd" in prep_cmd:
            parts = prep_cmd.split()
            idx = parts.index("-osd")
            if idx + 1 < len(parts):
                prepared_output_relative = parts[idx + 1]

        script_path = _create_gnina_script(
            ligands_dir,
            cmd,
            output_sdf,
            activate_cmd,
            ld_library_path,
            prep_cmd,
            prepared_output_relative,
            protein_prep_cmd,
        )
    except (OSError, ValueError, RuntimeError, FileNotFoundError) as e:
        logger.exception("Failed to setup GNINA: %s", e)
        return None
    else:
        logger.info("GNINA configuration prepared")
        return script_path


def run_docking(config: dict[str, Any]) -> bool:  # noqa: C901, PLR0911, PLR0912, PLR0915
    """Run docking orchestration function."""
    cfg = load_config(config["config_docking"])
    if not cfg.get("run", False):
        logger.info("Docking disabled in config")
        return False

    base_folder = Path(config["folder_to_save"]).resolve()
    source = _find_latest_input_source(base_folder)
    if source is None:
        logger.warning("No pass*SMILES.csv or sampledMols.csv found for docking input")
        return False

    try:
        df = pd.read_csv(source)
    except (OSError, IOError, pd.errors.ParserError) as e:
        logger.exception("Failed to read docking input %s: %s", source, e)
        return False

    ligands_dir = base_folder / "Docking"
    ligands_csv = ligands_dir / "ligands.csv"

    try:
        ligands_stats = _prepare_ligands_dataframe(df, ligands_csv)
    except ValueError as e:
        logger.error(f"Ligand preparation failed: {e}")
        return False

    ligand_preparation_tool = config.get("ligand_preparation_tool")
    protein_preparation_tool = config.get("protein_preparation_tool")
    tools_list = _parse_tools_config(cfg)
    logger.info(f"Docking tools configured: {tools_list}")

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
                        import subprocess  # noqa: PLC0415
                        import time  # noqa: PLC0415

                        result = subprocess.run(  # noqa: S602
                            prep_cmd,
                            check=False,
                            shell=True,
                            capture_output=True,
                            text=True,
                            cwd=str(ligands_dir),
                        )
                        if result.returncode != 0:
                            logger.error(
                                f"Protein preparation command failed: {result.stderr}"
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
                                    f"Protein prepared successfully: {prepared_receptor_path}"
                                )
                                break
                            time.sleep(wait_interval)
                            waited += wait_interval
                            if waited % 60 == 0:
                                logger.info(
                                    f"Waiting for protein preparation... ({waited}s)"
                                )

                        if not prepared_path.exists():
                            logger.error(
                                f"Protein preparation failed - output file not found after {max_wait}s: {prepared_receptor_path}"
                            )
                            logger.error(f"Command output: {result.stdout}")
                            logger.error(f"Command error: {result.stderr}")
                            return False
                        if prepared_path.stat().st_size == 0:
                            logger.error(
                                f"Protein preparation failed - output file is empty: {prepared_receptor_path}"
                            )
                            return False

                        cfg["receptor_pdb"] = prepared_receptor_path
                    except (
                        OSError,
                        IOError,
                        subprocess.CalledProcessError,
                        RuntimeError,
                    ) as e:
                        logger.exception("Failed to prepare protein: %s", e)
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
        except (OSError, IOError, ValueError, RuntimeError, FileNotFoundError) as e:
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
        logger.info(f"Docking job ID: {overall_job_id}")
    except (OSError, IOError, json.JSONEncodeError) as e:
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
            except (OSError, IOError, subprocess.CalledProcessError, RuntimeError) as e:
                logger.exception("GNINA execution failed: %s", e)
                run_status["gnina"] = {"status": "failed", "error": str(e)}

        try:
            _update_metadata_with_run_status(ligands_dir, run_status)
        except (OSError, IOError, json.JSONEncodeError) as e:
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
            if len(completed_tools) > 0:
                logger.warning(
                    f"Only {len(completed_tools)}/{len(selected_tools)} docking tools completed successfully"
                )
                return False
            logger.error("All docking tools failed")
            return False

    return True
