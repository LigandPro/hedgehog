import json
import os
import shlex
import subprocess
import time
import uuid
from datetime import datetime
from pathlib import Path

import pandas as pd

from hedgehog.configs.logger import load_config, logger
from hedgehog.stages.docking.binary import (
    _auto_detect_cudnn_path,
    _get_gnina_environment,
    _get_gnina_output_directory,
    _is_real_binary,
    _join_existing_library_paths,
    _resolve_autobox_path,
    _resolve_docking_binary,
    _resolve_path,
    _resolve_receptor_path,
    _validate_optional_tool_path,
)
from hedgehog.stages.docking.config import (
    _create_docking_config_file,
    _create_gnina_config_file,
    _create_per_molecule_configs,
    _create_smina_config_file,
)
from hedgehog.stages.docking.scripts import (
    _build_gnina_command_template,
    _create_gnina_per_molecule_script,
    _create_gnina_script,
    _create_smina_per_molecule_script,
    _create_smina_script,
    _extract_pdb_output_from_cmd,
    _extract_prepared_output_from_cmd,
)
from hedgehog.utils.datamol_import import import_datamol_quietly
from hedgehog.utils.input_paths import find_latest_input_source

dm = import_datamol_quietly()

__all__ = [
    "_auto_detect_cudnn_path",
    "_build_gnina_command_template",
    "_create_docking_config_file",
    "_create_gnina_config_file",
    "_create_gnina_per_molecule_script",
    "_create_gnina_script",
    "_create_per_molecule_configs",
    "_create_smina_config_file",
    "_create_smina_per_molecule_script",
    "_create_smina_script",
    "_extract_pdb_output_from_cmd",
    "_extract_prepared_output_from_cmd",
    "_get_gnina_environment",
    "_get_gnina_output_directory",
    "_is_real_binary",
    "_join_existing_library_paths",
    "_resolve_autobox_path",
    "_resolve_docking_binary",
    "_resolve_path",
    "_resolve_receptor_path",
    "_validate_optional_tool_path",
    "run_docking",
]


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
        except (ValueError, TypeError):
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
            except (ValueError, IndexError):
                continue
            coords.append((x, y, z))
    except OSError:
        return []
    return coords


def _sdf_center(sdf_path: Path) -> tuple[float, float, float] | None:
    """Compute center of the first conformer in an SDF file (mean of atom positions)."""
    try:
        from rdkit import Chem
    except ImportError:
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
    except Exception:  # noqa: BLE001 — intentional: warning-only heuristic must never fail docking
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
    except ImportError:
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
            except (ValueError, TypeError):
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
    except Exception:  # noqa: BLE001 — intentional: warning-only logic must never fail docking
        return


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
    except (ValueError, TypeError):
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
    except (ValueError, TypeError):
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
                except (ValueError, TypeError):
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
        except (OSError, ValueError, RuntimeError) as e:
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

    logger.info("Converted %d molecules to SDF using RDKit", written_count)
    return str(sdf_path.resolve()), None


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
    except OSError as e:
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
    except (OSError, json.JSONDecodeError) as e:
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
        except (OSError, RuntimeError, ValueError) as e:
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
        except (OSError, RuntimeError, ValueError) as e:
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
    except (OSError, ValueError, RuntimeError, FileNotFoundError) as e:
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
        activate_cmd, ld_library_path = _get_gnina_environment(cfg)

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
    except (OSError, ValueError, RuntimeError, FileNotFoundError) as e:
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
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as e:
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
                    except (OSError, subprocess.SubprocessError) as e:
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
        except (OSError, ValueError, RuntimeError, FileNotFoundError) as e:
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
    except (OSError, ValueError) as e:
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
            except (OSError, subprocess.SubprocessError, RuntimeError) as e:
                logger.error("GNINA execution failed: %s", e)
                run_status["gnina"] = {"status": "failed", "error": str(e)}

        try:
            _update_metadata_with_run_status(ligands_dir, run_status)
        except (OSError, json.JSONDecodeError) as e:
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
        except Exception:  # noqa: BLE001 — intentional: post-docking warnings must never fail
            pass

    return True
