import os
import shlex
import subprocess
import time
from pathlib import Path

import pandas as pd

from hedgehog.configs.logger import logger
from hedgehog.docking.binaries import _validate_optional_tool_path
from hedgehog.utils.nvmolkit import maybe_enable_nvmolkit


def _project_root() -> Path:
    # src/hedgehog/docking/ligand_prep.py -> project root
    return Path(__file__).resolve().parents[3]


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

    # Try to get model_name and mol_idx if present
    has_model_name = "model_name" in df.columns
    has_mol_idx = "mol_idx" in df.columns
    model_name_series = df["model_name"] if has_model_name else None
    mol_idx_series = df["mol_idx"] if has_mol_idx else None

    maybe_enable_nvmolkit(
        project_root=_project_root(),
        context="docking.ligand_prep",
    )

    writer = Chem.SDWriter(str(sdf_path))
    written_count = 0

    for idx, (smi, name) in enumerate(
        zip(smiles_series.astype(str), name_series.astype(str), strict=False)
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

            # Preserve model_name and mol_idx as SDF properties
            if has_model_name and model_name_series is not None:
                model_name = str(model_name_series.iloc[idx])
                if model_name and model_name != "nan":
                    mol.SetProp("model_name", model_name)

            if has_mol_idx and mol_idx_series is not None:
                mol_idx = str(mol_idx_series.iloc[idx])
                if mol_idx and mol_idx != "nan":
                    mol.SetProp("mol_idx", mol_idx)

            writer.write(mol)
            written_count += 1
        except Exception:
            continue

    writer.close()

    if written_count == 0:
        raise RuntimeError("RDKit conversion produced 0 molecules for GNINA SDF")

    logger.info("Converted %d molecules to SDF using RDKit", written_count)
    return str(sdf_path.resolve()), None


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
