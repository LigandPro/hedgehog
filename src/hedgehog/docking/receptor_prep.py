import shlex
import subprocess
import time
from pathlib import Path

from hedgehog.configs.logger import logger
from hedgehog.docking.binaries import _validate_optional_tool_path


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

    project_root = Path(__file__).parent.parent.parent.parent

    if not receptor_path.is_absolute():
        candidate = (project_root / receptor_pdb).resolve()
        if candidate.exists():
            return candidate

    candidate = Path(receptor_pdb).resolve()
    if candidate.exists():
        return candidate

    return None


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


def _prepare_receptor_if_needed(
    cfg, ligands_dir, protein_preparation_tool, base_folder=None
):
    """Resolve receptor path and update config. Actual preparation happens in script."""
    original_receptor = cfg.get("receptor_pdb")
    if not original_receptor:
        return

    receptor_path = Path(original_receptor)
    if not receptor_path.is_absolute():
        project_root = Path(__file__).parent.parent.parent.parent
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


def _restore_gnina_receptor(cfg):
    """Validate and restore GNINA receptor_pdb config entry.

    Returns the receptor path string, or None if unrecoverable.
    """
    original_receptor = cfg.get("receptor_pdb")
    if not original_receptor:
        logger.error("GNINA: receptor_pdb is missing in config")
        return None

    if "protein_prepared.pdb" not in original_receptor:
        return original_receptor

    logger.warning(
        "GNINA: Config has prepared path: %s, this should have been restored",
        original_receptor,
    )
    project_root = Path(__file__).parent.parent.parent.parent
    possible_originals = [
        project_root / "src/hedgehog/configs/examples/7EW9_apo.pdb",
    ]
    for possible in possible_originals:
        if possible.exists():
            cfg["receptor_pdb"] = str(possible.resolve())
            return cfg["receptor_pdb"]

    logger.error("GNINA: Could not find original receptor file")
    return None


def _execute_protein_preparation(cfg, ligands_dir, protein_preparation_tool) -> bool:
    """Run protein preparation subprocess, wait for output, update cfg.

    Returns True on success (or graceful skip), False on hard failure.
    """
    _prepare_receptor_if_needed(cfg, ligands_dir, protein_preparation_tool)
    original_receptor = cfg.get("receptor_pdb")
    if not original_receptor:
        return True

    original_receptor_path = Path(original_receptor)
    if not original_receptor_path.exists():
        return True

    prepared_receptor_path, prep_cmd = _prepare_protein_for_docking(
        original_receptor, ligands_dir, protein_preparation_tool
    )
    if prepared_receptor_path == original_receptor or not prep_cmd:
        return True

    try:
        result = subprocess.run(
            prep_cmd,
            shell=False,
            capture_output=True,
            text=True,
            cwd=str(ligands_dir),
            timeout=600,
        )
    except FileNotFoundError:
        logger.warning(
            "Protein preparation tool is unavailable at runtime. "
            "Skipping receptor preprocessing and using original receptor: %s",
            original_receptor,
        )
        cfg["receptor_pdb"] = str(original_receptor_path)
        return True
    except Exception as e:
        logger.error("Failed to prepare protein: %s", e)
        import traceback

        logger.debug(traceback.format_exc())
        return False

    if result.returncode != 0:
        return _handle_prep_failure(
            result, cfg, original_receptor, original_receptor_path
        )

    if not _wait_for_prepared_file(prepared_receptor_path):
        _log_prep_wait_failure(prepared_receptor_path, result)
        return False

    logger.info("Protein prepared successfully: %s", prepared_receptor_path)
    cfg["receptor_pdb"] = prepared_receptor_path
    return True


def _handle_prep_failure(result, cfg, original_receptor, original_path) -> bool:
    """Handle non-zero return code from protein preparation.

    Returns True if preparation was skipped gracefully (tool unavailable),
    False if preparation hard-failed.
    """
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
        cfg["receptor_pdb"] = str(original_path)
        return True
    logger.error("Protein preparation command failed: %s", stderr or stdout)
    return False


def _log_prep_wait_failure(prepared_receptor_path, result):
    """Log details when prepared protein file fails to appear or is empty."""
    prepared_path = Path(prepared_receptor_path)
    if not prepared_path.exists():
        logger.error(
            "Protein preparation failed - output file not found after 300s: %s",
            prepared_receptor_path,
        )
        logger.error("Command output: %s", result.stdout)
        logger.error("Command error: %s", result.stderr)
    else:
        logger.error(
            "Protein preparation failed - output file is empty: %s",
            prepared_receptor_path,
        )


def _wait_for_prepared_file(path: str, max_wait: int = 300) -> bool:
    """Poll-wait for a prepared protein file to appear and be non-empty."""
    prepared_path = Path(path)
    wait_interval = 2
    waited = 0
    while waited < max_wait:
        if prepared_path.exists() and prepared_path.stat().st_size > 0:
            return True
        time.sleep(wait_interval)
        waited += wait_interval
        if waited % 60 == 0:
            logger.info("Waiting for protein preparation... (%ds)", waited)
    return False
