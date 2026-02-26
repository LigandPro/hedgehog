from pathlib import Path

import pandas as pd

from hedgehog._constants import CFG_DOCKING, KEY_FOLDER_TO_SAVE
from hedgehog.configs.logger import load_config, logger
from hedgehog.docking.aggregation import _aggregate_docking_results  # noqa: F401
from hedgehog.docking.binaries import _validate_optional_tool_path
from hedgehog.docking.config_writer import _create_docking_config_file  # noqa: F401
from hedgehog.docking.execution import _execute_auto_run
from hedgehog.docking.input import _find_latest_input_source, _prepare_ligands_dataframe
from hedgehog.docking.ligand_prep import _prepare_ligands_for_docking  # noqa: F401
from hedgehog.docking.metadata import (
    _generate_job_id,
    _parse_tools_config,
    _save_job_ids,
    _save_job_metadata,
)
from hedgehog.docking.monitoring import _create_progress_tracker  # noqa: F401
from hedgehog.docking.paths import _warn_if_autobox_far_from_receptor
from hedgehog.docking.receptor_prep import _execute_protein_preparation
from hedgehog.docking.scripts import _emit_manual_mode_warnings, _setup_docking_tools


def run(config, reporter=None):
    """Main docking orchestration function."""
    # 1. Load config and validate input
    cfg = load_config(config[CFG_DOCKING])
    if not cfg.get("run", False):
        logger.info("Docking disabled in config")
        return False

    base_folder = Path(config[KEY_FOLDER_TO_SAVE]).resolve()
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

    for tool in tools_list:
        _warn_if_autobox_far_from_receptor(cfg, tool)

    # 2. Protein preparation
    if protein_preparation_tool:
        if not _execute_protein_preparation(cfg, ligands_dir, protein_preparation_tool):
            return False

    # 3. Tool setup
    scripts_prepared, job_ids = _setup_docking_tools(
        cfg, tools_list, base_folder, ligands_dir, ligands_csv, ligand_preparation_tool
    )
    if not scripts_prepared:
        logger.error("No docking tools were successfully configured")
        return False

    overall_job_id = _generate_job_id("dock")
    try:
        _save_job_metadata(
            ligands_dir,
            source,
            len(df),
            cfg.get("receptor_pdb"),
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

    # 4. Execution
    if cfg.get("auto_run", True):
        return _execute_auto_run(
            cfg, tools_list, job_ids, ligands_dir, base_folder, reporter
        )

    # 5. Post-run warnings for manual mode
    _emit_manual_mode_warnings(cfg, tools_list, ligands_dir, base_folder)
    return True
