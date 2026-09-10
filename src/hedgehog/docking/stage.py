import json
from pathlib import Path

import pandas as pd

from hedgehog._constants import CFG_DOCKING, KEY_FOLDER_TO_SAVE
from hedgehog.alignment_runtime import alignment_stage_thresholds
from hedgehog.configs.logger import load_config, logger
from hedgehog.docking.aggregation import _collect_docking_stage_results
from hedgehog.docking.binaries import _validate_optional_tool_path
from hedgehog.docking.configuration import (
    DockingConfigError,
    normalize_docking_config,
    validate_docking_config,
)
from hedgehog.docking.execution import _execute_auto_run
from hedgehog.docking.input import _find_latest_input_source, _prepare_ligands_dataframe
from hedgehog.docking.metadata import (
    _generate_job_id,
    _parse_tools_config,
    _save_job_ids,
    _save_job_metadata,
)
from hedgehog.docking.paths import _warn_if_autobox_far_from_receptor
from hedgehog.docking.receptor_prep import _execute_protein_preparation
from hedgehog.docking.scripts import _emit_manual_mode_warnings, _setup_docking_tools

DOCKING_COMPLETED_EMPTY_MARKER = "completed_empty.marker"


def _mark_docking_completed_empty(
    ligands_dir: Path, source: Path, ligands_stats: dict, tools_list: list[str]
) -> Path:
    """Persist a marker for successful docking run with zero valid ligands."""
    ligands_dir.mkdir(parents=True, exist_ok=True)
    marker = ligands_dir / DOCKING_COMPLETED_EMPTY_MARKER
    payload = {
        "status": "completed_empty",
        "reason": "no_valid_ligands",
        "source_file": str(source),
        "tools": list(tools_list),
        "ligands_counts": dict(ligands_stats),
    }
    marker.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return marker


def run(config, reporter=None):
    """Main docking orchestration function."""
    # 1. Load config and validate input
    docking_config_path = Path(config[CFG_DOCKING]).expanduser().resolve()
    cfg = load_config(docking_config_path)
    if not cfg.get("run", False):
        logger.info("Docking disabled in config")
        return False

    try:
        tools_list = _parse_tools_config(cfg)
        cfg = normalize_docking_config(cfg, docking_config_path, tools_list)
        aligned_thresholds = alignment_stage_thresholds(
            config,
            CFG_DOCKING,
            "docking",
        )
        score_thresholds = (
            aligned_thresholds.get("score_thresholds")
            if isinstance(aligned_thresholds, dict)
            else None
        )
        if isinstance(score_thresholds, dict):
            cfg["score_thresholds"] = score_thresholds
    except (DockingConfigError, ValueError) as exc:
        logger.error("%s", exc)
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
    ligands_dir.mkdir(parents=True, exist_ok=True)
    df[["smiles", "model_name", "mol_idx"]].drop_duplicates(
        subset=["mol_idx"], keep="first"
    ).to_csv(ligands_dir / "input_molecules.csv", index=False)
    empty_marker = ligands_dir / DOCKING_COMPLETED_EMPTY_MARKER
    if empty_marker.exists():
        try:
            empty_marker.unlink()
        except OSError:
            pass

    try:
        ligands_stats = _prepare_ligands_dataframe(df, ligands_csv)
    except ValueError as e:
        logger.error("Ligand preparation failed: %s", e)
        return False

    ligand_preparation_tool = None
    if cfg.get("prepare_ligands", False):
        ligand_preparation_tool = _validate_optional_tool_path(
            config.get("ligand_preparation_tool"), "Ligand preparation tool"
        )
    protein_preparation_tool = _validate_optional_tool_path(
        config.get("protein_preparation_tool"), "Protein preparation tool"
    )
    logger.info("Docking tools configured: %s", tools_list)

    source_sdf = None
    cfg_sdf = config.get("docking_source_sdf")
    if cfg_sdf:
        candidate = Path(str(cfg_sdf)).expanduser()
        if candidate.exists():
            source_sdf = candidate.resolve()
    if source_sdf is None:
        candidate = base_folder / "input" / "ligands.sdf"
        if candidate.exists():
            source_sdf = candidate.resolve()
    if source_sdf is not None:
        for tool_name in tools_list:
            cfg[f"{tool_name}_ligands"] = str(source_sdf)
        logger.info("Using precomputed docking ligands SDF: %s", source_sdf)

    if int(ligands_stats.get("written", 0)) == 0:
        marker = _mark_docking_completed_empty(
            ligands_dir, source, ligands_stats, tools_list
        )
        logger.info(
            "No valid ligands for docking (%d/%d written). Marked run as completed-empty: %s",
            int(ligands_stats.get("written", 0)),
            int(ligands_stats.get("total", 0)),
            marker,
        )
        _collect_docking_stage_results(
            ligands_dir,
            tools_list,
            {},
            score_thresholds=cfg.get("score_thresholds"),
        )
        return True

    try:
        validate_docking_config(cfg, tools_list)
    except DockingConfigError as exc:
        logger.error("%s", exc)
        return False

    for tool in tools_list:
        _warn_if_autobox_far_from_receptor(cfg, tool)

    # 2. Protein preparation
    if protein_preparation_tool:
        if not _execute_protein_preparation(cfg, ligands_dir, protein_preparation_tool):
            return False

    # 3. Tool setup
    config_dir = docking_config_path.parent
    scripts_prepared, job_ids = _setup_docking_tools(
        cfg,
        tools_list,
        base_folder,
        ligands_dir,
        ligands_csv,
        ligand_preparation_tool,
        config_dir=config_dir,
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
