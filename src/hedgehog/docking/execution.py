import subprocess
import time
from pathlib import Path

from hedgehog._constants import TOOL_GNINA, TOOL_SMINA
from hedgehog.configs.logger import logger
from hedgehog.docking.aggregation import _aggregate_docking_results
from hedgehog.docking.metadata import _update_metadata_with_run_status
from hedgehog.docking.monitoring import _create_progress_tracker, _evaluate_run_results
from hedgehog.docking.paths import _emit_post_docking_warnings
from hedgehog.docking.scripts import _get_gnina_output_directory


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


def _execute_auto_run(cfg, tools_list, job_ids, ligands_dir, base_folder, reporter):
    """Execute docking tools and evaluate results. Returns True on success."""
    background = bool(cfg.get("run_in_background", False))
    run_status = {}
    selected_tools = [t for t in tools_list if t in job_ids]
    gnina_output_sdf = _get_gnina_output_directory(cfg, base_folder) / "gnina_out.sdf"

    progress_total = 0
    make_tick = None
    if reporter is not None and not background and selected_tools:
        progress_total, make_tick = _create_progress_tracker(
            reporter, selected_tools, ligands_dir, gnina_output_sdf
        )

    if TOOL_SMINA in job_ids:
        logger.info("Running SMINA docking")
        tick = make_tick(TOOL_SMINA) if make_tick and not background else None
        run_status[TOOL_SMINA] = _run_smina(
            ligands_dir, background, job_ids[TOOL_SMINA], tick=tick
        )
        if not background and run_status[TOOL_SMINA].get("status") == "completed":
            _emit_post_docking_warnings(
                TOOL_SMINA, ligands_dir / "_workdir" / "smina_run.log"
            )

    if TOOL_GNINA in job_ids:
        logger.info("Running GNINA docking")
        try:
            output_sdf = gnina_output_sdf
            tick = make_tick(TOOL_GNINA) if make_tick and not background else None
            run_status[TOOL_GNINA] = _run_gnina(
                ligands_dir, output_sdf, background, job_ids[TOOL_GNINA], tick=tick
            )
            if not background and run_status[TOOL_GNINA].get("status") == "completed":
                _emit_post_docking_warnings(
                    TOOL_GNINA,
                    ligands_dir / "_workdir" / "gnina_run.log",
                    output_sdf=output_sdf,
                )
        except Exception as e:
            logger.error("GNINA execution failed: %s", e)
            run_status[TOOL_GNINA] = {"status": "failed", "error": str(e)}

    try:
        _update_metadata_with_run_status(ligands_dir, run_status)
    except Exception as e:
        logger.warning("Failed to update metadata with run status: %s", e)

    if background:
        return True

    return _evaluate_run_results(selected_tools, run_status, reporter, progress_total)
