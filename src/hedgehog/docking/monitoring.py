from pathlib import Path

from hedgehog._constants import TOOL_GNINA, TOOL_SMINA
from hedgehog.configs.logger import logger


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


def _create_progress_tracker(reporter, selected_tools, ligands_dir, gnina_output_sdf):
    """Build progress tracking closures for docking execution.

    Returns (progress_total, make_tick_fn) where make_tick_fn(tool_name)
    returns a tick callable, or (0, None) if tracking is not needed.
    """
    configs_dir = ligands_dir / "_workdir" / "configs"
    tool_totals: dict[str, int] = {}
    for tool in selected_tools:
        count = 0
        if configs_dir.exists():
            count = len(list(configs_dir.glob(f"{tool}_*.ini")))
        tool_totals[tool] = count if count > 0 else 1
    progress_total = sum(tool_totals.values()) or 1

    def _count_done() -> int:
        done = 0
        if TOOL_SMINA in tool_totals:
            done += min(_count_smina_done(ligands_dir), tool_totals[TOOL_SMINA])
        if TOOL_GNINA in tool_totals:
            done += min(
                _count_gnina_done(ligands_dir, gnina_output_sdf),
                tool_totals[TOOL_GNINA],
            )
        return done

    def _make_tick(tool_name: str):
        def _tick() -> None:
            reporter.progress(
                _count_done(), progress_total, message=f"Docking ({tool_name})"
            )

        return _tick

    reporter.progress(0, progress_total, message="Docking")
    return progress_total, _make_tick


def _evaluate_run_results(selected_tools, run_status, reporter, progress_total) -> bool:
    """Evaluate docking results and return True if all tools succeeded."""
    completed_tools = [
        t for t in selected_tools if run_status.get(t, {}).get("status") == "completed"
    ]
    failed_tools = [
        t for t in selected_tools if run_status.get(t, {}).get("status") == "failed"
    ]

    if reporter is not None and progress_total:
        reporter.progress(progress_total, progress_total, message="Docking complete")

    if failed_tools:
        logger.error("Docking tools failed: %s", ", ".join(failed_tools))

    if len(completed_tools) == len(selected_tools):
        return True
    if completed_tools:
        logger.warning(
            "Only %d/%d docking tools completed successfully",
            len(completed_tools),
            len(selected_tools),
        )
    else:
        logger.error("All docking tools failed")
    return False
