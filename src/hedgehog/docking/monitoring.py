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


def _count_csv_rows(path: Path) -> int:
    """Count non-header rows in a CSV file; returns 0 if missing/unreadable."""
    try:
        if not path.exists():
            return 0
        with path.open("r", encoding="utf-8", errors="ignore") as handle:
            # header + data lines
            return max(sum(1 for _ in handle) - 1, 0)
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
    return 0


def _count_gnina_done(ligands_dir: Path, output_sdf: Path) -> int:
    status_dir = ligands_dir / "_workdir" / "gnina" / "status"
    if status_dir.exists():
        return _count_lines(status_dir / "success.txt") + _count_lines(
            status_dir / "failed.txt"
        )
    return 0


def _create_progress_tracker(reporter, selected_tools, ligands_dir, gnina_output_sdf):
    """Build per-tool molecule progress tracking closures for docking execution.

    Returns (tool_totals, make_tick_fn) where make_tick_fn(tool_name)
    returns a tick callable updating ``done/total`` in molecule units.
    """
    configs_dir = ligands_dir / "_workdir" / "configs"
    molecule_total = _count_csv_rows(ligands_dir / "ligands.csv")
    tool_totals: dict[str, int] = {}
    for tool in selected_tools:
        count = 0
        if configs_dir.exists():
            count = len(list(configs_dir.glob(f"{tool}_*.ini")))
        if count > 1:
            tool_totals[tool] = count
        else:
            tool_totals[tool] = molecule_total if molecule_total > 0 else 1

    def _count_done(tool_name: str) -> int:
        if tool_name == TOOL_SMINA:
            return min(_count_smina_done(ligands_dir), tool_totals.get(TOOL_SMINA, 1))
        if tool_name == TOOL_GNINA:
            return min(
                _count_gnina_done(ligands_dir, gnina_output_sdf),
                tool_totals.get(TOOL_GNINA, 1),
            )
        return 0

    def _make_tick(tool_name: str):
        total = tool_totals.get(tool_name, 1)

        def _tick() -> None:
            reporter.progress(
                _count_done(tool_name), total, message=f"Docking ({tool_name})"
            )

        return _tick

    return tool_totals, _make_tick


def _evaluate_run_results(selected_tools, run_status, reporter, tool_totals) -> bool:
    """Evaluate docking results and return True if all tools succeeded."""
    completed_tools = [
        t for t in selected_tools if run_status.get(t, {}).get("status") == "completed"
    ]
    failed_tools = [
        t for t in selected_tools if run_status.get(t, {}).get("status") == "failed"
    ]

    if reporter is not None and tool_totals:
        final_total = max(tool_totals.values()) if tool_totals else 1
        reporter.progress(final_total, final_total, message="Docking complete")

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
