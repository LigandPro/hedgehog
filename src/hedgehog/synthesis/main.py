from pathlib import Path

import pandas as pd

from hedgehog.configs.logger import load_config, logger
from hedgehog.synthesis.utils import (
    apply_synthesis_score_filters,
    calculate_synthesis_scores,
    get_input_path,
    merge_retrosynthesis_results,
    parse_retrosynthesis_results,
    prepare_input_smiles,
    run_aizynthfinder,
)
from hedgehog.utils.dataframe import (
    IDENTITY_COLUMNS,
    save_ordered_csv,
)
from hedgehog.utils.paths import process_path


def _get_aizynthfinder_config() -> Path:
    """Get default path to AiZynthFinder config file."""
    project_root = Path(__file__).resolve().parents[3]
    return (
        project_root
        / "modules"
        / "retrosynthesis"
        / "aizynthfinder"
        / "public"
        / "config.yml"
    )


def _get_aizynthfinder_root() -> Path:
    """Get path to AiZynthFinder project directory."""
    project_root = Path(__file__).resolve().parents[3]
    return project_root / "modules" / "retrosynthesis" / "aizynthfinder"


def _resolve_retrosynthesis_config(config: dict) -> Path:
    """Resolve retrosynthesis config path from runtime config or default."""
    custom_path = config.get("config_retrosynthesis")
    if custom_path:
        return Path(custom_path)
    return _get_aizynthfinder_config()


def _log_aizynthfinder_setup_instructions(config_path: Path) -> None:
    """Log instructions for setting up AiZynthFinder."""
    module_dir = config_path.parent.parent
    logger.error("AiZynthFinder config file not found: %s", config_path)
    logger.error("To set up retrosynthesis, run:")
    logger.error("  cd %s", module_dir)
    logger.error("  mkdir -p public")
    logger.error("  uv run python -m aizynthfinder.tools.download_public_data ./public")
    logger.error(
        "Synthesis stage will be skipped - continuing pipeline without retrosynthesis"
    )


def main(config: dict, reporter=None) -> None:
    """Main entry point for synthesis stage (retrosynthesis analysis).

    Args:
        config: Configuration dictionary containing pipeline settings
    """
    folder_to_save = process_path(config["folder_to_save"])
    output_folder = Path(folder_to_save) / "stages" / "04_synthesis"
    output_folder.mkdir(parents=True, exist_ok=True)
    config_synthesis = load_config(config["config_synthesis"])
    filtered_output = output_folder / "filtered_molecules.csv"

    input_path = get_input_path(config, folder_to_save)
    if not Path(input_path).exists():
        logger.error("Input file not found: %s", input_path)
        raise FileNotFoundError(f"Input file not found: {input_path}")

    input_df = pd.read_csv(input_path)
    logger.info("Processing %d filtered molecules", len(input_df))

    if len(input_df) == 0:
        logger.warning("No molecules to process for synthesis analysis")
        empty = pd.DataFrame(columns=IDENTITY_COLUMNS)
        save_ordered_csv(empty, filtered_output)
        logger.info("Saved 0 molecules to %s", filtered_output)
        return

    stage_total = 400

    def _progress_scores(phase: str, done: int, total: int) -> None:
        if reporter is None:
            return
        if total <= 0:
            pct = 0
        else:
            pct = int(round((done / total) * 100))
        pct = max(0, min(100, pct))
        base = 0 if phase == "sa_score" else 100 if phase == "syba_score" else 0
        reporter.progress(base + pct, stage_total, message=f"Computing {phase}")

    if reporter is not None:
        reporter.progress(0, stage_total, message="Computing synthesis scores")

    scored_df = calculate_synthesis_scores(
        input_df,
        folder_to_save,
        config_synthesis,
        progress_cb=_progress_scores if reporter is not None else None,
    )
    if reporter is not None:
        reporter.progress(200, stage_total, message="Applying score filters")
    save_ordered_csv(scored_df, output_folder / "synthesis_scores.csv")
    score_filtered_df = apply_synthesis_score_filters(scored_df, config_synthesis)

    if len(score_filtered_df) == 0:
        logger.warning("No molecules passed synthesis score filters")
        save_ordered_csv(score_filtered_df, filtered_output)
        logger.info("Saved 0 molecules to %s", filtered_output)
        if reporter is not None:
            reporter.progress(stage_total, stage_total, message="Synthesis complete")
        return

    # Check if retrosynthesis is enabled
    run_retrosynthesis = config_synthesis.get("run_retrosynthesis", True)
    if not run_retrosynthesis:
        logger.info("Retrosynthesis disabled (run_retrosynthesis=false)")
        save_ordered_csv(score_filtered_df, filtered_output)
        logger.info(
            "Saved %d molecules (scores only) to %s",
            len(score_filtered_df),
            filtered_output,
        )
        if reporter is not None:
            reporter.progress(stage_total, stage_total, message="Synthesis complete")
        return

    aizynthfinder_root = _get_aizynthfinder_root()
    aizynth_config = _resolve_retrosynthesis_config(config)
    if not aizynth_config.exists():
        from hedgehog.setup import ensure_aizynthfinder

        project_root = Path(__file__).resolve().parents[3]
        try:
            aizynth_config = ensure_aizynthfinder(project_root)
        except Exception:  # noqa: BLE001 — intentional: auto-setup can fail in many ways
            _log_aizynthfinder_setup_instructions(aizynth_config)
            save_ordered_csv(score_filtered_df, filtered_output)
            return

    prepare_input_smiles(score_filtered_df, output_folder / "input_smiles.smi")
    output_json = output_folder / "retrosynthesis_results.json"
    if reporter is not None:
        reporter.progress(
            300, stage_total, message="Running retrosynthesis (AiZynthFinder)"
        )
    if not run_aizynthfinder(
        output_folder / "input_smiles.smi",
        output_json,
        aizynth_config,
        aizynthfinder_dir=aizynthfinder_root,
    ):
        logger.error("Retrosynthesis analysis failed")
        raise RuntimeError("Retrosynthesis analysis failed")
    if reporter is not None:
        reporter.progress(stage_total, stage_total, message="Retrosynthesis complete")

    retrosynth_df = parse_retrosynthesis_results(output_json)
    if len(retrosynth_df) == 0:
        logger.warning("No retrosynthesis results found in JSON file")
        return

    merged_df = merge_retrosynthesis_results(score_filtered_df, retrosynth_df)
    save_ordered_csv(merged_df, output_folder / "synthesis_extended.csv")

    filter_solved_only = config_synthesis.get("filter_solved_only", True)
    if filter_solved_only:
        filtered_df = merged_df[merged_df["solved"] == 1].copy()
    else:
        filtered_df = merged_df.copy()
        logger.info("Keeping all molecules (filter_solved_only=False)")

    save_ordered_csv(filtered_df, filtered_output)

    if "search_time" in merged_df.columns:
        avg_time = merged_df["search_time"].mean()
        logger.info("Average retrosynthesis search time: %.2f s", avg_time)
