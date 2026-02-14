from pathlib import Path

import pandas as pd

from hedgehog.configs.logger import load_config, logger
from hedgehog.stages.synthesis.utils import *

IDENTITY_COLUMNS = ["smiles", "model_name", "mol_idx"]


def _order_identity_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Order dataframe columns with identity columns first."""
    existing_id_cols = [c for c in IDENTITY_COLUMNS if c in df.columns]
    ordered_cols = existing_id_cols + [
        c for c in df.columns if c not in IDENTITY_COLUMNS
    ]
    return df[ordered_cols]


def _save_ordered_csv(df: pd.DataFrame, path: Path) -> None:
    """Save DataFrame with identity columns ordered first."""
    path.parent.mkdir(parents=True, exist_ok=True)
    _order_identity_columns(df).to_csv(path, index=False)


def _get_aizynthfinder_config() -> Path:
    """Get path to AiZynthFinder config file."""
    project_root = Path(__file__).resolve().parents[4]
    return (
        project_root
        / "modules"
        / "retrosynthesis"
        / "aizynthfinder"
        / "public"
        / "config.yml"
    )


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
    _save_ordered_csv(scored_df, output_folder / "synthesis_scores.csv")
    score_filtered_df = apply_synthesis_score_filters(scored_df, config_synthesis)

    if len(score_filtered_df) == 0:
        logger.warning("No molecules passed synthesis score filters")
        _save_ordered_csv(score_filtered_df, filtered_output)
        logger.info("Saved 0 molecules to %s", filtered_output)
        if reporter is not None:
            reporter.progress(stage_total, stage_total, message="Synthesis complete")
        return

    # Check if retrosynthesis is enabled
    run_retrosynthesis = config_synthesis.get("run_retrosynthesis", True)
    if not run_retrosynthesis:
        logger.info("Retrosynthesis disabled (run_retrosynthesis=false)")
        _save_ordered_csv(score_filtered_df, filtered_output)
        logger.info(
            "Saved %d molecules (scores only) to %s",
            len(score_filtered_df),
            filtered_output,
        )
        if reporter is not None:
            reporter.progress(stage_total, stage_total, message="Synthesis complete")
        return

    aizynth_config = _get_aizynthfinder_config()
    if not aizynth_config.exists():
        from hedgehog.setup import ensure_aizynthfinder

        project_root = Path(__file__).resolve().parents[4]
        try:
            aizynth_config = ensure_aizynthfinder(project_root)
        except RuntimeError:
            _log_aizynthfinder_setup_instructions(aizynth_config)
            _save_ordered_csv(score_filtered_df, filtered_output)
            return

    prepare_input_smiles(score_filtered_df, output_folder / "input_smiles.smi")
    output_json = output_folder / "retrosynthesis_results.json"
    if reporter is not None:
        reporter.progress(
            300, stage_total, message="Running retrosynthesis (AiZynthFinder)"
        )
    if not run_aizynthfinder(
        output_folder / "input_smiles.smi", output_json, aizynth_config
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
    _save_ordered_csv(merged_df, output_folder / "synthesis_extended.csv")

    filter_solved_only = config_synthesis.get("filter_solved_only", True)
    if filter_solved_only:
        filtered_df = merged_df[merged_df["solved"] == 1].copy()
    else:
        filtered_df = merged_df.copy()
        logger.info("Keeping all molecules (filter_solved_only=False)")

    _save_ordered_csv(filtered_df, filtered_output)

    if "search_time" in merged_df.columns:
        avg_time = merged_df["search_time"].mean()
        logger.info("Average retrosynthesis search time: %.2f s", avg_time)
