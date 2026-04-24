from pathlib import Path

import pandas as pd

from hedgehog._constants import KEY_FOLDER_TO_SAVE
from hedgehog.configs.logger import load_config, logger
from hedgehog.synthesis.utils import *

IDENTITY_COLUMNS = ["smiles", "model_name", "mol_idx"]


def _project_root() -> Path:
    """Return project root path for this module."""
    return Path(__file__).resolve().parents[3]


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
    """Get default path to AiZynthFinder config file."""
    return _project_root() / "modules" / "aizynthfinder" / "public" / "config.yml"


def _resolve_retrosynthesis_config(config: dict) -> Path:
    """Resolve retrosynthesis config path from runtime config or default."""
    custom_path = config.get("config_retrosynthesis")
    if custom_path:
        resolved = Path(custom_path).expanduser()
        if not resolved.is_absolute():
            resolved = _project_root() / resolved
        return resolved.resolve(strict=False)
    return _get_aizynthfinder_config()


def _log_aizynthfinder_setup_instructions(config_path: Path) -> None:
    """Log instructions for setting up AiZynthFinder."""
    logger.error("AiZynthFinder config file not found: %s", config_path)
    logger.error("To set up retrosynthesis, run:")
    logger.error("  uv sync --extra retrosynthesis")
    logger.error("  uv run hedgehog setup aizynthfinder")
    logger.error(
        "Synthesis stage will be skipped - continuing pipeline without retrosynthesis"
    )


def main(config: dict, reporter=None) -> None:
    """Main entry point for synthesis stage (retrosynthesis analysis).

    Args:
        config: Configuration dictionary containing pipeline settings
    """
    folder_to_save = process_path(config[KEY_FOLDER_TO_SAVE])
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
        _save_ordered_csv(empty, filtered_output)
        logger.info("Saved 0 molecules to %s", filtered_output)
        return

    total_input_mols = max(len(input_df), 1)

    def _report_progress(done: int, total: int, message: str) -> None:
        if reporter is not None:
            reporter.progress(done, total, message=message)

    def _progress_scores(phase: str, done: int, total: int) -> None:
        phase_total = total if total > 0 else total_input_mols
        phase_done = max(0, min(done, phase_total))
        _report_progress(phase_done, phase_total, f"Computing {phase}")

    _report_progress(0, total_input_mols, "Computing synthesis scores")

    scored_df = calculate_synthesis_scores(
        input_df,
        folder_to_save,
        config_synthesis,
        progress_cb=_progress_scores if reporter is not None else None,
    )
    _report_progress(0, total_input_mols, "Applying score filters")
    _save_ordered_csv(scored_df, output_folder / "synthesis_scores.csv")
    score_filtered_df = apply_synthesis_score_filters(scored_df, config_synthesis)
    _report_progress(total_input_mols, total_input_mols, "Applying score filters")

    if len(score_filtered_df) == 0:
        logger.warning("No molecules passed synthesis score filters")
        _save_ordered_csv(score_filtered_df, filtered_output)
        logger.info("Saved 0 molecules to %s", filtered_output)
        _report_progress(total_input_mols, total_input_mols, "Synthesis complete")
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
        _report_progress(total_input_mols, total_input_mols, "Synthesis complete")
        return

    aizynth_config = _resolve_retrosynthesis_config(config)
    if not aizynth_config.exists():
        from hedgehog.setup import ensure_aizynthfinder

        def _handle_setup_failure() -> None:
            _log_aizynthfinder_setup_instructions(aizynth_config)
            _save_ordered_csv(score_filtered_df, filtered_output)

        project_root = _project_root()
        try:
            aizynth_config = ensure_aizynthfinder(project_root)
        except Exception as exc:
            if isinstance(exc, RuntimeError) and "supports Python 3.10-3.12" in str(
                exc
            ):
                raise
            _handle_setup_failure()
            return

    input_smiles_file = output_folder / "input_smiles.smi"
    prepare_input_smiles(score_filtered_df, input_smiles_file)
    output_json = output_folder / "retrosynthesis_results.json"

    retrosynthesis_total = max(len(score_filtered_df), 1)

    def _progress_retrosynthesis(
        done: int, total: int, detail: str | None = None
    ) -> None:
        phase_total = total if total > 0 else retrosynthesis_total
        phase_done = max(0, min(done, phase_total))
        message = "Running retrosynthesis (AiZynthFinder)"
        if detail:
            message = f"{message}: {detail}"
        _report_progress(phase_done, phase_total, message)

    _progress_retrosynthesis(0, retrosynthesis_total, "starting")

    if not run_aizynthfinder(
        input_smiles_file,
        output_json,
        aizynth_config,
        synthesis_config=config_synthesis,
        progress_cb=_progress_retrosynthesis if reporter is not None else None,
    ):
        logger.error("Retrosynthesis analysis failed")
        raise RuntimeError("Retrosynthesis analysis failed")
    _report_progress(
        retrosynthesis_total,
        retrosynthesis_total,
        "Retrosynthesis complete",
    )

    retrosynth_df = parse_retrosynthesis_results(output_json)
    if len(retrosynth_df) == 0:
        logger.warning("No retrosynthesis results found in JSON file")
        empty = score_filtered_df.iloc[0:0].copy()
        _save_ordered_csv(empty, output_folder / "synthesis_extended.csv")
        _save_ordered_csv(empty, filtered_output)
        logger.info("Saved 0 molecules to %s", filtered_output)
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
