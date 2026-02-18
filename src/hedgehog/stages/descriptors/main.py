from pathlib import Path

from hedgehog.configs.logger import load_config, logger
from hedgehog.stages.descriptors.utils import (
    compute_metrics,
    draw_filtered_mols,
    filter_molecules,
    process_path,
)


def main(data, config, subfolder=None, reporter=None):
    """
    Computes default set of 22 physicochemical descriptors per molecule using RDKit,
    filters molecules based on configurable thresholds, and generates
    distribution plots.

    Args:
        data: DataFrame with molecules (must have 'smiles' column)
        config: Configuration file
        subfolder: Optional subfolder for output (e.g., 'stages/01_descriptors_initial' or 'stages/07_descriptors_final')
    """
    if data is None or len(data) == 0:
        logger.warning("No molecules provided for descriptor calculation. Skipping.")
        return None

    folder_to_save = Path(process_path(config["folder_to_save"]))
    subfolder = subfolder or str(Path("stages") / "01_descriptors_initial")
    descriptors_folder = folder_to_save / subfolder

    metrics_folder = descriptors_folder / "metrics"
    filtered_folder = descriptors_folder / "filtered"
    plots_folder = descriptors_folder / "plots"

    for folder in [descriptors_folder, metrics_folder, filtered_folder, plots_folder]:
        Path(folder).mkdir(parents=True, exist_ok=True)

    stage_total = 1000
    compute_base = 30
    compute_span = 720
    filter_start = compute_base + compute_span + 10
    filter_done = filter_start + 80
    plots_start = filter_done + 10
    plots_end = 990

    if reporter is not None:
        reporter.progress(10, stage_total, message="Loading descriptor config")

    config_descriptors = load_config(config["config_descriptors"])
    metrics_df = compute_metrics(
        data,
        metrics_folder,
        config=config,
        config_descriptors=config_descriptors,
        reporter=reporter,
        progress_stage_total=stage_total,
        progress_completed_base=compute_base,
        progress_completed_span=compute_span,
    )

    if config_descriptors["filter_data"]:
        if reporter is not None:
            reporter.progress(filter_start, stage_total, message="Applying descriptor filters")
        filter_molecules(metrics_df, config_descriptors["borders"], filtered_folder)
        if reporter is not None:
            reporter.progress(filter_done, stage_total, message="Saving descriptor outputs")

            def _plot_progress(done: int, total: int) -> None:
                if total <= 0:
                    mapped = plots_start
                else:
                    ratio = min(1.0, max(0.0, done / total))
                    mapped = plots_start + int(round((plots_end - plots_start) * ratio))
                reporter.progress(mapped, stage_total, message="Rendering descriptor plots")

        else:
            _plot_progress = None

        draw_filtered_mols(
            metrics_df,
            plots_folder,
            config,
            progress_cb=_plot_progress,
        )

    if reporter is not None:
        reporter.progress(stage_total, stage_total, message="Descriptors complete")

    return metrics_df
