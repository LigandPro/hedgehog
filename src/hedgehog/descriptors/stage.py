"""Descriptors stage entry point."""

from pathlib import Path

from hedgehog._constants import CFG_DESCRIPTORS, KEY_FOLDER_TO_SAVE
from hedgehog.configs.logger import load_config, logger
from hedgehog.descriptors.compute import compute_metrics
from hedgehog.descriptors.filtering import filter_molecules
from hedgehog.descriptors.io import process_path
from hedgehog.descriptors.plotting import draw_filtered_mols
from hedgehog.descriptors.waves.registry import run_waves


def run(data, config, subfolder=None, reporter=None):
    """Compute physicochemical descriptors, filter, and plot distributions.

    Computes default set of 22 physicochemical descriptors per molecule using RDKit,
    filters molecules based on configurable thresholds, and generates
    distribution plots.

    Args:
        data: DataFrame with molecules (must have 'smiles' column)
        config: Configuration file
        subfolder: Optional subfolder for output (e.g., 'stages/01_descriptors_initial'
                   or 'stages/07_descriptors_final')
        reporter: Optional progress reporter instance
    """
    if data is None or len(data) == 0:
        logger.warning("No molecules provided for descriptor calculation. Skipping.")
        return None

    folder_to_save = Path(process_path(config[KEY_FOLDER_TO_SAVE]))
    subfolder = subfolder or str(Path("stages") / "01_descriptors_initial")
    descriptors_folder = folder_to_save / subfolder

    metrics_folder = descriptors_folder / "metrics"
    filtered_folder = descriptors_folder / "filtered"
    plots_folder = descriptors_folder / "plots"

    for folder in [descriptors_folder, metrics_folder, filtered_folder, plots_folder]:
        Path(folder).mkdir(parents=True, exist_ok=True)

    molecule_total = max(len(data), 1)

    if reporter is not None:
        reporter.progress(0, molecule_total, message="Loading descriptor config")

    config_descriptors = load_config(config[CFG_DESCRIPTORS])
    metrics_df = compute_metrics(
        data,
        metrics_folder,
        config=config,
        config_descriptors=config_descriptors,
        reporter=reporter,
    )

    if config_descriptors["filter_data"]:
        if reporter is not None:
            reporter.progress(0, molecule_total, message="Applying descriptor filters")
        filter_molecules(metrics_df, config_descriptors, filtered_folder)
        if reporter is not None:
            reporter.progress(
                molecule_total, molecule_total, message="Saving descriptor outputs"
            )

            def _plot_progress(done: int, total: int) -> None:
                if total <= 0:
                    mapped = 0
                else:
                    ratio = min(1.0, max(0.0, done / total))
                    mapped = int(round(ratio * molecule_total))
                reporter.progress(
                    mapped, molecule_total, message="Rendering descriptor plots"
                )

        else:
            _plot_progress = None

        draw_filtered_mols(
            metrics_df,
            plots_folder,
            config,
            progress_cb=_plot_progress,
        )

    if reporter is not None:
        reporter.progress(
            molecule_total, molecule_total, message="Descriptors complete"
        )

    # Run wave extensions
    context = {"metrics_df": metrics_df, "config": config, "subfolder": subfolder}
    run_waves(context)

    return metrics_df
