"""Descriptors stage entry point."""

from pathlib import Path

import pandas as pd

from hedgehog._constants import CFG_DESCRIPTORS, KEY_FOLDER_TO_SAVE
from hedgehog.configs.logger import load_config, logger
from hedgehog.descriptors.compute import compute_metrics
from hedgehog.descriptors.filtering import filter_molecules
from hedgehog.descriptors.io import process_path
from hedgehog.descriptors.large import run_large
from hedgehog.descriptors.plotting import draw_filtered_mols
from hedgehog.descriptors.waves.registry import run_waves
from hedgehog.large_dataset import is_large_dataset_mode
from hedgehog.utils.mol_index import assign_mol_idx


def _prepare_identity_columns(data, run_base):
    """Ensure descriptor stage input always contains model_name and mol_idx."""
    prepared = data.copy()
    default_model = "single"

    if "model_name" not in prepared.columns:
        prepared["model_name"] = default_model
    else:
        model_series = prepared["model_name"].astype("string").str.strip()
        prepared["model_name"] = model_series.mask(model_series.eq(""), pd.NA).fillna(
            default_model
        )

    needs_mol_idx = "mol_idx" not in prepared.columns
    if not needs_mol_idx:
        mol_idx_series = prepared["mol_idx"].astype("string").str.strip()
        missing_mask = mol_idx_series.isna() | mol_idx_series.eq("")
        needs_mol_idx = bool(missing_mask.any())
        if needs_mol_idx:
            prepared["mol_idx"] = mol_idx_series.mask(missing_mask, pd.NA)

    if needs_mol_idx:
        assigned = assign_mol_idx(prepared, run_base=run_base, logger=logger)
        if "mol_idx" not in prepared.columns:
            prepared = assigned
        else:
            prepared["mol_idx"] = prepared["mol_idx"].fillna(assigned["mol_idx"])

    return prepared


def run(data, config, subfolder=None, reporter=None):
    """Compute physicochemical descriptors, filter, and plot distributions.

    Computes default set of 22 physicochemical descriptors per molecule using RDKit,
    filters molecules based on configurable thresholds, and generates
    distribution plots.

    Args:
        data: DataFrame with molecules (must have 'smiles' column)
        config: Configuration file
        subfolder: Optional subfolder for output (e.g., 'stages/02_descriptors_initial'
                   or 'stages/07_descriptors_final')
        reporter: Optional progress reporter instance
    """
    folder_to_save = Path(process_path(config[KEY_FOLDER_TO_SAVE]))
    subfolder = subfolder or str(Path("stages") / "02_descriptors_initial")
    descriptors_folder = folder_to_save / subfolder

    metrics_folder = descriptors_folder / "metrics"
    filtered_folder = descriptors_folder / "filtered"
    plots_folder = descriptors_folder / "plots"

    for folder in [descriptors_folder, metrics_folder, filtered_folder, plots_folder]:
        Path(folder).mkdir(parents=True, exist_ok=True)

    molecule_total = max(len(data), 1) if data is not None else 1

    if reporter is not None:
        reporter.progress(0, molecule_total, message="Loading descriptor config")

    config_descriptors = load_config(config[CFG_DESCRIPTORS])
    if is_large_dataset_mode(config):
        return run_large(
            data,
            config,
            config_descriptors,
            descriptors_folder,
            reporter=reporter,
        )

    if data is None or len(data) == 0:
        logger.warning("No molecules provided for descriptor calculation. Skipping.")
        return None

    data = _prepare_identity_columns(data, run_base=folder_to_save)
    molecule_total = max(len(data), 1)
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
        filter_molecules(
            metrics_df,
            config_descriptors["borders"],
            filtered_folder,
            structural_constraints=config_descriptors.get("structural_constraints"),
        )
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
