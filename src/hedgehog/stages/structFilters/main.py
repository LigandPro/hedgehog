from pathlib import Path

import pandas as pd

from hedgehog.configs.logger import load_config, logger
from hedgehog.stages.structFilters.utils import (
    filter_data,
    filter_function_applier,
    get_basic_stats,
    inject_identity_columns_to_all_csvs,
    plot_calculated_stats,
    plot_filter_failures_analysis,
    plot_restriction_ratios,
    process_one_dataframe,
    process_one_file,
    process_path,
)
from hedgehog.utils.input_paths import find_sampled_molecules

IDENTITY_COLUMNS = ["smiles", "model_name", "mol_idx"]


def _is_post_descriptors_stage(stage_dir):
    """Check if stage directory indicates post-descriptors processing."""
    return "03_structural_filters_post" in stage_dir or stage_dir == "StructFilters"


def _order_identity_columns(df):
    """Order dataframe columns with identity columns first."""
    existing_id_cols = [c for c in IDENTITY_COLUMNS if c in df.columns]
    ordered_cols = existing_id_cols + [
        c for c in df.columns if c not in IDENTITY_COLUMNS
    ]
    return df[ordered_cols]


def _find_existing_file(paths):
    """Return first existing non-empty file path from list, or None."""
    for path in paths:
        p = Path(path)
        if p.exists() and p.stat().st_size > 0:
            return path
    return None


def _get_sampled_molecules_path(folder_to_save):
    """Get path to sampled molecules file, checking both new and legacy locations."""
    path = find_sampled_molecules(Path(folder_to_save))
    return str(path) if path else None


def _get_input_path(config, stage_dir, folder_to_save):
    """Determine input path based on stage directory.

    Args:
        config: Configuration dictionary
        stage_dir: Stage directory path (e.g., 'stages/02_structural_filters_pre' or 'stages/03_structural_filters_post')
        folder_to_save: Base output folder

    Returns:
        Path to input file
    """
    descriptors_enabled = False
    if "config_descriptors" in config:
        try:
            desc_config = load_config(config["config_descriptors"])
            descriptors_enabled = desc_config.get("run", False)
        except Exception:
            pass

    is_pre_descriptors = (
        "02_structural_filters_pre" in stage_dir or "pre" in stage_dir.lower()
    )
    is_post_descriptors = _is_post_descriptors_stage(stage_dir)

    if is_post_descriptors:
        base = Path(folder_to_save)
        descriptors_candidates = [
            str(
                base
                / "stages"
                / "01_descriptors_initial"
                / "filtered"
                / "filtered_molecules.csv"
            ),
            str(base / "Descriptors" / "passDescriptorsSMILES.csv"),
        ]
        descriptors_path = _find_existing_file(descriptors_candidates)
        if descriptors_path:
            return descriptors_path

        if descriptors_enabled:
            logger.error(
                "Descriptors are enabled but output not found. Cannot proceed with post-descriptors filters."
            )
            raise FileNotFoundError(
                f"Descriptors output not found at {descriptors_candidates}"
            )

        pre_filters_path = (
            base / "stages" / "02_structural_filters_pre" / "filtered_molecules.csv"
        )
        if pre_filters_path.exists() and pre_filters_path.stat().st_size > 0:
            return str(pre_filters_path)

        sampled_path = _get_sampled_molecules_path(folder_to_save)
        if sampled_path:
            return sampled_path

    if is_pre_descriptors:
        sampled_path = _get_sampled_molecules_path(folder_to_save)
        if sampled_path:
            return sampled_path

    return config["generated_mols_path"]


def _load_input_data(input_path):
    """Load and prepare input data with required columns.

    Args:
        input_path: Path to input CSV file

    Returns:
        Tuple of (input_df, model_name) where model_name is a single name or list

    Raises:
        Exception: If input data cannot be loaded
    """
    input_df = pd.read_csv(input_path)

    if "model_name" not in input_df.columns:
        inferred_model = Path(input_path).stem
        input_df["model_name"] = inferred_model
        logger.info("model_name column missing; using '%s'", inferred_model)

    if "mol_idx" not in input_df.columns:
        input_df["mol_idx"] = range(len(input_df))

    model_names = sorted(input_df["model_name"].dropna().unique().tolist())
    model_name = model_names[0] if len(model_names) == 1 else model_names

    return input_df, model_name


def _ensure_pass_column(df, filter_name):
    """Ensure DataFrame has a 'pass' column, creating from alternatives if needed."""
    if "pass" in df.columns:
        return df

    if "pass_filter" in df.columns:
        df["pass"] = df["pass_filter"]
    else:
        logger.warning(
            "'pass' column not found in %s extended results. Assuming all molecules pass.",
            filter_name,
        )
        df["pass"] = True

    return df


def _save_filter_results(output_dir, filter_name, metrics_df, extended_df):
    """Save filter results to subdirectory with consistent formatting.

    Args:
        output_dir: Base output directory
        filter_name: Name of the filter
        metrics_df: DataFrame with filter metrics
        extended_df: DataFrame with extended results
    """
    filter_subdir = Path(output_dir) / filter_name
    filter_subdir.mkdir(parents=True, exist_ok=True)

    metrics_df = _order_identity_columns(metrics_df)
    metrics_df.to_csv(filter_subdir / "metrics.csv", index=False)

    extended_df = _order_identity_columns(extended_df)
    extended_df.to_csv(filter_subdir / "extended.csv", index=False)

    extended_df = _ensure_pass_column(extended_df, filter_name)
    filtered_mols = extended_df[extended_df["pass"]].copy()
    filtered_mols = _order_identity_columns(filtered_mols)
    filtered_mols.to_csv(filter_subdir / "filtered_molecules.csv", index=False)


def _get_enabled_filters(config_structFilters):
    """Extract enabled filters from configuration."""
    return {
        k.replace("calculate_", ""): v
        for k, v in config_structFilters.items()
        if "calculate_" in k and v
    }


def main(config, stage_dir, reporter=None):
    """Main entry point for structural filters stage.

    Args:
        config: Configuration dictionary
        stage_dir: Stage directory path (e.g., 'stages/02_structural_filters_pre' or 'stages/03_structural_filters_post')
    """
    sample_size = config["sample_size"]
    folder_to_save = Path(process_path(config["folder_to_save"]))
    output_dir = folder_to_save / stage_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    input_path = _get_input_path(config, stage_dir, folder_to_save)

    try:
        input_df, model_name = _load_input_data(input_path)
    except Exception as e:
        logger.error("Could not load input data from %s: %s", input_path, e)
        raise

    config_structFilters = load_config(config["config_structFilters"])
    filters_to_calculate = _get_enabled_filters(config_structFilters)

    # Determine if input is CSV (use DataFrame directly to avoid re-reading)
    is_csv = input_path.lower().endswith(".csv")

    filter_names = list(filters_to_calculate)
    stage_total = max(1, len(filter_names) * 100)

    for idx, filter_name in enumerate(filter_names):
        base = idx * 100
        if reporter is not None:
            reporter.progress(
                base, stage_total, message=f"StructFilters: {filter_name}"
            )

        apply_func = filter_function_applier(filter_name)
        progress_cb = None
        if reporter is not None and filter_name == "common_alerts":

            def _alerts_progress(
                done: int, total: int, base_value: int = base, name: str = filter_name
            ) -> None:
                if total <= 0:
                    pct = 0
                else:
                    pct = int(round((done / total) * 100))
                pct = max(0, min(100, pct))
                reporter.progress(
                    base_value + pct,
                    stage_total,
                    message=f"StructFilters: {name}",
                )

            progress_cb = _alerts_progress

        if is_csv:
            filter_results = process_one_dataframe(
                config, input_df, apply_func, sample_size, progress_cb=progress_cb
            )
        else:
            filter_results = process_one_file(
                config, input_path, apply_func, sample_size, progress_cb=progress_cb
            )

        if filter_results is None:
            logger.warning("No molecules to process for filter: %s", filter_name)
            continue

        final_res, final_extended = get_basic_stats(
            config_structFilters, filter_results, model_name, filter_name=filter_name
        )
        _save_filter_results(output_dir, filter_name, final_res, final_extended)
        if reporter is not None:
            reporter.progress(
                base + 100, stage_total, message=f"StructFilters: {filter_name}"
            )

    plot_calculated_stats(config, stage_dir)
    plot_restriction_ratios(config, stage_dir)

    if _is_post_descriptors_stage(stage_dir):
        plot_filter_failures_analysis(config, stage_dir)

    is_single_stage = config.get("_run_single_stage_override") == "struct_filters"
    if config_structFilters.get("filter_data", False) or is_single_stage:
        filter_data(config, stage_dir)

    inject_identity_columns_to_all_csvs(config, stage_dir)

    if reporter is not None:
        reporter.progress(stage_total, stage_total, message="StructFilters complete")
