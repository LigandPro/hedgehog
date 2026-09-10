from pathlib import Path
from time import perf_counter

import pandas as pd

from hedgehog._constants import CFG_DESCRIPTORS, CFG_STRUCT_FILTERS, KEY_FOLDER_TO_SAVE
from hedgehog.configs.logger import load_config, logger
from hedgehog.large_dataset import is_large_dataset_mode
from hedgehog.struct_filters.common_alert_diagnostics import (
    write_common_alert_diagnostics,
)
from hedgehog.struct_filters.large import run_large
from hedgehog.struct_filters.utils import (
    attach_structural_enforcement_pass,
    build_structural_policy_pass_masks,
    combine_filter_results_in_memory,
    filter_function_applier,
    finalize_structural_liability_profile,
    get_basic_stats,
    initialize_structural_liability_profile,
    inject_identity_columns_to_all_csvs,
    merge_structural_liability_profile,
    merge_structural_policy_aliases,
    plot_calculated_stats,
    plot_filter_failures_analysis,
    plot_restriction_ratios,
    prepare_structfilters_input,
    process_one_dataframe,
    process_one_file,
    process_path,
    process_prepared_payload,
)
from hedgehog.struct_filters.waves.registry import (
    get_aligned_enforced_filters,
    get_calculated_policy_names,
    policy_calculation_filter,
)
from hedgehog.utils.input_paths import find_sampled_molecules
from hedgehog.utils.mol_index import assign_mol_idx
from hedgehog.utils.parallel import resolve_n_jobs

IDENTITY_COLUMNS = ["smiles", "model_name", "mol_idx"]

_FILTER_DESCRIPTIONS: dict[str, str] = {
    "common_alerts": "SMARTS-based structural alert screening using curated rule sets (PAINS, Dundee, BMS, Glaxo, etc.).",
    "molgraph_stats": "Generator graph-sanity patterns with catalog severities 5, 8, and 10.",
    "molcomplexity": "Molecular complexity heuristics (medchem complexity filters).",
    "NIBR": "Novartis in-house structural filters (severity-based developability heuristics).",
    "bredt": "Bredt's rule violations at bridgehead positions in small bicyclic systems.",
    "lilly": "Eli Lilly Medchem Rules demerit scoring; higher demerits indicate less desirable structures.",
    "protecting_groups": "Curated terminal Fmoc, N-tert-butoxymethyl, and Boc-family motifs that are undesirable in final compounds.",
    "ring_infraction": "Strained/unusual ring systems and ring infractions (configurable heterocycle minimum size).",
    "stereo_center": "Total stereocenter count diagnostic with an optional independent hard gate.",
    "undefined_stereo_center": "Undefined stereocenters above the configured inclusive maximum.",
    "halogenicity": "Excessive halogen content (configurable thresholds for F/Cl/Br).",
    "symmetry": "Highly symmetric molecules (optional; configurable symmetry threshold).",
}


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
        stage_dir: Stage directory path (e.g., 'stages/03_structural_filters_post')
        folder_to_save: Base output folder

    Returns:
        Path to input file
    """
    descriptors_enabled = False
    if CFG_DESCRIPTORS in config:
        try:
            desc_config = load_config(config[CFG_DESCRIPTORS])
            descriptors_enabled = desc_config.get("run", False)
        except Exception:
            pass

    is_post_descriptors = _is_post_descriptors_stage(stage_dir)
    is_single_struct_filters = (
        config.get("_run_single_stage_override") == "struct_filters"
    )

    if is_post_descriptors:
        base = Path(folder_to_save)
        descriptors_candidates = [
            str(base / "stages" / "02_descriptors_initial" / "filtered_molecules.csv"),
            str(
                base
                / "stages"
                / "02_descriptors_initial"
                / "filtered"
                / "filtered_molecules.csv"
            ),
            str(base / "Descriptors" / "passDescriptorsSMILES.csv"),
        ]
        descriptors_path = _find_existing_file(descriptors_candidates)
        if descriptors_path:
            return descriptors_path

        if descriptors_enabled and not is_single_struct_filters:
            logger.error(
                "Descriptors are enabled but output not found. Cannot proceed with post-descriptors filters."
            )
            raise FileNotFoundError(
                f"Descriptors output not found at {descriptors_candidates}"
            )
        if descriptors_enabled and is_single_struct_filters:
            logger.warning(
                "Single-stage struct_filters: descriptors output not found; falling back to mol_prep/sampled input."
            )

        mol_prep_path = base / "stages" / "01_mol_prep" / "filtered_molecules.csv"
        if mol_prep_path.exists() and mol_prep_path.stat().st_size > 0:
            return str(mol_prep_path)

        sampled_path = _get_sampled_molecules_path(folder_to_save)
        if sampled_path:
            return sampled_path

    return config["generated_mols_path"]


def _load_input_data(input_path, run_base=None):
    """Load and prepare input data with required columns.

    Args:
        input_path: Path to input CSV file
        run_base: Base path for assign_mol_idx mapping persistence

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
    else:
        model_series = input_df["model_name"].astype("string").str.strip()
        missing_models = model_series.isna() | model_series.eq("")
        if missing_models.any():
            inferred_model = Path(input_path).stem
            input_df["model_name"] = model_series.mask(
                missing_models, inferred_model
            ).astype(object)
            logger.info(
                "model_name had %d empty value(s); using '%s' for those rows",
                int(missing_models.sum()),
                inferred_model,
            )

    run_base_path = Path(run_base) if run_base is not None else Path(input_path).parent
    needs_assignment = "mol_idx" not in input_df.columns
    if not needs_assignment:
        mol_idx_series = input_df["mol_idx"].astype("string").str.strip()
        missing_mask = mol_idx_series.isna() | mol_idx_series.eq("")
        needs_assignment = bool(missing_mask.any())
        if needs_assignment:
            input_df["mol_idx"] = mol_idx_series.mask(missing_mask, pd.NA)
    if needs_assignment:
        assigned = assign_mol_idx(input_df, run_base=run_base_path, logger=logger)
        if "mol_idx" not in input_df.columns:
            input_df = assigned
        else:
            input_df["mol_idx"] = input_df["mol_idx"].fillna(assigned["mol_idx"])

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

    if filter_name == "common_alerts":
        write_common_alert_diagnostics(filter_subdir, extended_df)


def _get_enabled_filters(config_struct_filters):
    """Extract enabled filters from configuration."""
    return {
        k.replace("calculate_", ""): v
        for k, v in config_struct_filters.items()
        if "calculate_" in k and v
    }


def _write_stage_readme(output_dir: Path, config_struct_filters: dict) -> None:
    """Write a README.md into the stage output directory describing available filters and outputs."""
    enabled = sorted(_get_enabled_filters(config_struct_filters).keys())
    hard_filters = sorted(
        key.removeprefix("filter_")
        for key, value in config_struct_filters.items()
        if key.startswith("filter_") and key != "filter_data" and value is True
    )
    all_known = sorted(_FILTER_DESCRIPTIONS.keys(), key=lambda x: x.lower())
    disabled = [f for f in all_known if f not in enabled]

    lines: list[str] = []
    lines.append("# Structural Filters (Stage Output)")
    lines.append("")
    lines.append(
        "This directory is generated by the HEDGEHOG pipeline and contains per-filter results and combined outputs."
    )
    lines.append("")
    lines.append("## Stage Position")
    lines.append("")
    lines.append(
        "- Post-descriptors: runs after descriptor-based filtering and writes results to `03_structural_filters_post/`."
    )
    lines.append("")
    lines.append("## Filters")
    lines.append("")
    lines.append(
        "`calculate_<filter_name>` enables diagnostics; `filter_<filter_name>` independently includes that method in Stage 3 survival."
    )
    lines.append("")
    lines.append("### Enabled in this run")
    lines.append("")
    if enabled:
        for name in enabled:
            desc = _FILTER_DESCRIPTIONS.get(name, "").strip()
            lines.append(f"- `{name}` — {desc}" if desc else f"- `{name}`")
    else:
        lines.append("- (none)")
    lines.append("")
    lines.append("### Hard filters in this run")
    lines.append("")
    if hard_filters:
        for name in hard_filters:
            lines.append(f"- `{name}`")
    else:
        lines.append("- (none; diagnostics only)")
    lines.append("")
    lines.append("### Available filters (reference)")
    lines.append("")
    for name in all_known:
        desc = _FILTER_DESCRIPTIONS.get(name, "").strip()
        lines.append(f"- `{name}` — {desc}" if desc else f"- `{name}`")
    lines.append("")
    lines.append("## Output Structure")
    lines.append("")
    lines.append("Per-filter subdirectories:")
    lines.append("")
    lines.append("- `{filter_name}/metrics.csv` — per-filter summary statistics")
    lines.append("- `{filter_name}/extended.csv` — detailed per-molecule results")
    lines.append(
        "- `{filter_name}/filtered_molecules.csv` — molecules passing that filter"
    )
    lines.append("")
    lines.append("Combined outputs (stage root):")
    lines.append("")
    lines.append(
        "- `filtered_molecules.csv` — molecules passing every method with `filter_<name>: true`"
    )
    lines.append(
        "- `failed_molecules.csv` — molecules failing at least one enforced hard filter (best effort)"
    )
    lines.append(
        "- `structural_liability_profile.csv` — hard decision plus all molecule-level diagnostics"
    )
    lines.append("")
    lines.append("Plots (if generated):")
    lines.append("")
    lines.append("- `plots/molecule_counts_comparison.png`")
    lines.append("- `plots/restriction_ratios_comparison.png`")
    lines.append("")
    lines.append("## Notes")
    lines.append("")
    lines.append(
        "Filter folder names match the config keys without the `calculate_` prefix (e.g., `calculate_common_alerts` -> `common_alerts/`)."
    )
    if disabled:
        preview = ", ".join(f"`{n}`" for n in disabled[:5])
        suffix = "..." if len(disabled) > 5 else ""
        lines.append(
            f"Some available filters may be disabled in this run (e.g., {preview}{suffix})."
        )
    lines.append("")

    (output_dir / "README.md").write_text("\n".join(lines), encoding="utf-8")


def _log_stage_timings(timings: dict[str, float]) -> None:
    """Log structured stage timings."""
    if not timings:
        return
    logger.info("StructFilters timings (seconds):")
    for name, value in timings.items():
        logger.info("  %-28s %.3f", name, value)


def _resolve_stage_options(config_struct_filters, config):
    """Extract and validate stage-level options from config."""
    n_jobs = resolve_n_jobs(config_struct_filters, config)
    write_per_filter_outputs = bool(
        config_struct_filters.get("write_per_filter_outputs", True)
    )
    generate_plots = bool(config_struct_filters.get("generate_plots", True))
    generate_failure_analysis = bool(
        config_struct_filters.get("generate_failure_analysis", True)
    )
    logger.info(
        "StructFilters mode: n_jobs=%s, write_per_filter_outputs=%s, "
        "generate_plots=%s, generate_failure_analysis=%s",
        n_jobs,
        write_per_filter_outputs,
        generate_plots,
        generate_failure_analysis,
    )
    return {
        "n_jobs": n_jobs,
        "write_per_filter_outputs": write_per_filter_outputs,
        "generate_plots": generate_plots,
        "generate_failure_analysis": generate_failure_analysis,
    }


def _make_filter_progress_cb(reporter, molecule_total, filter_name):
    """Build a per-filter progress callback in molecule units."""

    def _alerts_progress(done: int, total: int, name: str = filter_name) -> None:
        progress_total = total if total > 0 else molecule_total
        progress_done = max(0, min(done, progress_total))
        reporter.progress(
            progress_done,
            progress_total,
            message=f"StructFilters: {name}",
        )

    return _alerts_progress


def _compute_filter(
    config,
    apply_func,
    prepared_payload,
    is_csv,
    input_df,
    input_path,
    sample_size,
    progress_cb,
):
    """Dispatch a single filter computation to the appropriate processing path."""
    if prepared_payload is not None:
        return process_prepared_payload(
            config, prepared_payload, apply_func, progress_cb=progress_cb
        )
    if is_csv:
        return process_one_dataframe(
            config, input_df, apply_func, sample_size, progress_cb=progress_cb
        )
    return process_one_file(
        config, input_path, apply_func, sample_size, progress_cb=progress_cb
    )


def _run_post_filter_phases(
    config,
    config_struct_filters,
    stage_dir,
    output_dir,
    input_df,
    pass_mask_by_filter,
    opts,
    timings,
):
    """Execute combine, plot, failure-analysis, and inject phases."""
    combine_started = perf_counter()
    is_single_stage = config.get("_run_single_stage_override") == "struct_filters"
    if config_struct_filters.get("filter_data", False) or is_single_stage:
        combine_filter_results_in_memory(output_dir, input_df, pass_mask_by_filter)
    timings["combine"] = perf_counter() - combine_started

    plot_started = perf_counter()
    if opts["generate_plots"] and opts["write_per_filter_outputs"]:
        plot_calculated_stats(config, stage_dir)
        plot_restriction_ratios(config, stage_dir)
    elif opts["generate_plots"]:
        logger.info("Skipping plots because write_per_filter_outputs is disabled.")
    timings["plots"] = perf_counter() - plot_started

    fail_analysis_started = perf_counter()
    is_post_desc = _is_post_descriptors_stage(stage_dir)
    if (
        is_post_desc
        and opts["generate_failure_analysis"]
        and opts["write_per_filter_outputs"]
    ):
        plot_filter_failures_analysis(config, stage_dir)
    elif is_post_desc and opts["generate_failure_analysis"]:
        logger.info(
            "Skipping failure analysis because write_per_filter_outputs is disabled."
        )
    timings["failure_analysis"] = perf_counter() - fail_analysis_started

    inject_started = perf_counter()
    inject_identity_columns_to_all_csvs(config, stage_dir)
    timings["inject_identity"] = perf_counter() - inject_started


def main(config, stage_dir, reporter=None):
    """Main entry point for structural filters stage.

    Args:
        config: Configuration dictionary
        stage_dir: Stage directory path (e.g., 'stages/03_structural_filters_post')
    """
    sample_size = config.get("sample_size")
    folder_to_save = Path(process_path(config[KEY_FOLDER_TO_SAVE]))
    output_dir = folder_to_save / stage_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    config_struct_filters = load_config(config[CFG_STRUCT_FILTERS])
    _write_stage_readme(output_dir, config_struct_filters)
    if is_large_dataset_mode(config):
        return run_large(
            config,
            stage_dir,
            config_struct_filters,
            output_dir,
            reporter=reporter,
        )

    input_path = _get_input_path(config, stage_dir, folder_to_save)

    try:
        input_df, model_name = _load_input_data(input_path, run_base=folder_to_save)
    except Exception as e:
        logger.error("Could not load input data from %s: %s", input_path, e)
        raise

    if input_df.empty:
        identity_cols = [col for col in IDENTITY_COLUMNS if col in input_df.columns]
        empty_output = input_df[identity_cols].copy()
        empty_output.to_csv(output_dir / "filtered_molecules.csv", index=False)
        empty_output.to_csv(output_dir / "failed_molecules.csv", index=False)
        if config_struct_filters.get("write_structural_liability_profile", True):
            empty_profile = empty_output.copy()
            empty_profile["stage3_hard_pass"] = pd.Series(dtype=bool)
            empty_profile["hard_failed_filters"] = pd.Series(dtype=str)
            empty_profile.to_csv(
                output_dir / "structural_liability_profile.csv", index=False
            )
        logger.info("No molecules available for structural filters; stage is empty.")
        if reporter is not None:
            reporter.progress(1, 1, message="StructFilters complete (empty input)")
        return empty_output

    filters_to_calculate = _get_enabled_filters(config_struct_filters)
    enforced_filters = get_aligned_enforced_filters(config, config_struct_filters)
    if enforced_filters is not None:
        missing_enforced = {
            name
            for name in enforced_filters
            if policy_calculation_filter(name) not in filters_to_calculate
        }
        if missing_enforced:
            names = ", ".join(sorted(missing_enforced))
            raise ValueError(
                f"Structural filter flag enabled but calculation disabled: {names}"
            )
    opts = _resolve_stage_options(config_struct_filters, config)

    is_csv = input_path.lower().endswith(".csv")
    filter_names = list(filters_to_calculate)
    policy_filter_names = get_calculated_policy_names(
        filter_names, config_struct_filters
    )
    molecule_total = max(1, len(input_df))
    timings: dict[str, float] = {}
    stage_started = perf_counter()

    parse_started = perf_counter()
    prepared_payload = None
    if is_csv:
        prepared_payload = prepare_structfilters_input(
            input_df,
            sample_size,
            opts["n_jobs"],
        )
    timings["input_parse"] = perf_counter() - parse_started

    pass_mask_by_filter: dict[str, pd.DataFrame] = {}
    liability_profile = initialize_structural_liability_profile(input_df)

    for filter_name in filter_names:
        if reporter is not None:
            reporter.progress(
                0, molecule_total, message=f"StructFilters: {filter_name}"
            )

        apply_func = filter_function_applier(filter_name)
        progress_cb = None
        if reporter is not None:
            progress_cb = _make_filter_progress_cb(
                reporter, molecule_total, filter_name
            )

        compute_started = perf_counter()
        filter_results = _compute_filter(
            config,
            apply_func,
            prepared_payload,
            is_csv,
            input_df,
            input_path,
            sample_size,
            progress_cb,
        )
        timings[f"filter_compute:{filter_name}"] = perf_counter() - compute_started

        if filter_results is None:
            logger.warning("No molecules to process for filter: %s", filter_name)
            continue

        post_started = perf_counter()
        final_res, final_extended = get_basic_stats(
            config_struct_filters, filter_results, model_name, filter_name=filter_name
        )
        final_extended, enforcement_mask = attach_structural_enforcement_pass(
            config_struct_filters, filter_name, final_extended
        )
        liability_profile = merge_structural_liability_profile(
            liability_profile, filter_name, final_extended
        )
        liability_profile = merge_structural_policy_aliases(
            liability_profile, filter_name
        )
        policy_masks = build_structural_policy_pass_masks(
            config_struct_filters,
            filter_name,
            final_extended,
            default_mask=enforcement_mask,
        )
        for policy_name, policy_mask in policy_masks.items():
            if enforced_filters is None or policy_name in enforced_filters:
                pass_mask_by_filter[policy_name] = policy_mask
            else:
                logger.info(
                    "Calculated structural policy %s; current enforcement policy excludes it from rejection.",
                    policy_name,
                )
        if opts["write_per_filter_outputs"]:
            _save_filter_results(output_dir, filter_name, final_res, final_extended)
        timings[f"filter_post:{filter_name}"] = perf_counter() - post_started
        if reporter is not None:
            reporter.progress(
                molecule_total,
                molecule_total,
                message=f"StructFilters: {filter_name}",
            )

    if config_struct_filters.get("write_structural_liability_profile", True):
        liability_profile = finalize_structural_liability_profile(
            liability_profile, enforced_filters, policy_filter_names
        )
        liability_profile.to_csv(
            output_dir / "structural_liability_profile.csv", index=False
        )

    _run_post_filter_phases(
        config,
        config_struct_filters,
        stage_dir,
        output_dir,
        input_df,
        pass_mask_by_filter,
        opts,
        timings,
    )

    timings["total"] = perf_counter() - stage_started
    _log_stage_timings(timings)

    if reporter is not None:
        reporter.progress(
            molecule_total, molecule_total, message="StructFilters complete"
        )
