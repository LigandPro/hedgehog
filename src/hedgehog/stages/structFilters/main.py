from pathlib import Path
from time import perf_counter

import pandas as pd

from hedgehog.configs.logger import load_config, logger
from hedgehog.stages.structFilters.utils import (
    combine_filter_results_in_memory,
    filter_data,
    filter_function_applier,
    get_basic_stats,
    inject_identity_columns_to_all_csvs,
    plot_calculated_stats,
    plot_filter_failures_analysis,
    plot_restriction_ratios,
    prepare_structfilters_input,
    process_prepared_payload,
    process_one_dataframe,
    process_one_file,
    process_path,
)
from hedgehog.utils.input_paths import find_sampled_molecules
from hedgehog.utils.parallel import resolve_n_jobs

IDENTITY_COLUMNS = ["smiles", "model_name", "mol_idx"]

_FILTER_DESCRIPTIONS: dict[str, str] = {
    "common_alerts": "SMARTS-based structural alert screening using curated rule sets (PAINS, Dundee, BMS, Glaxo, etc.).",
    "molgraph_stats": "Graph-theoretic severity metrics (levels 1–11) highlighting topological anomalies.",
    "molcomplexity": "Molecular complexity heuristics (medchem complexity filters).",
    "NIBR": "Novartis in-house structural filters (severity-based developability heuristics).",
    "bredt": "Bredt's rule violations at bridgehead positions in small bicyclic systems.",
    "lilly": "Eli Lilly Medchem Rules demerit scoring; higher demerits indicate less desirable structures.",
    "protecting_groups": "Common protecting group motifs (e.g., Boc/Fmoc/Cbz) that are undesirable in final compounds.",
    "ring_infraction": "Strained/unusual ring systems and ring infractions (configurable heterocycle minimum size).",
    "stereo_center": "Excessive or undefined stereocenters (configurable max count / max undefined).",
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
    if "config_descriptors" in config:
        try:
            desc_config = load_config(config["config_descriptors"])
            descriptors_enabled = desc_config.get("run", False)
        except Exception:
            pass

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

        mol_prep_path = base / "stages" / "00_mol_prep" / "filtered_molecules.csv"
        if mol_prep_path.exists() and mol_prep_path.stat().st_size > 0:
            return str(mol_prep_path)

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


def _write_stage_readme(
    output_dir: Path, stage_dir: str, config_structFilters: dict
) -> None:
    """Write a README.md into the stage output directory describing available filters and outputs."""
    enabled = sorted(_get_enabled_filters(config_structFilters).keys())
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
        "Each filter is controlled via `config_structFilters.yml` keys like `calculate_<filter_name>: true`."
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
    lines.append("- `filtered_molecules.csv` — molecules passing all enabled filters")
    lines.append(
        "- `failed_molecules.csv` — molecules failing at least one filter (best effort)"
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


def _resolve_parse_input_n_jobs(config_structFilters: dict, config: dict) -> int:
    """Resolve workers for one-time SMILES parsing."""
    raw = config_structFilters.get("parse_input_n_jobs", -1)
    try:
        value = int(raw)
    except (TypeError, ValueError):
        logger.warning(
            "Invalid parse_input_n_jobs=%r. Falling back to auto.", raw
        )
        value = -1

    if value > 0:
        return value
    return resolve_n_jobs({"n_jobs": value}, config)


def _build_filter_pass_mask(filter_extended: pd.DataFrame) -> pd.DataFrame:
    """Build compact pass mask table for in-memory combine."""
    id_cols = [c for c in IDENTITY_COLUMNS if c in filter_extended.columns]
    if "pass" in filter_extended.columns:
        pass_col = "pass"
    elif "pass_filter" in filter_extended.columns:
        pass_col = "pass_filter"
    else:
        pass_col = None

    if pass_col is None or not id_cols:
        return pd.DataFrame(columns=id_cols + ["pass"])

    out = filter_extended[id_cols + [pass_col]].copy()
    if pass_col != "pass":
        out = out.rename(columns={pass_col: "pass"})
    out["pass"] = out["pass"].fillna(False).astype(bool)
    return out.drop_duplicates(subset=id_cols, keep="last")


def _log_stage_timings(timings: dict[str, float]) -> None:
    """Log structured stage timings."""
    if not timings:
        return
    logger.info("StructFilters timings (seconds):")
    for name, value in timings.items():
        logger.info("  %-28s %.3f", name, value)


def main(config, stage_dir, reporter=None):
    """Main entry point for structural filters stage.

    Args:
        config: Configuration dictionary
        stage_dir: Stage directory path (e.g., 'stages/03_structural_filters_post')
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
    _write_stage_readme(output_dir, stage_dir, config_structFilters)
    filters_to_calculate = _get_enabled_filters(config_structFilters)

    parse_input_n_jobs = _resolve_parse_input_n_jobs(config_structFilters, config)
    write_per_filter_outputs = bool(
        config_structFilters.get("write_per_filter_outputs", True)
    )
    generate_plots = bool(config_structFilters.get("generate_plots", True))
    generate_failure_analysis = bool(
        config_structFilters.get("generate_failure_analysis", True)
    )
    combine_in_memory = bool(config_structFilters.get("combine_in_memory", True))
    if not write_per_filter_outputs and not combine_in_memory:
        logger.warning(
            "write_per_filter_outputs=false with combine_in_memory=false is invalid. "
            "Forcing combine_in_memory=true."
        )
        combine_in_memory = True

    logger.info(
        "StructFilters mode: parse_input_n_jobs=%s, write_per_filter_outputs=%s, "
        "generate_plots=%s, generate_failure_analysis=%s, combine_in_memory=%s",
        parse_input_n_jobs,
        write_per_filter_outputs,
        generate_plots,
        generate_failure_analysis,
        combine_in_memory,
    )

    is_csv = input_path.lower().endswith(".csv")
    filter_names = list(filters_to_calculate)
    stage_total = max(1, len(filter_names) * 100)
    timings: dict[str, float] = {}
    stage_started = perf_counter()

    parse_started = perf_counter()
    prepared_payload = None
    if is_csv:
        prepared_payload = prepare_structfilters_input(
            input_df,
            sample_size,
            parse_input_n_jobs,
        )
    timings["input_parse"] = perf_counter() - parse_started

    pass_mask_by_filter: dict[str, pd.DataFrame] = {}

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

        compute_started = perf_counter()
        if prepared_payload is not None:
            filter_results = process_prepared_payload(
                config, prepared_payload, apply_func, progress_cb=progress_cb
            )
        elif is_csv:
            filter_results = process_one_dataframe(
                config, input_df, apply_func, sample_size, progress_cb=progress_cb
            )
        else:
            filter_results = process_one_file(
                config, input_path, apply_func, sample_size, progress_cb=progress_cb
            )
        timings[f"filter_compute:{filter_name}"] = perf_counter() - compute_started

        if filter_results is None:
            logger.warning("No molecules to process for filter: %s", filter_name)
            continue

        post_started = perf_counter()
        final_res, final_extended = get_basic_stats(
            config_structFilters, filter_results, model_name, filter_name=filter_name
        )
        pass_mask_by_filter[filter_name] = _build_filter_pass_mask(final_extended)
        if write_per_filter_outputs:
            _save_filter_results(output_dir, filter_name, final_res, final_extended)
        timings[f"filter_post:{filter_name}"] = perf_counter() - post_started
        if reporter is not None:
            reporter.progress(
                base + 100, stage_total, message=f"StructFilters: {filter_name}"
            )

    combine_started = perf_counter()
    is_single_stage = config.get("_run_single_stage_override") == "struct_filters"
    if config_structFilters.get("filter_data", False) or is_single_stage:
        if combine_in_memory:
            combine_filter_results_in_memory(output_dir, input_df, pass_mask_by_filter)
        else:
            filter_data(config, stage_dir)
    timings["combine"] = perf_counter() - combine_started

    plot_started = perf_counter()
    if generate_plots and write_per_filter_outputs:
        plot_calculated_stats(config, stage_dir)
        plot_restriction_ratios(config, stage_dir)
    elif generate_plots and not write_per_filter_outputs:
        logger.info(
            "Skipping plots because write_per_filter_outputs is disabled."
        )
    timings["plots"] = perf_counter() - plot_started

    fail_analysis_started = perf_counter()
    if (
        _is_post_descriptors_stage(stage_dir)
        and generate_failure_analysis
        and write_per_filter_outputs
    ):
        plot_filter_failures_analysis(config, stage_dir)
    elif _is_post_descriptors_stage(stage_dir) and generate_failure_analysis:
        logger.info(
            "Skipping failure analysis because write_per_filter_outputs is disabled."
        )
    timings["failure_analysis"] = perf_counter() - fail_analysis_started

    inject_started = perf_counter()
    inject_identity_columns_to_all_csvs(config, stage_dir)
    timings["inject_identity"] = perf_counter() - inject_started
    timings["total"] = perf_counter() - stage_started
    _log_stage_timings(timings)

    if reporter is not None:
        reporter.progress(stage_total, stage_total, message="StructFilters complete")
