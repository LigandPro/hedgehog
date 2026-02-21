from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

from hedgehog.configs.logger import load_config, logger
from hedgehog.utils.paths import process_path


def format_number(x, pos=None):
    """Format number for display. pos parameter is for matplotlib FuncFormatter compatibility."""
    if x >= 1e6:
        return f"{x / 1e6:.1f}M"
    elif x >= 1e3:
        return f"{x / 1e3:.1f}K"
    return f"{x:.0f}"


def get_model_colors(model_names, cmap=None):
    """Generate color map for models."""
    if cmap is None:
        colors = plt.cm.YlOrRd(np.linspace(1, 0, len(model_names) + 1))
    else:
        colors = plt.colormaps.get_cmap(cmap)(np.linspace(1, 0, len(model_names) + 1))

    return dict(zip(model_names, colors, strict=False))


def clean_name(name):
    """Clean metric names for display."""
    for pattern in ["metrics", "_", ".csv"]:
        name = name.replace(pattern, "")
    return name.strip()


def check_paths(config, paths):
    all_filters = {}
    for k, v in config.items():
        if "calculate_" in k:
            k = k.replace("calculate_", "")
            all_filters[k] = v

    path_folders = []
    for path in paths:
        parts = path.split("/")
        if len(parts) >= 2:
            folder_name = parts[-2]
            path_folders.append(folder_name.lower())

    missing_filters = []
    for filter_name, enabled in all_filters.items():
        if enabled:
            filter_name_lower = filter_name.lower()
            filter_name_no_underscore = filter_name_lower.replace("_", "")
            found = any(
                filter_name_lower == folder
                or filter_name_no_underscore == folder.replace("_", "")
                or filter_name_lower in folder
                or filter_name_no_underscore in folder.replace("_", "")
                for folder in path_folders
            )
            if not found:
                missing_filters.append(filter_name_no_underscore)

    if len(missing_filters) > 0:
        raise AssertionError(
            f"Invalid filter name(s) missing: {', '.join(missing_filters)}"
        )
    return True


def plot_calculated_stats(config, stage_dir):
    """Plot calculated statistics for structural filters."""
    folder_to_save = Path(process_path(config["folder_to_save"]))
    config_structFilters = load_config(config["config_structFilters"])

    struct_folder = folder_to_save / stage_dir
    paths = list(struct_folder.glob("*/metrics.csv"))
    if not paths:
        paths = list(struct_folder.glob("*metrics.csv"))
    paths = [str(p) for p in paths]
    check_paths(config_structFilters, paths)

    datas = []
    filter_names = []
    all_model_names = set()

    for path in paths:
        data = pd.read_csv(path)

        all_model_names.update(data["model_name"].dropna().unique())
        data.set_index("model_name", inplace=True)

        banned_cols = [col for col in data.columns if "banned_ratio" in col]
        data_filtered = data[banned_cols + ["num_mol"]].copy()
        for banned_col in banned_cols:
            data_filtered.loc[:, f"num_banned_{banned_col}"] = (
                data_filtered[banned_col] * data_filtered["num_mol"]
            )
        datas.append(data_filtered)

        filter_name = Path(path).parent.name
        filter_names.append(filter_name)

    model_name_set = sorted(list(all_model_names))

    filter_results = {}
    filters_to_find = list(struct_folder.glob("*/filtered_molecules.csv"))
    if not filters_to_find:
        filters_to_find = list(struct_folder.glob("*filteredMols.csv"))
    filters_to_find = [str(p) for p in filters_to_find]

    for path in filters_to_find:
        try:
            filter_data = pd.read_csv(path)
            filter_name = Path(path).parent.name

            num_passed_by_model = None
            if "pass" in filter_data.columns:
                passed = filter_data[filter_data["pass"]]
                if len(passed) > 0:
                    num_passed_by_model = passed.groupby("model_name").size().to_dict()

            if num_passed_by_model is not None:
                filter_results[filter_name] = num_passed_by_model
            else:
                default_models = filter_data["model_name"].unique().tolist()
                filter_results[filter_name] = {m: 0 for m in default_models}

        except (IndexError, FileNotFoundError) as e:
            logger.warning("Could not process %s: %s", path, e)
            filter_results[filter_name] = {}

    all_models = model_name_set

    for filter_name, values in filter_results.items():
        if len(values) != len(all_models):
            for model in all_models:
                if model not in values:
                    filter_results[filter_name][model] = 0
        filter_results[filter_name] = dict(sorted(filter_results[filter_name].items()))

    n_plots = len(datas)
    n_cols = 2
    n_rows = (n_plots + n_cols - 1) // n_cols

    plt.figure(figsize=(40, 5 * n_rows))
    for idx, (data, filter_name) in enumerate(zip(datas, filter_names, strict=False)):
        ax = plt.subplot(n_rows, n_cols, idx + 1)
        models = data.index
        x = np.arange(len(models))
        width = 0.8
        total_mols = data["num_mol"].sum()
        total = ax.barh(
            x,
            data.loc[models, "num_mol"],
            width,
            label=f"Total Molecules ({format_number(total_mols)})",
            color="#E5E5E5",
            alpha=0.5,
        )

        clean_filter_name = filter_name.split("/")[-1].lower()
        for known_filter in filter_results.keys():
            if known_filter.lower() in clean_filter_name:
                for i, (model, passed) in enumerate(
                    filter_results[known_filter].items()
                ):
                    bar_center_x = data.loc[models, "num_mol"].values[0] / 2
                    bar_center_y = x[i]
                    model_total = (
                        data.loc[model, "num_mol"]
                        if model in data.index
                        else data["num_mol"].iloc[0]
                    )
                    if model_total != 0:
                        text = f"Passed molecules: {passed} ({(passed / model_total * 100):.1f}%)"
                    else:
                        text = f"Passed molecules: {passed} (0%)"
                    ax.annotate(
                        text,
                        (bar_center_x, bar_center_y),
                        ha="center",
                        va="center",
                        fontsize=12,
                        color="black",
                        fontweight="bold",
                        bbox=dict(
                            facecolor="white", alpha=0.7, edgecolor="none", pad=3
                        ),
                        zorder=1000,
                    )

        for i, model in enumerate(models):
            count = data.loc[model, "num_mol"]
            max_bar_width = data["num_mol"].max()
            text_x_position = max_bar_width
            ax.text(
                text_x_position,
                i,
                int(count),
                va="center",
                ha="left",
                fontsize=12,
                color="black",
                fontweight="bold",
            )

        banned_bars = []
        banned_percentages = []
        ratio_cols = [
            col
            for col in data.columns
            if "banned_ratio" in col and "num_banned" not in col
        ]
        colors = get_model_colors(model_names=ratio_cols, cmap="Paired")
        for col, color in zip(ratio_cols, colors.values(), strict=False):
            num_banned_col = f"num_banned_{col}"
            ratio_name = col.replace("banned_ratio", "").strip("_")
            ratio_name = clean_name(ratio_name)

            banned_count = data[num_banned_col]
            total_banned = banned_count.sum()
            banned_percent = (total_banned / total_mols) * 100 if total_mols > 0 else 0
            banned_percentages.append(banned_percent)

            if banned_percent == 0.0:
                label = f"{ratio_name} (0%)"
            else:
                label = f"{ratio_name} ({format_number(total_banned)}, {banned_percent:.1f}%)"

            bar = ax.barh(x, banned_count, width, label=label, color=color, alpha=0.8)
            banned_bars.append(bar)

        clean_filter_name = filter_name.split("/")[-1]
        ax.set_title(
            clean_name(clean_filter_name), fontsize=14, pad=20, fontweight="bold"
        )
        ax.set_yticks(x)
        ax.set_yticklabels(models, fontsize=12)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(format_number))
        ax.set_xlim(left=0)

    ax.set_xlabel(
        f"Number of Molecules (Total: {format_number(total_mols)})",
        fontsize=12,
        labelpad=10,
    )
    ax.set_ylabel("Models", fontsize=12, labelpad=10)

    ax.grid(True, axis="x", alpha=0.2, linestyle="--")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    sorted_indices = np.argsort(banned_percentages)[::-1]
    sorted_handles = [total] + [banned_bars[i] for i in sorted_indices]

    legend = ax.legend(
        handles=sorted_handles,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=11,
        ncol=1,
    )
    legend.get_frame().set_alpha(0.9)
    legend.get_frame().set_edgecolor("lightgray")

    plt.subplots_adjust(right=0.85, hspace=0.6, wspace=0.5)

    plots_dir = struct_folder / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(
        plots_dir / "molecule_counts_comparison.png",
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
        edgecolor="none",
    )
    plt.close()


def plot_restriction_ratios(config, stage_dir):
    """Plot restriction ratios for structural filters."""
    folder_to_save = Path(process_path(config["folder_to_save"]))
    folder_name = config["folder_to_save"].split("/")[-1]

    config_structFilters = load_config(config["config_structFilters"])

    struct_folder = folder_to_save / stage_dir
    paths = list(struct_folder.glob("*/metrics.csv"))
    if not paths:
        paths = list(struct_folder.glob("*metrics.csv"))
    paths = [str(p) for p in paths]
    check_paths(config_structFilters, paths)

    if not paths:
        logger.error("No data files found in %s", str(folder_to_save))
        return

    filter_data = {}
    model_names_filters = {}
    for path in paths:
        filter_name = path.split(f"{folder_name}/")[-1].split("_metrics.csv")[0]
        data = pd.read_csv(path)

        ratio_cols = [col for col in data.columns if "banned_ratio" in col]
        model_names_filters[filter_name] = dict(
            zip(data["model_name"].tolist(), data["num_mol"].tolist(), strict=False)
        )

        if not ratio_cols:
            continue

        clean_cols = {
            col: col.replace("_banned_ratio", "")
            .replace("banned_ratio", "")
            .replace("_s", "s")
            for col in ratio_cols
        }
        ratios = data[ratio_cols].rename(columns=clean_cols)
        actual_model_names = data["model_name"].tolist()
        ratios.index = actual_model_names

        row = ratios.iloc[0]
        if row.isna().all():
            continue

        all_value = None
        if "all" in row.index:
            all_value = row["all"]
            row = row.drop("all")

        sorted_values = row.sort_values(ascending=False)

        if all_value is not None:
            if all_value >= sorted_values.iloc[0]:
                sorted_index = pd.Index(["all"]).append(sorted_values.index)
            else:
                sorted_index = sorted_values.index.append(pd.Index(["all"]))
            ratios = ratios[sorted_index]
        else:
            ratios = ratios[sorted_values.index]

        filter_data[filter_name] = ratios
    if not filter_data:
        logger.error("No valid data to plot")
        return

    model_names_filters = pd.DataFrame(model_names_filters).reset_index()
    model_names_filters = model_names_filters.rename(columns={"index": "model_name"})

    plt.style.use("default")
    sns.set_style("white")
    sns.set_context("talk")

    n_filters = len(filter_data)
    n_cols = min(2, n_filters)
    n_rows = (n_filters + n_cols - 1) // n_cols

    fig = plt.figure(figsize=(16, 7 * n_rows))
    fig.suptitle(
        "Comparison of Restriction Ratios Across Different Filters",
        fontsize=16,
        y=0.98,
        fontweight="bold",
    )

    for idx, (filter_name, data) in enumerate(filter_data.items()):
        number_of_mols = np.array(model_names_filters[filter_name].tolist())
        ax = plt.subplot(n_rows, n_cols, idx + 1)
        for col in data.columns:
            if col not in ["any", "all", "model_name", "num_mol"]:
                data[col] = number_of_mols * (1 - np.array(data[col].tolist()))
        if not data.empty and data.notna().any().any():
            if "all" in data.columns:
                data.drop(columns=["all"], inplace=True)
            if "any" in data.columns:
                data.drop(columns=["any"], inplace=True)

            custom_cmap = LinearSegmentedColormap.from_list(
                "custom", ["white", "#B29EEE"]
            )
            show_cbar = idx == 1
            sns.heatmap(
                data.T,
                cmap=custom_cmap,
                ax=ax,
                cbar_kws={"label": "Passed Molecules", "format": "%d"},
                vmin=0,
                vmax=max(data.max()),
                fmt=".0f",
                annot=True,
                annot_kws={"size": 12, "rotation": 0, "color": "black"},
                cbar=show_cbar,
            )

            ax.set_title(
                f"{clean_name(filter_name)} Filter", fontsize=12, fontweight="bold"
            )
            plt.setp(ax.get_yticklabels(), rotation=0, ha="right", fontsize=12)
            plt.setp(ax.get_xticklabels(), rotation=0, ha="right", fontsize=12)
            ax.set_xlabel("Model")

            actual_model_names = data.index.tolist()
            if len(actual_model_names) == len(ax.get_xticklabels()):
                ax.set_xticklabels(actual_model_names)

        else:
            ax.text(
                0.5,
                0.5,
                "No data available",
                horizontalalignment="center",
                verticalalignment="center",
                transform=ax.transAxes,
            )
            ax.set_title(
                f"{clean_name(filter_name)} Filter",
                pad=10,
                fontsize=11,
                fontweight="bold",
            )

    plt.tight_layout()
    plots_dir = struct_folder / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(
        plots_dir / "restriction_ratios_comparison.png",
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
        edgecolor="none",
    )
    plt.close()


def _get_breakdown_folder(file_path):
    """Get the CommonAlertsBreakdown folder path, creating it if necessary."""
    path_to_save = Path(file_path).parent / "CommonAlertsBreakdown"
    path_to_save.mkdir(parents=True, exist_ok=True)
    return str(path_to_save)


def _save_plot(file_path, filename, dpi=600):
    """Save plot to CommonAlertsBreakdown folder and close it."""
    output_path = Path(_get_breakdown_folder(file_path)) / filename
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close()


def _aggregate_reasons(filter_reasons):
    """Aggregate all reasons across filters into a single dictionary."""
    all_reasons = {}
    for reasons in filter_reasons.values():
        for reason, count in reasons:
            all_reasons[reason] = all_reasons.get(reason, 0) + count
    return all_reasons


def _truncate_string(s, max_length, suffix="..."):
    """Truncate string if longer than max_length."""
    if len(s) > max_length:
        return s[: max_length - len(suffix)] + suffix
    return s


def _calculate_grid_layout(num_items):
    """Calculate rows and columns for grid layout."""
    if num_items <= 3:
        return 1, num_items
    if num_items <= 6:
        return 2, 3
    if num_items <= 9:
        return 3, 3
    if num_items <= 12:
        return 3, 4
    cols = 4
    rows = (num_items + cols - 1) // cols
    return rows, cols


def analyze_filter_failures(file_path):
    """Analyze filter failures from extended CSV file and generate visualizations."""
    logger.debug("Analyzing filter failures from: %s", file_path)
    df = pd.read_csv(file_path, low_memory=False)

    filter_columns = [
        col
        for col in df.columns
        if col.startswith("pass_") and col not in ("pass", "pass_any")
    ]

    if not filter_columns:
        return None, None, None

    filter_failures = {}
    filter_reasons = {}
    all_detailed_reasons = {}

    for col in filter_columns:
        filter_name = col.replace("pass_", "")
        failures = (~df[col]).sum()
        total = len(df)

        filter_failures[filter_name] = {
            "failures": failures,
            "total": total,
            "percentage": (failures / total) * 100,
        }

        reasons_col = f"reasons_{filter_name}"
        if reasons_col in df.columns:
            failed_molecules = df[~df[col]]
            reasons_data = failed_molecules[reasons_col].dropna()

            reason_counts = {}
            for reasons_str in reasons_data:
                if pd.notna(reasons_str) and str(reasons_str).strip():
                    for reason in str(reasons_str).split(";"):
                        reason = reason.strip()
                        if reason:
                            reason_counts[reason] = reason_counts.get(reason, 0) + 1

            filter_reasons[filter_name] = sorted(
                reason_counts.items(), key=lambda x: x[1], reverse=True
            )
            all_detailed_reasons[filter_name] = reason_counts

    _create_main_filter_plot(filter_failures, file_path)
    _create_individual_filter_plots(filter_failures, filter_reasons, file_path)
    _create_multi_panel_filter_plot(filter_failures, filter_reasons, file_path)

    all_reasons = _aggregate_reasons(filter_reasons)
    top_reasons = sorted(all_reasons.items(), key=lambda x: x[1], reverse=True)[:5]

    if top_reasons:
        logger.info(
            "Top 5 most common filter failure reasons (molecules may have multiple reasons):"
        )
        for i, (reason, count) in enumerate(top_reasons, 1):
            logger.info("  %d. %s: %d failures", i, reason, count)

    _create_complete_reasons_breakdown(all_detailed_reasons, filter_failures, file_path)
    _create_comprehensive_overview(filter_reasons, file_path)
    _create_summary_table(filter_failures, filter_reasons, file_path)

    return filter_failures, filter_reasons, all_detailed_reasons


def _create_main_filter_plot(filter_failures, file_path):
    """Create main filter failures bar chart."""
    plot_data = [
        {
            "filter": name,
            "failures": stats["failures"],
            "percentage": stats["percentage"],
        }
        for name, stats in filter_failures.items()
    ]
    plot_df = pd.DataFrame(plot_data).sort_values("failures", ascending=False)

    plt.figure(figsize=(max(16, len(plot_df) * 0.6), 16))
    plt.bar(
        range(len(plot_df)),
        plot_df["failures"],
        color="steelblue",
        alpha=0.8,
        width=0.3,
    )

    plt.xlabel("Filters", fontsize=20)
    plt.ylabel("Number of Molecules Failed", fontsize=20)
    plt.title(
        "Number of Molecules Failed by Each Filter", fontsize=26, fontweight="bold"
    )
    plt.xticks(
        range(len(plot_df)), plot_df["filter"], rotation=45, ha="right", fontsize=16
    )

    max_failures = max(plot_df["failures"])
    for i, (_, row) in enumerate(plot_df.iterrows()):
        plt.text(
            i,
            row["failures"] + max_failures * 0.01,
            f"{row['failures']}\n({row['percentage']:.1f}%)",
            ha="center",
            va="bottom",
            fontsize=14,
        )

    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    _save_plot(file_path, "filter_failures_plot.png")


def _create_individual_filter_plots(filter_failures, filter_reasons, file_path):
    """Create individual plots for each filter showing failure reasons."""
    for filter_name, stats in filter_failures.items():
        if stats["failures"] == 0:
            continue

        reasons_data = filter_reasons.get(filter_name, [])
        if not reasons_data:
            continue

        plot_data = [
            {
                "Reason": reason,
                "Count": count,
                "Percentage": (count / stats["failures"]) * 100,
            }
            for reason, count in reasons_data
        ]
        plot_df = pd.DataFrame(plot_data).sort_values("Count", ascending=False)

        plt.figure(figsize=(max(16, len(plot_df) * 0.6), 20))
        plt.bar(
            range(len(plot_df)),
            plot_df["Count"],
            color="steelblue",
            alpha=0.8,
            width=0.3,
        )

        plt.xlabel("Failure Reasons", fontsize=20)
        plt.ylabel("Number of Molecules Failed", fontsize=20)
        plt.title(
            f"{filter_name.upper()} - Failure Reasons ({len(plot_df)} reasons, {stats['failures']} total failures)",
            fontsize=26,
            fontweight="bold",
        )
        plt.xticks(
            range(len(plot_df)),
            plot_df["Reason"],
            rotation=45,
            ha="right",
            fontsize=max(10, min(16, 300 // len(plot_df))),
        )

        max_count = max(plot_df["Count"])
        for i, (_, row) in enumerate(plot_df.iterrows()):
            plt.text(
                i,
                row["Count"] + max_count * 0.01,
                f"{row['Count']}\n({row['Percentage']:.1f}%)",
                ha="center",
                va="bottom",
                fontsize=12,
            )

        plt.grid(axis="y", alpha=0.3)
        plt.tight_layout()
        _save_plot(file_path, f"{filter_name}_reasons_plot.png")


def _create_multi_panel_filter_plot(filter_failures, filter_reasons, file_path):
    """Create multi-panel plot showing all filters with reasons."""
    sorted_filters = [
        (name, stats)
        for name, stats in sorted(
            filter_failures.items(), key=lambda x: x[1]["failures"], reverse=True
        )
        if stats["failures"] > 0
    ]

    if not sorted_filters:
        return

    rows, cols = _calculate_grid_layout(len(sorted_filters))
    plt.figure(figsize=(cols * 6, rows * 6))

    for i, (filter_name, _stats) in enumerate(sorted_filters):
        all_reasons_data = filter_reasons.get(filter_name, [])
        reason_names = [r[0] for r in all_reasons_data]
        reason_counts = [r[1] for r in all_reasons_data]

        if len(reason_names) > 10:
            title_suffix = f"(Top 10 of {len(all_reasons_data)} reasons)"
            reason_names = reason_names[:10]
            reason_counts = reason_counts[:10]
        else:
            title_suffix = f"({len(all_reasons_data)} reasons)"

        plt.subplot(rows, cols, i + 1)
        plt.bar(
            range(len(reason_names)),
            reason_counts,
            color="steelblue",
            alpha=0.8,
            width=0.3,
        )

        plt.xlabel("Reasons", fontsize=14)
        plt.ylabel("Molecules Failed", fontsize=14)
        plt.title(
            f"{filter_name.upper()}\n{title_suffix}", fontsize=16, fontweight="bold"
        )

        truncated_names = [_truncate_string(name, 15, "...") for name in reason_names]
        plt.xticks(
            range(len(truncated_names)),
            truncated_names,
            rotation=45,
            ha="right",
            fontsize=12,
        )

        max_count = max(reason_counts) if reason_counts else 0
        for j, count in enumerate(reason_counts):
            plt.text(
                j,
                count + max_count * 0.01,
                f"{count}",
                ha="center",
                va="bottom",
                fontsize=11,
            )

        plt.grid(axis="y", alpha=0.3)

    plt.tight_layout(h_pad=1.5, w_pad=1.0)
    _save_plot(file_path, "all_filters_reasons_plot.png")


def _create_complete_reasons_breakdown(
    all_detailed_reasons, filter_failures, file_path
):
    """Create complete CSV breakdown of all reasons."""
    breakdown_data = []
    for filter_name, reasons_dict in all_detailed_reasons.items():
        total_failures = filter_failures[filter_name]["failures"]
        for reason, count in sorted(
            reasons_dict.items(), key=lambda x: x[1], reverse=True
        ):
            breakdown_data.append(
                {
                    "Ruleset": filter_name,
                    "Reason": reason,
                    "Count": count,
                    "Percentage_of_Filter_Failures": (count / total_failures) * 100
                    if total_failures > 0
                    else 0,
                    "Total_Filter_Failures": total_failures,
                }
            )

    breakdown_df = pd.DataFrame(breakdown_data)
    output_path = (
        Path(_get_breakdown_folder(file_path)) / "complete_reasons_breakdown.csv"
    )
    breakdown_df.to_csv(output_path, index=False)
    return breakdown_df


def _create_comprehensive_overview(filter_reasons, file_path):
    """Create comprehensive overview plot of most common failure reasons."""
    all_reasons = _aggregate_reasons(filter_reasons)
    top_reasons = sorted(all_reasons.items(), key=lambda x: x[1], reverse=True)

    if not top_reasons:
        return

    display_count = min(30, len(top_reasons))
    reason_names = [
        _truncate_string(r[0], 30, "...") for r in top_reasons[:display_count]
    ]
    reason_counts = [r[1] for r in top_reasons[:display_count]]

    plt.figure(figsize=(max(16, display_count * 0.6), 16))
    plt.bar(
        range(len(reason_names)), reason_counts, color="darkgreen", alpha=0.7, width=0.3
    )

    plt.xlabel("Failure Reasons", fontsize=20)
    plt.ylabel("Total Number of Molecules Failed", fontsize=20)
    plt.title(
        f"Most Common Molecular Filter Failure Reasons (Top {display_count} of {len(top_reasons)})",
        fontsize=26,
        fontweight="bold",
    )
    plt.xticks(
        range(len(reason_names)), reason_names, rotation=45, ha="right", fontsize=16
    )

    max_count = max(reason_counts)
    for i, count in enumerate(reason_counts):
        plt.text(
            i,
            count + max_count * 0.01,
            f"{count}",
            ha="center",
            va="bottom",
            fontsize=14,
        )

    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    _save_plot(file_path, "comprehensive_reasons_overview.png")

    # Save CSV summary
    all_reasons_df = pd.DataFrame(top_reasons, columns=["Reason", "Total_Count"])
    output_path = Path(_get_breakdown_folder(file_path)) / "all_reasons_summary.csv"
    all_reasons_df.to_csv(output_path, index=False)


def _create_summary_table(filter_failures, filter_reasons, file_path):
    """Create summary table CSV with filter statistics."""
    summary_data = []
    for filter_name, stats in filter_failures.items():
        row = {
            "Ruleset": filter_name,
            "Total_Failures": stats["failures"],
            "Failure_Percentage": stats["percentage"],
            "Total_Molecules": stats["total"],
            "Unique_Reasons_Count": len(filter_reasons.get(filter_name, [])),
        }

        reasons = filter_reasons.get(filter_name, [])
        for i, (reason, count) in enumerate(reasons[:5], 1):
            row[f"Top_Reason_{i}"] = reason
            row[f"Top_Reason_{i}_Count"] = count
            row[f"Top_Reason_{i}_Percentage"] = (
                (count / stats["failures"]) * 100 if stats["failures"] > 0 else 0
            )

        summary_data.append(row)

    summary_df = pd.DataFrame(summary_data).sort_values(
        "Total_Failures", ascending=False
    )
    output_path = Path(_get_breakdown_folder(file_path)) / "filter_summary_table.csv"
    summary_df.to_csv(output_path, index=False)
    logger.debug("Summary table saved to: %s", output_path)
    return summary_df


def plot_filter_failures_analysis(config, stage_dir):
    """Analyze and plot filter failures for extended CSV files."""
    folder_to_save = Path(process_path(config["folder_to_save"]))
    struct_folder = folder_to_save / stage_dir

    if not struct_folder.exists():
        return

    extended_files = [str(p) for p in struct_folder.glob("*_extended.csv")]

    if not extended_files:
        logger.debug("No extended CSV files found for failure analysis")
        return

    for file_path in extended_files:
        try:
            analyze_filter_failures(file_path)
        except (OSError, pd.errors.ParserError, KeyError, ValueError) as e:
            logger.debug("Error analyzing filter failures for %s: %s", file_path, e)
