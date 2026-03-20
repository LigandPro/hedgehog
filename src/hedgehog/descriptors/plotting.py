"""Chart and graph generation for descriptor distributions."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from hedgehog._constants import CFG_DESCRIPTORS
from hedgehog.configs.logger import load_config
from hedgehog.descriptors.constants import _DESCRIPTOR_KEY_MAP
from hedgehog.descriptors.filtering import (
    _build_structural_constraint_borders,
    _extract_borders_and_constraints,
    _get_border_values,
    _parse_chars_in_mol_column,
    _parse_ring_size_column,
)
from hedgehog.descriptors.io import process_path


def _get_model_colors(model_names):
    """Get color mapping for models.

    Args:
        model_names: List of model name strings

    Returns:
        dict: Mapping of model names to colors
    """
    distinct_colors = [
        "brown",
        "green",
        "blue",
        "cyan",
        "yellow",
        "pink",
        "orange",
        "#dd37fa",
        "#ad5691",
        "#f46fa1",
        "#89cff0",
        "#93c83e",
    ]
    if len(model_names) <= 12:
        colors = distinct_colors
    else:
        colors = plt.cm.YlOrRd(np.linspace(1, 0, len(model_names) + 1))
    return dict(zip(sorted(model_names), colors, strict=False))


def _get_column_values(model_df, col):
    """Extract values from a column, handling special columns.

    Args:
        model_df: DataFrame for a specific model
        col: Column name

    Returns:
        list: Extracted values
    """
    if col == "chars":
        return _parse_chars_in_mol_column(model_df[col].dropna())
    if col == "ring_size":
        return _parse_ring_size_column(model_df[col].dropna())
    return model_df[col].dropna().tolist()


def _filter_values_by_bounds(values, min_val, max_val):
    """Filter values within min/max bounds.

    Args:
        values: List of values
        min_val: Minimum value (or None)
        max_val: Maximum value (or None, or 'inf')

    Returns:
        list: Filtered values
    """
    result = values.copy()
    if min_val is not None:
        result = [v for v in result if v >= min_val]
    if max_val is not None and max_val != "inf":
        result = [v for v in result if v <= max_val]
    return result


def _plot_discrete_chars(
    ax, values, offset, bar_width, color, label, borders, model_index
):
    """Plot discrete character counts as bar chart."""
    value_counts = pd.Series(values).value_counts()
    desired_order = ["C", "N", "S", "O", "F", "Cl", "Br", "H"]
    all_chars = borders.get("allowed_chars", desired_order)
    sorted_chars = [c for c in desired_order if c in all_chars] + [
        c for c in all_chars if c not in desired_order
    ]

    complete_counts = pd.Series(0, index=sorted_chars)
    complete_counts.update(value_counts)
    x_positions = [i + offset for i in range(len(complete_counts.index))]
    ax.bar(
        x_positions,
        complete_counts.values,
        width=bar_width,
        alpha=0.4,
        color=color,
        edgecolor="black",
        linewidth=0.3,
        label=label,
    )

    if model_index == 0:
        ax.set_xticks(list(range(len(complete_counts.index))))
        ax._discrete_tick_values = complete_counts.index
        ax.set_xticklabels(complete_counts.index)


def _plot_discrete_numeric(
    ax, values, col, offset, bar_width, color, label, max_val, model_index
):
    """Plot discrete numeric values as bar chart."""
    value_counts = pd.Series(values).value_counts().sort_index()

    if max_val is not None and max_val != "inf":
        extended_max = int(max_val) + 5
        full_range = list(range(0, extended_max + 1))
        complete_counts = pd.Series(0, index=full_range)
        for val in value_counts.index:
            if val in complete_counts.index:
                complete_counts[val] = value_counts[val]
    else:
        complete_counts = value_counts

    x_positions = [i + offset for i in range(len(complete_counts.index))]
    ax.bar(
        x_positions,
        complete_counts.values,
        width=bar_width,
        alpha=0.4,
        color=color,
        edgecolor="black",
        linewidth=0.3,
        label=label,
    )

    if model_index == 0:
        if col == "n_rigid_bonds":
            all_values = list(complete_counts.index)
            tick_values = [x for x in all_values if x % 5 == 0]
            if max(all_values) not in tick_values:
                tick_values.append(max(all_values))
            tick_positions = [
                all_values.index(val) for val in tick_values if val in all_values
            ]
            ax.set_xticks(tick_positions)
            ax.set_xticklabels([str(int(val)) for val in tick_values])
        else:
            ax.set_xticks(list(range(len(complete_counts.index))))
            ax.set_xticklabels([str(int(x)) for x in complete_counts.index])
        ax._discrete_tick_values = complete_counts.index


def _plot_continuous(ax, values, col, color, label):
    """Plot continuous distribution using KDE."""
    clip = (0, 1.0) if col in ["fsp3", "qed"] else (0, None)
    sns.kdeplot(
        values, label=label, fill=True, alpha=0.4, ax=ax, color=color, clip=clip
    )


def _get_boundary_position(ax, val, col, discrete_feats, fallback_offset=0):
    """Calculate boundary position for vertical lines and spans.

    Args:
        ax: Matplotlib axis
        val: Boundary value
        col: Column name
        discrete_feats: List of discrete feature names
        fallback_offset: Offset to add for continuous features

    Returns:
        float or None: Position for the boundary
    """
    if col not in discrete_feats:
        return val

    if not hasattr(ax, "_discrete_tick_values"):
        return val

    try:
        return list(ax._discrete_tick_values).index(val)
    except ValueError:
        return None


def _draw_boundary_lines(ax, col, min_val, max_val, discrete_feats):
    """Draw vertical lines at min/max boundaries."""
    # Small nudge offsets for continuous max boundaries that need
    # a tiny shift to avoid sitting exactly on the data edge.
    _MAX_NUDGE = {"fsp3": 0.01, "n_rigid_bonds": 0.0000001}

    boundaries = [
        (min_val, "red", "min", -0.5, 0),
        (max_val, "blue", "max", +0.5, None),
    ]
    for val, color, label_prefix, discrete_offset, _fallback_pos in boundaries:
        if val is None or val == "inf":
            continue

        pos = _get_boundary_position(ax, val, col, discrete_feats)
        if pos is not None:
            if col in discrete_feats:
                x = pos + discrete_offset
            elif label_prefix == "max" and col in _MAX_NUDGE:
                x = val + _MAX_NUDGE[col]
            else:
                x = val
        else:
            # Fallback position when value is not found among discrete ticks
            if label_prefix == "min":
                x = 0
            else:
                x = (
                    len(ax._discrete_tick_values) - 1
                    if hasattr(ax, "_discrete_tick_values")
                    else val
                )

        ax.axvline(
            x,
            color=color,
            linestyle="--",
            linewidth=1.5,
            label=f"{label_prefix}: {val}",
        )


def _draw_boundary_spans(ax, col, min_val, max_val, discrete_feats):
    """Draw shaded spans for excluded regions."""
    x_min, x_max = ax.get_xlim()

    # Draw min boundary span (shade area below min)
    if min_val is not None:
        pos = _get_boundary_position(ax, min_val, col, discrete_feats)
        if col in discrete_feats:
            if pos is not None:
                ax.axvspan(x_min, pos - 0.5, color="grey", alpha=0.2, zorder=0)
            else:
                ax.axvspan(x_min, 0, color="grey", alpha=0.2, zorder=0)
        else:
            ax.axvspan(x_min, min_val, color="grey", alpha=0.2, zorder=0)

    # Draw max boundary span (shade area above max)
    if max_val is not None and max_val != "inf":
        pos = _get_boundary_position(ax, max_val, col, discrete_feats)
        if col in discrete_feats:
            if pos is not None:
                ax.axvspan(pos + 0.5, x_max, color="grey", alpha=0.2, zorder=0)
            else:
                tick_len = (
                    len(ax._discrete_tick_values) - 0.5
                    if hasattr(ax, "_discrete_tick_values")
                    else x_max
                )
                ax.axvspan(tick_len, x_max, color="grey", alpha=0.2, zorder=0)
        else:
            ax.axvspan(max_val, x_max, color="grey", alpha=0.2, zorder=0)


def _style_axis(ax):
    """Apply consistent styling to axis."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.spines["left"].set_linewidth(0.5)
    ax.tick_params(
        axis="both",
        which="both",
        bottom=True,
        top=False,
        left=True,
        right=False,
        labelbottom=True,
        labeltop=False,
        labelleft=True,
        labelright=False,
        length=4,
        width=0.5,
        colors="black",
        labelsize=10,
    )


def _plot_single_column(
    ax, df, col, model_names, colors, borders, discrete_feats, renamer, is_multi
):
    """Plot distribution for a single column across all models.

    Args:
        ax: Matplotlib axis to plot on
        df: Full DataFrame
        col: Column name to plot
        model_names: List of model names
        colors: Color mapping for models
        borders: Border configurations
        discrete_feats: List of discrete feature names
        renamer: Column name mapping for display
        is_multi: Whether multiple models exist
    """
    min_val, max_val = _get_border_values(col, borders)
    minmax_str = f"min: {min_val}, max: {max_val}"
    bar_width = 0.1
    name_map = {
        str(m): str(m).upper() for m in sorted(df["model_name"].dropna().unique())
    }
    name_map_lc = {str(k).lower(): v for k, v in name_map.items()}

    for model_index, model in enumerate(model_names):
        model_df = df[df["model_name"].str.lower() == model] if is_multi else df
        values_before = _get_column_values(model_df, col)

        if not values_before:
            continue

        values_after = _filter_values_by_bounds(values_before, min_val, max_val)
        total = len(values_before)
        mols_passed = len(values_after) / total * 100 if total > 0 else 0

        label_name = name_map_lc.get(str(model).lower(), str(model).upper())
        label = f"{label_name}, pass: {mols_passed:.1f}%"
        color = colors[model]
        offset = model_index * bar_width

        if len(values_before) <= 1:
            ax.scatter(
                values_before,
                [0.01] * len(values_before),
                label=label,
                alpha=0.4,
                color=color,
            )
        elif col in discrete_feats:
            if col == "chars":
                _plot_discrete_chars(
                    ax,
                    values_before,
                    offset,
                    bar_width,
                    color,
                    label,
                    borders,
                    model_index,
                )
            else:
                _plot_discrete_numeric(
                    ax,
                    values_before,
                    col,
                    offset,
                    bar_width,
                    color,
                    label,
                    max_val,
                    model_index,
                )
        else:
            _plot_continuous(ax, values_before, col, color, label)

    # Set title and labels
    display_name = renamer.get(col, col)
    title = display_name if col == "chars" else f"{display_name} ({minmax_str})"
    ax.set_title(title, fontsize=12)
    ax.set_xlabel(display_name, fontsize=10)

    # Set tick locators
    if col not in discrete_feats:
        if col == "fsp3":
            ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        elif col != "qed":
            ax.xaxis.set_major_locator(plt.MaxNLocator(integer=False))
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

    # Draw boundaries and styling
    _draw_boundary_lines(ax, col, min_val, max_val, discrete_feats)
    _draw_boundary_spans(ax, col, min_val, max_val, discrete_feats)

    # Add sorted legend
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        sorted_pairs = sorted(zip(labels, handles, strict=False), key=lambda t: t[0])
        sorted_labels, sorted_handles = zip(*sorted_pairs, strict=False)
        ax.legend(sorted_handles, sorted_labels, fontsize=8, loc="upper right")

    _style_axis(ax)


def draw_filtered_mols(df, folder_to_save, config, progress_cb=None):
    """Generate distribution plots for descriptor filters.

    Args:
        df: DataFrame with computed descriptors
        folder_to_save: Output folder path (should already include 'Descriptors' subfolder)
        config: Configuration dictionary
    """
    folder_to_save = Path(process_path(folder_to_save))

    descriptors_config = load_config(config[CFG_DESCRIPTORS])
    normalized_borders, constraints = _extract_borders_and_constraints(
        descriptors_config.get("borders", {}),
        descriptors_config.get("structural_constraints"),
    )
    borders = dict(normalized_borders)
    borders.update(_build_structural_constraint_borders(constraints))
    if "charged_mol_allowed" in borders:
        borders["charged_mol_allowed"] = int(borders["charged_mol_allowed"])
    else:
        borders["charged_mol_allowed"] = False

    cols_to_plot = descriptors_config["filtered_cols_to_plot"]
    discrete_feats = descriptors_config["discrete_features_to_plot"]
    renamer = descriptors_config["renamer"]

    model_names = sorted(
        [m.lower() for m in df["model_name"].dropna().unique().tolist()]
    )
    is_multi = len(model_names) > 1
    colors = _get_model_colors(model_names)

    # Calculate grid dimensions
    n_cols = min(5, len(cols_to_plot))
    n_rows = (len(cols_to_plot) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(
        nrows=n_rows, ncols=n_cols, figsize=(n_rows * n_cols, n_rows * n_cols)
    )
    axes = axes.flatten()

    # Plot each column
    plotted_count = 0
    total_to_plot = max(len(cols_to_plot), 1)
    for i, col in enumerate(cols_to_plot):
        col_lower = col.lower()
        canonical = _DESCRIPTOR_KEY_MAP.get(col_lower, col).lower()
        relevant_keys = [
            k
            for k in borders
            if k.endswith(("_min", "_max"))
            and k.rsplit("_", 1)[0].lower() in (col_lower, canonical)
        ]
        if not relevant_keys:
            if progress_cb is not None:
                progress_cb(i + 1, total_to_plot)
            continue
        _plot_single_column(
            axes[i],
            df,
            col,
            model_names,
            colors,
            borders,
            discrete_feats,
            renamer,
            is_multi,
        )
        plotted_count += 1
        if progress_cb is not None:
            progress_cb(i + 1, total_to_plot)

    # Remove unused axes
    for j in range(len(cols_to_plot), len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.subplots_adjust(top=0.93, bottom=0.05, left=0.05, right=0.98)
    plt.savefig(
        folder_to_save / "descriptors_distribution.png",
        dpi=300,
        bbox_inches="tight",
        format="png",
    )
    if progress_cb is not None and plotted_count == 0:
        progress_cb(total_to_plot, total_to_plot)
