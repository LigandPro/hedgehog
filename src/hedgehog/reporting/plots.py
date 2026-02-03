"""Plot generation functions for HEDGEHOG reports using Plotly."""

import statistics
from typing import Any

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Default color palette for models
MODEL_COLORS = px.colors.qualitative.Set2


def plot_funnel(funnel_data: list[dict[str, Any]]) -> str:
    """Create a horizontal pipeline funnel chart showing molecule counts at each stage.

    Args:
        funnel_data: List of dicts with 'stage' and 'count' keys

    Returns:
        HTML string of the plotly figure
    """
    if not funnel_data:
        return _empty_plot("No funnel data available")

    stages = [d["stage"] for d in funnel_data]
    counts = [d["count"] for d in funnel_data]

    # Purple gradient for funnel segments
    n = len(stages)
    colors = [f"rgba(139, 92, 246, {0.95 - i * 0.1})" for i in range(n)]

    fig = go.Figure(
        go.Funnel(
            x=stages,
            y=counts,
            orientation="h",
            textposition="inside",
            textinfo="value+percent initial",
            textfont={"size": 12, "color": "white"},
            marker={"color": colors},
            connector={"line": {"color": "rgba(139, 92, 246, 0.3)", "width": 1}},
        )
    )

    fig.update_layout(
        font={"family": "-apple-system, BlinkMacSystemFont, sans-serif", "size": 11},
        height=320,
        margin={"l": 20, "r": 20, "t": 20, "b": 60},
        plot_bgcolor="white",
        paper_bgcolor="white",
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_sankey(funnel_data: list[dict[str, Any]]) -> str:
    """Create a Sankey diagram showing molecule flow through pipeline stages.

    Shows the flow of molecules through each stage with "Lost" nodes
    indicating molecules filtered out at each stage.

    Args:
        funnel_data: List of dicts with 'stage' and 'count' keys

    Returns:
        HTML string of the plotly figure
    """
    if not funnel_data or len(funnel_data) < 2:
        return _empty_plot("No funnel data available for Sankey diagram")

    stages = [d["stage"] for d in funnel_data]
    counts = [d["count"] for d in funnel_data]

    # Build Sankey nodes and links
    # Nodes: stage names + "Lost at <stage>" for each transition
    labels = list(stages)
    node_colors = []

    # Colors for main flow stages (soft purple palette)
    stage_color = "rgba(139, 92, 246, 0.85)"  # Soft purple
    lost_color = "rgba(156, 163, 175, 0.6)"  # Gray for lost molecules

    for _ in stages:
        node_colors.append(stage_color)

    # Add "Lost" nodes for each transition
    lost_nodes = []
    for i in range(len(stages) - 1):
        lost_count = counts[i] - counts[i + 1]
        if lost_count > 0:
            lost_label = f"Lost ({stages[i + 1]})"
            labels.append(lost_label)
            node_colors.append(lost_color)
            lost_nodes.append(
                {
                    "from_idx": i,
                    "to_idx": len(labels) - 1,
                    "count": lost_count,
                    "from_stage": stages[i],
                    "at_stage": stages[i + 1],
                }
            )

    # Build links
    sources = []
    targets = []
    values = []
    link_colors = []

    main_flow_color = "rgba(139, 92, 246, 0.35)"  # Soft purple
    lost_flow_color = "rgba(156, 163, 175, 0.35)"  # Gray for lost flow

    # Main flow links (between consecutive stages)
    for i in range(len(stages) - 1):
        sources.append(i)
        targets.append(i + 1)
        values.append(counts[i + 1])
        link_colors.append(main_flow_color)

    # Lost flow links
    for lost in lost_nodes:
        sources.append(lost["from_idx"])
        targets.append(lost["to_idx"])
        values.append(lost["count"])
        link_colors.append(lost_flow_color)

    # Calculate node positions
    n_stages = len(stages)
    n_lost = len(lost_nodes)

    # X positions: stages evenly spaced
    x_positions = []
    for i in range(n_stages):
        x_positions.append(i / max(n_stages - 1, 1))

    # Lost nodes positioned at the same x as target stage, but lower y
    for lost in lost_nodes:
        # Position lost node at the same x as the "at_stage" (target stage)
        target_idx = lost["from_idx"] + 1
        x_positions.append(target_idx / max(n_stages - 1, 1))

    # Y positions: main stages at top, lost nodes at bottom
    y_positions = []
    for _ in range(n_stages):
        y_positions.append(0.3)  # Main flow at top

    for _ in lost_nodes:
        y_positions.append(0.7)  # Lost nodes below (closer to main)

    # Hover text
    customdata = []
    for i, stage in enumerate(stages):
        if i == 0:
            customdata.append(f"{stage}: {counts[i]} molecules (100%)")
        else:
            pct = 100 * counts[i] / counts[0] if counts[0] > 0 else 0
            customdata.append(f"{stage}: {counts[i]} molecules ({pct:.1f}% of initial)")

    for lost in lost_nodes:
        pct = 100 * lost["count"] / counts[0] if counts[0] > 0 else 0
        customdata.append(
            f"Lost at {lost['at_stage']}: {lost['count']} molecules ({pct:.1f}% of initial)"
        )

    fig = go.Figure(
        data=[
            go.Sankey(
                arrangement="snap",
                node={
                    "pad": 20,
                    "thickness": 20,
                    "line": {"color": "black", "width": 0.5},
                    "label": labels,
                    "color": node_colors,
                    "x": x_positions,
                    "y": y_positions,
                    "customdata": customdata,
                    "hovertemplate": "%{customdata}<extra></extra>",
                },
                link={
                    "source": sources,
                    "target": targets,
                    "value": values,
                    "color": link_colors,
                    "hovertemplate": (
                        "From %{source.label}<br>"
                        "To %{target.label}<br>"
                        "Molecules: %{value}<extra></extra>"
                    ),
                },
            )
        ]
    )

    fig.update_layout(
        font={"size": 12},
        height=450,
        margin={"l": 20, "r": 20, "t": 20, "b": 20},
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_sankey_json(funnel_data: list[dict[str, Any]]) -> dict:
    """Generate Sankey diagram data as JSON for client-side rendering.

    This is used for dynamic model filtering where the chart needs
    to be redrawn with JavaScript.

    Args:
        funnel_data: List of dicts with 'stage' and 'count' keys

    Returns:
        Dictionary with nodes, links, and layout for Plotly.js
    """
    if not funnel_data or len(funnel_data) < 2:
        return {}

    stages = [d["stage"] for d in funnel_data]
    counts = [d["count"] for d in funnel_data]

    labels = list(stages)
    node_colors = []

    stage_color = "rgba(139, 92, 246, 0.85)"
    lost_color = "rgba(156, 163, 175, 0.6)"  # Gray for lost molecules

    for _ in stages:
        node_colors.append(stage_color)

    lost_nodes = []
    for i in range(len(stages) - 1):
        lost_count = counts[i] - counts[i + 1]
        if lost_count > 0:
            lost_label = f"Lost ({stages[i + 1]})"
            labels.append(lost_label)
            node_colors.append(lost_color)
            lost_nodes.append(
                {
                    "from_idx": i,
                    "to_idx": len(labels) - 1,
                    "count": lost_count,
                }
            )

    sources = []
    targets = []
    values = []
    link_colors = []

    main_flow_color = "rgba(139, 92, 246, 0.35)"
    lost_flow_color = "rgba(156, 163, 175, 0.35)"  # Gray for lost flow

    for i in range(len(stages) - 1):
        sources.append(i)
        targets.append(i + 1)
        values.append(counts[i + 1])
        link_colors.append(main_flow_color)

    for lost in lost_nodes:
        sources.append(lost["from_idx"])
        targets.append(lost["to_idx"])
        values.append(lost["count"])
        link_colors.append(lost_flow_color)

    n_stages = len(stages)
    x_positions = [i / max(n_stages - 1, 1) for i in range(n_stages)]
    for lost in lost_nodes:
        target_idx = lost["from_idx"] + 1
        x_positions.append(target_idx / max(n_stages - 1, 1))

    y_positions = [0.3] * n_stages + [0.7] * len(lost_nodes)

    return {
        "labels": labels,
        "node_colors": node_colors,
        "x_positions": x_positions,
        "y_positions": y_positions,
        "sources": sources,
        "targets": targets,
        "values": values,
        "link_colors": link_colors,
    }


def plot_sankey_compare_json(
    funnel_by_model: dict[str, list[dict[str, Any]]],
    models: list[str],
) -> dict:
    """Generate comparison Sankey data with each model in different color.

    Args:
        funnel_by_model: Dict mapping model names to funnel data
        models: List of model names

    Returns:
        Dictionary with Sankey data for Plotly.js
    """
    if not funnel_by_model or not models:
        return {}

    # Purple palette for models (different shades)
    model_colors_solid = [
        "rgba(139, 92, 246, 0.9)",  # Purple 500
        "rgba(167, 139, 250, 0.9)",  # Purple 400
        "rgba(196, 181, 253, 0.9)",  # Purple 300
        "rgba(109, 40, 217, 0.9)",  # Purple 700
        "rgba(124, 58, 237, 0.9)",  # Purple 600
        "rgba(221, 214, 254, 0.9)",  # Purple 200
        "rgba(91, 33, 182, 0.9)",  # Purple 800
        "rgba(76, 29, 149, 0.9)",  # Purple 900
    ]

    model_colors_light = [
        "rgba(139, 92, 246, 0.5)",
        "rgba(167, 139, 250, 0.5)",
        "rgba(196, 181, 253, 0.6)",
        "rgba(109, 40, 217, 0.5)",
        "rgba(124, 58, 237, 0.5)",
        "rgba(221, 214, 254, 0.7)",
        "rgba(91, 33, 182, 0.5)",
        "rgba(76, 29, 149, 0.5)",
    ]

    # Gray shades for lost molecules
    lost_colors = [
        "rgba(156, 163, 175, 0.3)",  # Gray 400
        "rgba(107, 114, 128, 0.3)",  # Gray 500
        "rgba(75, 85, 99, 0.3)",  # Gray 600
        "rgba(209, 213, 219, 0.4)",  # Gray 300
        "rgba(55, 65, 81, 0.3)",  # Gray 700
        "rgba(229, 231, 235, 0.5)",  # Gray 200
        "rgba(31, 41, 55, 0.3)",  # Gray 800
        "rgba(17, 24, 39, 0.3)",  # Gray 900
    ]

    first_model_data = funnel_by_model.get(models[0], [])
    if not first_model_data:
        return {}

    stages = [d["stage"] for d in first_model_data]
    n_stages = len(stages)

    # Nodes: stages + lost nodes for each stage transition
    labels = list(stages)
    node_colors = ["rgba(107, 114, 128, 0.7)"] * n_stages  # Gray for stage nodes

    # Add "Lost" nodes
    for i in range(n_stages - 1):
        labels.append(f"Lost ({stages[i + 1]})")
        node_colors.append("rgba(107, 114, 128, 0.3)")

    # Calculate positions
    x_positions = [i / max(n_stages - 1, 1) for i in range(n_stages)]
    y_positions = [0.3] * n_stages

    # Lost nodes at bottom (closer to main flow)
    for i in range(n_stages - 1):
        x_positions.append((i + 1) / max(n_stages - 1, 1))
        y_positions.append(0.7)

    sources = []
    targets = []
    values = []
    link_colors = []
    link_labels = []

    for model_idx, model in enumerate(models):
        model_funnel = funnel_by_model.get(model, [])
        if len(model_funnel) < 2:
            continue

        counts = [d["count"] for d in model_funnel]
        color_idx = model_idx % len(model_colors_light)

        # Main flow links
        for i in range(len(counts) - 1):
            if counts[i + 1] > 0:
                sources.append(i)
                targets.append(i + 1)
                values.append(counts[i + 1])
                link_colors.append(model_colors_light[color_idx])
                link_labels.append(model)

        # Lost flow links
        for i in range(len(counts) - 1):
            lost = counts[i] - counts[i + 1]
            if lost > 0:
                lost_node_idx = n_stages + i
                sources.append(i)
                targets.append(lost_node_idx)
                values.append(lost)
                link_colors.append(lost_colors[color_idx])
                link_labels.append(model)

    return {
        "labels": labels,
        "node_colors": node_colors,
        "x_positions": x_positions,
        "y_positions": y_positions,
        "sources": sources,
        "targets": targets,
        "values": values,
        "link_colors": link_colors,
        "link_labels": link_labels,
        "models": models,
        "model_colors": model_colors_solid[: len(models)],
        "is_comparison": True,
    }


def plot_model_comparison(model_stats: list[dict[str, Any]]) -> str:
    """Create grouped bar chart comparing models by retention rate.

    Args:
        model_stats: List of dicts with model statistics

    Returns:
        HTML string of the plotly figure
    """
    if not model_stats:
        return _empty_plot("No model data available")

    models = [m["model_name"] for m in model_stats]
    initial = [m.get("initial", 0) for m in model_stats]
    final = [m.get("final", 0) for m in model_stats]
    retention = [
        100.0 * f / i if i > 0 else 0 for i, f in zip(initial, final, strict=False)
    ]

    fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=("Molecule Counts", "Retention Rate (%)"),
        horizontal_spacing=0.15,
    )

    # Bar chart for counts
    fig.add_trace(
        go.Bar(
            name="Initial",
            x=models,
            y=initial,
            marker_color=px.colors.qualitative.Set2[0],
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Bar(
            name="Final",
            x=models,
            y=final,
            marker_color=px.colors.qualitative.Set2[1],
        ),
        row=1,
        col=1,
    )

    # Retention rate bar
    fig.add_trace(
        go.Bar(
            name="Retention %",
            x=models,
            y=retention,
            marker_color=px.colors.qualitative.Set2[2],
            text=[f"{r:.1f}%" for r in retention],
            textposition="outside",
            showlegend=False,
        ),
        row=1,
        col=2,
    )

    fig.update_layout(
        title="Model Comparison",
        barmode="group",
        height=400,
        showlegend=True,
        legend={
            "orientation": "h",
            "yanchor": "bottom",
            "y": 1.02,
            "xanchor": "right",
            "x": 1,
        },
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_model_stacked_losses(model_stats: list[dict[str, Any]]) -> str:
    """Create stacked bar chart showing where molecules were lost per model.

    Args:
        model_stats: List of dicts with stage-wise loss data

    Returns:
        HTML string of the plotly figure
    """
    if not model_stats:
        return _empty_plot("No model loss data available")

    models = [m["model_name"] for m in model_stats]
    stages = ["descriptors", "struct_filters", "synthesis", "docking"]
    stage_labels = ["Descriptors", "Structural Filters", "Synthesis", "Docking"]

    fig = go.Figure()

    for i, (stage, label) in enumerate(zip(stages, stage_labels, strict=False)):
        losses = [m.get("losses", {}).get(stage, 0) for m in model_stats]
        fig.add_trace(
            go.Bar(
                name=label,
                x=models,
                y=losses,
                marker_color=MODEL_COLORS[i % len(MODEL_COLORS)],
            )
        )

    fig.update_layout(
        title="Molecule Losses by Stage per Model",
        barmode="stack",
        height=400,
        xaxis_title="Model",
        yaxis_title="Molecules Lost",
        legend={
            "orientation": "h",
            "yanchor": "bottom",
            "y": 1.02,
            "xanchor": "right",
            "x": 1,
        },
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_descriptor_distributions(
    descriptor_data: dict[str, list[float]],
    model_names: list[str] | None = None,
) -> str:
    """Create violin/box plots for descriptor distributions.

    Args:
        descriptor_data: Dict mapping descriptor names to values
        model_names: Optional list of model names for grouping

    Returns:
        HTML string of the plotly figure
    """
    if not descriptor_data:
        return _empty_plot("No descriptor data available")

    # Create subplots for each descriptor
    descriptors = list(descriptor_data.keys())[:6]  # Limit to 6 descriptors
    n_cols = min(3, len(descriptors))
    n_rows = (len(descriptors) + n_cols - 1) // n_cols

    fig = make_subplots(
        rows=n_rows,
        cols=n_cols,
        subplot_titles=descriptors,
        vertical_spacing=0.15,
        horizontal_spacing=0.1,
    )

    for i, desc in enumerate(descriptors):
        row = i // n_cols + 1
        col = i % n_cols + 1

        values = descriptor_data[desc]
        fig.add_trace(
            go.Violin(
                y=values,
                name=desc,
                box_visible=True,
                meanline_visible=True,
                fillcolor=MODEL_COLORS[i % len(MODEL_COLORS)],
                line_color="black",
                showlegend=False,
            ),
            row=row,
            col=col,
        )

    fig.update_layout(
        title="Descriptor Distributions",
        height=300 * n_rows,
        showlegend=False,
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_descriptor_by_model(
    data: list[dict[str, Any]],
    descriptor: str,
) -> str:
    """Create box plot for a single descriptor grouped by model.

    Args:
        data: List of dicts with 'model_name' and descriptor value
        descriptor: Name of the descriptor to plot

    Returns:
        HTML string of the plotly figure
    """
    if not data:
        return _empty_plot(f"No data for {descriptor}")

    models = [d.get("model_name", "Unknown") for d in data]
    values = [d.get(descriptor, 0) for d in data]

    fig = px.box(
        x=models,
        y=values,
        color=models,
        labels={"x": "Model", "y": descriptor},
    )

    fig.update_layout(
        title=f"{descriptor} by Model",
        showlegend=False,
        height=350,
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_filter_heatmap(filter_data: dict[str, dict[str, int]]) -> str:
    """Create heatmap showing filter failures by model.

    Args:
        filter_data: Dict mapping filter names to {model: count} dicts

    Returns:
        HTML string of the plotly figure
    """
    if not filter_data:
        return _empty_plot("No filter data available")

    filters = list(filter_data.keys())
    if not filters:
        return _empty_plot("No filter data available")

    # Get all models
    models = set()
    for counts in filter_data.values():
        models.update(counts.keys())
    models = sorted(models)

    if not models:
        return _empty_plot("No model data in filters")

    # Build matrix
    z = []
    for filter_name in filters:
        row = [filter_data[filter_name].get(model, 0) for model in models]
        z.append(row)

    fig = go.Figure(
        data=go.Heatmap(
            z=z,
            x=models,
            y=filters,
            colorscale="RdYlGn_r",
            text=z,
            texttemplate="%{text}",
            textfont={"size": 12},
            hoverongaps=False,
        )
    )

    fig.update_layout(
        title="Filter Failures by Model",
        xaxis_title="Model",
        yaxis_title="Filter",
        height=max(300, 40 * len(filters)),
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_top_filter_failures(filter_totals: dict[str, int], top_n: int = 10) -> str:
    """Create bar chart of top filter failure reasons.

    Args:
        filter_totals: Dict mapping filter names to total failure counts
        top_n: Number of top filters to show

    Returns:
        HTML string of the plotly figure
    """
    if not filter_totals:
        return _empty_plot("No filter failure data available")

    # Sort and get top N
    sorted_filters = sorted(filter_totals.items(), key=lambda x: x[1], reverse=True)[
        :top_n
    ]
    filters, counts = zip(*sorted_filters, strict=False) if sorted_filters else ([], [])

    fig = go.Figure(
        go.Bar(
            x=list(counts),
            y=list(filters),
            orientation="h",
            marker_color=px.colors.sequential.Reds_r[: len(filters)],
            text=list(counts),
            textposition="outside",
        )
    )

    fig.update_layout(
        title=f"Top {top_n} Filter Failure Reasons",
        xaxis_title="Molecules Failed",
        yaxis_title="Filter",
        height=max(300, 35 * len(filters)),
        yaxis={"categoryorder": "total ascending"},
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_synthesis_distributions(synthesis_data: dict[str, list[float]]) -> str:
    """Create distribution plots for synthesis scores (SA, SYBA, RA).

    Args:
        synthesis_data: Dict mapping score names to lists of values

    Returns:
        HTML string of the plotly figure
    """
    if not synthesis_data:
        return _empty_plot("No synthesis score data available")

    score_names = list(synthesis_data.keys())
    n_scores = len(score_names)

    fig = make_subplots(
        rows=1,
        cols=n_scores,
        subplot_titles=score_names,
        horizontal_spacing=0.1,
    )

    for i, score_name in enumerate(score_names):
        values = synthesis_data[score_name]
        fig.add_trace(
            go.Histogram(
                x=values,
                name=score_name,
                marker_color=MODEL_COLORS[i % len(MODEL_COLORS)],
                opacity=0.75,
                showlegend=False,
            ),
            row=1,
            col=i + 1,
        )

    fig.update_layout(
        title="Synthesis Score Distributions",
        height=350,
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_synthesis_scatter(
    sa_scores: list[float],
    syba_scores: list[float],
    model_names: list[str],
) -> str:
    """Create scatter plot of SA Score vs SYBA Score colored by model.

    Args:
        sa_scores: List of SA score values
        syba_scores: List of SYBA score values
        model_names: List of model names for coloring

    Returns:
        HTML string of the plotly figure
    """
    if not sa_scores or not syba_scores:
        return _empty_plot("No synthesis scatter data available")

    fig = px.scatter(
        x=sa_scores,
        y=syba_scores,
        color=model_names,
        labels={"x": "SA Score", "y": "SYBA Score", "color": "Model"},
        opacity=0.7,
    )

    fig.update_layout(
        title="SA Score vs SYBA Score by Model",
        height=400,
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_docking_distribution(scores: list[float], tool: str = "gnina") -> str:
    """Create distribution plot for docking scores.

    Args:
        scores: List of docking score values
        tool: Name of docking tool (gnina, smina)

    Returns:
        HTML string of the plotly figure
    """
    if not scores:
        return _empty_plot("No docking score data available")

    fig = go.Figure()

    fig.add_trace(
        go.Histogram(
            x=scores,
            name=f"{tool.upper()} Scores",
            marker_color=px.colors.qualitative.Set2[4],
            opacity=0.75,
        )
    )

    fig.add_trace(
        go.Scatter(
            x=scores,
            y=[0] * len(scores),
            mode="markers",
            marker={"symbol": "line-ns-open", "size": 10, "color": "black"},
            showlegend=False,
        )
    )

    fig.update_layout(
        title=f"{tool.upper()} Docking Score Distribution",
        xaxis_title="Docking Score (kcal/mol)",
        yaxis_title="Count",
        height=350,
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_stage_summary(stage_stats: list[dict[str, Any]]) -> str:
    """Create summary bar chart of molecules passing each stage.

    Args:
        stage_stats: List of dicts with 'stage', 'passed', 'failed' keys

    Returns:
        HTML string of the plotly figure
    """
    if not stage_stats:
        return _empty_plot("No stage summary data available")

    stages = [s["stage"] for s in stage_stats]
    passed = [s.get("passed", 0) for s in stage_stats]
    failed = [s.get("failed", 0) for s in stage_stats]

    fig = go.Figure()

    fig.add_trace(
        go.Bar(
            name="Passed",
            x=stages,
            y=passed,
            marker_color="#2ecc71",
        )
    )

    fig.add_trace(
        go.Bar(
            name="Failed",
            x=stages,
            y=failed,
            marker_color="#e74c3c",
        )
    )

    fig.update_layout(
        title="Molecules Passed/Failed per Stage",
        barmode="stack",
        height=350,
        xaxis_title="Stage",
        yaxis_title="Molecules",
        legend={
            "orientation": "h",
            "yanchor": "bottom",
            "y": 1.02,
            "xanchor": "right",
            "x": 1,
        },
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def _empty_plot(message: str) -> str:
    """Create a placeholder figure with a message.

    Args:
        message: Message to display

    Returns:
        HTML string of the empty figure
    """
    fig = go.Figure()
    fig.add_annotation(
        text=message,
        xref="paper",
        yref="paper",
        x=0.5,
        y=0.5,
        showarrow=False,
        font={"size": 16, "color": "gray"},
    )
    fig.update_layout(
        xaxis={"visible": False},
        yaxis={"visible": False},
        height=200,
    )
    return fig.to_html(full_html=False, include_plotlyjs=False)


# Purple color palette for consistent styling (single model views)
PURPLE_PALETTE = [
    "rgba(139, 92, 246, 0.85)",  # Purple 500
    "rgba(167, 139, 250, 0.85)",  # Purple 400
    "rgba(196, 181, 253, 0.85)",  # Purple 300
    "rgba(109, 40, 217, 0.85)",  # Purple 700
    "rgba(124, 58, 237, 0.85)",  # Purple 600
    "rgba(221, 214, 254, 0.85)",  # Purple 200
    "rgba(91, 33, 182, 0.85)",  # Purple 800
    "rgba(76, 29, 149, 0.85)",  # Purple 900
]

# Distinct color palette for Compare mode (easily distinguishable)
COMPARE_PALETTE = [
    "rgba(59, 130, 246, 0.75)",  # Blue
    "rgba(239, 68, 68, 0.75)",  # Red
    "rgba(34, 197, 94, 0.75)",  # Green
    "rgba(249, 115, 22, 0.75)",  # Orange
    "rgba(168, 85, 247, 0.75)",  # Purple
    "rgba(236, 72, 153, 0.75)",  # Pink
    "rgba(20, 184, 166, 0.75)",  # Teal
    "rgba(234, 179, 8, 0.75)",  # Yellow
]


# =============================================================================
# DESCRIPTORS (Stage 01) - Enhanced plots
# =============================================================================


def plot_descriptors_violin_by_model(
    data: list[dict[str, Any]],
    descriptors: list[str],
) -> str:
    """Create violin plots for key descriptors grouped by model.

    Args:
        data: List of dicts with model_name and descriptor values
        descriptors: List of descriptor names to plot (e.g., MolWt, LogP, TPSA, QED)

    Returns:
        HTML string of the plotly figure
    """
    if not data or not descriptors:
        return _empty_plot("No descriptor data available")

    # Filter to available descriptors
    available = [d for d in descriptors if any(d in row for row in data)]
    if not available:
        return _empty_plot("No descriptor data available")

    n_desc = len(available)
    n_cols = min(2, n_desc)
    n_rows = (n_desc + n_cols - 1) // n_cols

    fig = make_subplots(
        rows=n_rows,
        cols=n_cols,
        subplot_titles=available,
        vertical_spacing=0.12,
        horizontal_spacing=0.1,
    )

    # Get unique models
    models = sorted(set(row.get("model_name", "Unknown") for row in data))

    for i, desc in enumerate(available):
        row = i // n_cols + 1
        col = i % n_cols + 1

        for j, model in enumerate(models):
            values = [
                row_data.get(desc)
                for row_data in data
                if row_data.get("model_name") == model
                and row_data.get(desc) is not None
            ]
            if values:
                fig.add_trace(
                    go.Violin(
                        y=values,
                        name=model,
                        legendgroup=model,
                        scalegroup=desc,
                        box_visible=True,
                        meanline_visible=True,
                        fillcolor=PURPLE_PALETTE[j % len(PURPLE_PALETTE)],
                        line_color="rgba(0,0,0,0.5)",
                        showlegend=(i == 0),
                    ),
                    row=row,
                    col=col,
                )

    fig.update_layout(
        height=350 * n_rows,
        showlegend=True,
        legend={
            "orientation": "h",
            "yanchor": "bottom",
            "y": 1.02,
            "xanchor": "right",
            "x": 1,
        },
        font={"family": "-apple-system, BlinkMacSystemFont, sans-serif", "size": 11},
        paper_bgcolor="white",
        plot_bgcolor="white",
        violinmode="group",
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_descriptors_hbd_hba_box(data: list[dict[str, Any]]) -> str:
    """Create box plots comparing HBD/HBA by model.

    Args:
        data: List of dicts with model_name, NumHDonors, NumHAcceptors

    Returns:
        HTML string of the plotly figure
    """
    if not data:
        return _empty_plot("No HBD/HBA data available")

    fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=["H-Bond Donors", "H-Bond Acceptors"],
        horizontal_spacing=0.12,
    )

    models = sorted(set(row.get("model_name", "Unknown") for row in data))

    for j, model in enumerate(models):
        # HBD
        hbd_values = [
            row.get("NumHDonors")
            for row in data
            if row.get("model_name") == model and row.get("NumHDonors") is not None
        ]
        if hbd_values:
            fig.add_trace(
                go.Box(
                    y=hbd_values,
                    name=model,
                    legendgroup=model,
                    marker_color=PURPLE_PALETTE[j % len(PURPLE_PALETTE)],
                    showlegend=True,
                ),
                row=1,
                col=1,
            )

        # HBA
        hba_values = [
            row.get("NumHAcceptors")
            for row in data
            if row.get("model_name") == model and row.get("NumHAcceptors") is not None
        ]
        if hba_values:
            fig.add_trace(
                go.Box(
                    y=hba_values,
                    name=model,
                    legendgroup=model,
                    marker_color=PURPLE_PALETTE[j % len(PURPLE_PALETTE)],
                    showlegend=False,
                ),
                row=1,
                col=2,
            )

    fig.update_layout(
        height=350,
        boxmode="group",
        showlegend=True,
        legend={
            "orientation": "h",
            "yanchor": "bottom",
            "y": 1.02,
            "xanchor": "right",
            "x": 1,
        },
        font={"family": "-apple-system, BlinkMacSystemFont, sans-serif", "size": 11},
        paper_bgcolor="white",
        plot_bgcolor="white",
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_descriptors_summary_table(summary: dict[str, dict[str, float]]) -> str:
    """Generate HTML table with descriptor means by model.

    Args:
        summary: Dict mapping model -> {descriptor: mean_value}

    Returns:
        HTML table string
    """
    if not summary:
        return "<p>No descriptor summary available</p>"

    models = sorted(summary.keys())
    if not models:
        return "<p>No descriptor summary available</p>"

    # Get all descriptors
    all_descs = set()
    for model_data in summary.values():
        all_descs.update(model_data.keys())
    descriptors = sorted(all_descs)

    rows = []
    for model in models:
        model_data = summary.get(model, {})
        cells = [
            f"<td>{model_data.get(d, '-'):.2f}</td>"
            if d in model_data
            else "<td>-</td>"
            for d in descriptors
        ]
        rows.append(f"<tr><td><strong>{model}</strong></td>{''.join(cells)}</tr>")

    header = "<th>Model</th>" + "".join(f"<th>{d}</th>" for d in descriptors)

    return f"""
    <div style="overflow-x: auto; max-width: 100%;">
        <table class="descriptor-table" style="min-width: 600px;">
            <thead><tr>{header}</tr></thead>
            <tbody>{"".join(rows)}</tbody>
        </table>
    </div>
    """


# =============================================================================
# STRUCTURAL FILTERS (Stage 02) - Enhanced plots
# =============================================================================


def plot_filter_stacked_bar(filter_data: dict[str, dict[str, int]]) -> str:
    """Create stacked bar chart showing molecules rejected per filter by model.

    Args:
        filter_data: Dict mapping filter_name -> {model: rejected_count}

    Returns:
        HTML string of the plotly figure
    """
    if not filter_data:
        return _empty_plot("No filter data available")

    # Get all models
    models = set()
    for counts in filter_data.values():
        models.update(counts.keys())
    models = sorted(models)

    if not models:
        return _empty_plot("No model data in filters")

    filters = list(filter_data.keys())

    fig = go.Figure()

    for i, filter_name in enumerate(filters):
        counts = [filter_data[filter_name].get(model, 0) for model in models]
        fig.add_trace(
            go.Bar(
                name=filter_name,
                x=models,
                y=counts,
                marker_color=PURPLE_PALETTE[i % len(PURPLE_PALETTE)],
            )
        )

    fig.update_layout(
        barmode="stack",
        height=400,
        xaxis_title="Model",
        yaxis_title="Molecules Rejected",
        legend={
            "orientation": "h",
            "yanchor": "bottom",
            "y": 1.02,
            "xanchor": "right",
            "x": 1,
        },
        font={"family": "-apple-system, BlinkMacSystemFont, sans-serif", "size": 11},
        paper_bgcolor="white",
        plot_bgcolor="white",
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_filter_banned_ratio_heatmap(
    ratio_data: dict[str, dict[str, float]],
) -> str:
    """Create heatmap showing banned_ratio for each filter × model.

    Args:
        ratio_data: Dict mapping filter_name -> {model: banned_ratio}

    Returns:
        HTML string of the plotly figure
    """
    if not ratio_data:
        return _empty_plot("No filter ratio data available")

    filters = list(ratio_data.keys())
    models = set()
    for ratios in ratio_data.values():
        models.update(ratios.keys())
    models = sorted(models)

    if not models or not filters:
        return _empty_plot("No filter ratio data available")

    # Build matrix
    z = []
    text = []
    for filter_name in filters:
        row = [ratio_data[filter_name].get(model, 0) for model in models]
        z.append(row)
        text.append([f"{v:.1%}" for v in row])

    fig = go.Figure(
        data=go.Heatmap(
            z=z,
            x=models,
            y=filters,
            colorscale=[
                [0, "rgba(255, 255, 255, 1)"],
                [0.5, "rgba(196, 181, 253, 1)"],
                [1, "rgba(109, 40, 217, 1)"],
            ],
            text=text,
            texttemplate="%{text}",
            textfont={"size": 11},
            hoverongaps=False,
            colorbar={"title": "Banned Ratio"},
        )
    )

    fig.update_layout(
        xaxis_title="Model",
        yaxis_title="Filter",
        height=max(300, 45 * len(filters)),
        font={"family": "-apple-system, BlinkMacSystemFont, sans-serif", "size": 11},
        paper_bgcolor="white",
        plot_bgcolor="white",
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_filter_top_reasons_bar(reasons_data: dict[str, int], top_n: int = 10) -> str:
    """Create horizontal bar chart of top rejection reasons.

    Args:
        reasons_data: Dict mapping reason -> count
        top_n: Number of top reasons to show

    Returns:
        HTML string of the plotly figure
    """
    if not reasons_data:
        return _empty_plot("No rejection reasons data available")

    # Sort and get top N
    sorted_reasons = sorted(reasons_data.items(), key=lambda x: x[1], reverse=True)[
        :top_n
    ]
    if not sorted_reasons:
        return _empty_plot("No rejection reasons data available")

    reasons, counts = zip(*sorted_reasons, strict=False)

    # Create gradient colors
    n = len(reasons)
    colors = [f"rgba(139, 92, 246, {0.95 - i * 0.05})" for i in range(n)]

    fig = go.Figure(
        go.Bar(
            x=list(counts),
            y=list(reasons),
            orientation="h",
            marker_color=colors,
            text=list(counts),
            textposition="outside",
        )
    )

    fig.update_layout(
        xaxis_title="Molecules Rejected",
        height=max(300, 35 * len(reasons)),
        yaxis={"categoryorder": "total ascending"},
        font={"family": "-apple-system, BlinkMacSystemFont, sans-serif", "size": 11},
        paper_bgcolor="white",
        plot_bgcolor="white",
        margin={"l": 150},
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


# =============================================================================
# SYNTHESIS (Stage 04) - Enhanced plots
# =============================================================================


def plot_synthesis_sa_histogram(scores: list[float]) -> str:
    """Create histogram of SA (Synthetic Accessibility) scores.

    Args:
        scores: List of SA score values (1-10, lower is easier)

    Returns:
        HTML string of the plotly figure
    """
    if not scores:
        return _empty_plot("No SA score data available")

    fig = go.Figure()

    fig.add_trace(
        go.Histogram(
            x=scores,
            nbinsx=30,
            marker_color="rgba(139, 92, 246, 0.75)",
            marker_line_color="rgba(109, 40, 217, 1)",
            marker_line_width=1,
        )
    )

    # Add mean line
    mean_val = statistics.mean(scores)
    fig.add_vline(
        x=mean_val,
        line_dash="dash",
        line_color="rgba(109, 40, 217, 1)",
        annotation_text=f"Mean: {mean_val:.2f}",
        annotation_position="top",
    )

    fig.update_layout(
        xaxis_title="SA Score (lower = easier to synthesize)",
        yaxis_title="Count",
        height=350,
        font={"family": "-apple-system, BlinkMacSystemFont, sans-serif", "size": 11},
        paper_bgcolor="white",
        plot_bgcolor="white",
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_synthesis_syba_histogram(scores: list[float]) -> str:
    """Create histogram of SYBA (SYnthetic Bayesian Accessibility) scores.

    Args:
        scores: List of SYBA score values (higher is easier)

    Returns:
        HTML string of the plotly figure
    """
    if not scores:
        return _empty_plot("No SYBA score data available")

    fig = go.Figure()

    fig.add_trace(
        go.Histogram(
            x=scores,
            nbinsx=30,
            marker_color="rgba(167, 139, 250, 0.75)",
            marker_line_color="rgba(124, 58, 237, 1)",
            marker_line_width=1,
        )
    )

    # Add mean line
    mean_val = statistics.mean(scores)
    fig.add_vline(
        x=mean_val,
        line_dash="dash",
        line_color="rgba(124, 58, 237, 1)",
        annotation_text=f"Mean: {mean_val:.2f}",
        annotation_position="top",
    )

    fig.update_layout(
        xaxis_title="SYBA Score (higher = easier to synthesize)",
        yaxis_title="Count",
        height=350,
        font={"family": "-apple-system, BlinkMacSystemFont, sans-serif", "size": 11},
        paper_bgcolor="white",
        plot_bgcolor="white",
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_synthesis_solved_pie(solved_count: int, unsolved_count: int) -> str:
    """Create pie/donut chart showing % of molecules with synthesis route found.

    Args:
        solved_count: Number of molecules with synthesis route
        unsolved_count: Number of molecules without synthesis route

    Returns:
        HTML string of the plotly figure
    """
    if solved_count == 0 and unsolved_count == 0:
        return _empty_plot("No synthesis route data available")

    fig = go.Figure(
        data=[
            go.Pie(
                labels=["Route Found", "No Route"],
                values=[solved_count, unsolved_count],
                hole=0.5,
                marker_colors=["rgba(139, 92, 246, 0.85)", "rgba(229, 231, 235, 0.85)"],
                textinfo="percent+label",
                textfont={"size": 12},
                hoverinfo="label+value+percent",
            )
        ]
    )

    total = solved_count + unsolved_count
    pct = 100 * solved_count / total if total > 0 else 0

    fig.update_layout(
        height=350,
        font={"family": "-apple-system, BlinkMacSystemFont, sans-serif", "size": 11},
        paper_bgcolor="white",
        annotations=[
            {
                "text": f"{pct:.1f}%<br>Solved",
                "x": 0.5,
                "y": 0.5,
                "font_size": 16,
                "showarrow": False,
            }
        ],
        showlegend=False,
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_synthesis_time_box(data: list[dict[str, Any]]) -> str:
    """Create box plot of search_time by model.

    Args:
        data: List of dicts with model_name and search_time

    Returns:
        HTML string of the plotly figure
    """
    if not data:
        return _empty_plot("No synthesis time data available")

    models = sorted(set(row.get("model_name", "Unknown") for row in data))

    fig = go.Figure()

    for i, model in enumerate(models):
        times = [
            row.get("search_time")
            for row in data
            if row.get("model_name") == model and row.get("search_time") is not None
        ]
        if times:
            fig.add_trace(
                go.Box(
                    y=times,
                    name=model,
                    marker_color=PURPLE_PALETTE[i % len(PURPLE_PALETTE)],
                )
            )

    fig.update_layout(
        yaxis_title="Search Time (s)",
        height=350,
        showlegend=False,
        font={"family": "-apple-system, BlinkMacSystemFont, sans-serif", "size": 11},
        paper_bgcolor="white",
        plot_bgcolor="white",
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


# =============================================================================
# DOCKING (Stage 05) - Enhanced plots
# =============================================================================


def plot_docking_affinity_histogram(scores: list[float], tool: str = "gnina") -> str:
    """Create histogram of docking binding affinity scores.

    Args:
        scores: List of binding affinity values (kcal/mol, more negative = better)
        tool: Docking tool name

    Returns:
        HTML string of the plotly figure
    """
    if not scores:
        return _empty_plot("No docking affinity data available")

    fig = go.Figure()

    fig.add_trace(
        go.Histogram(
            x=scores,
            nbinsx=30,
            marker_color="rgba(139, 92, 246, 0.75)",
            marker_line_color="rgba(109, 40, 217, 1)",
            marker_line_width=1,
        )
    )

    # Add mean and best lines
    mean_val = statistics.mean(scores)
    best_val = min(scores)

    fig.add_vline(
        x=mean_val,
        line_dash="dash",
        line_color="rgba(109, 40, 217, 0.8)",
        annotation_text=f"Mean: {mean_val:.2f}",
        annotation_position="top left",
    )

    fig.add_vline(
        x=best_val,
        line_dash="solid",
        line_color="rgba(76, 29, 149, 1)",
        annotation_text=f"Best: {best_val:.2f}",
        annotation_position="top right",
    )

    fig.update_layout(
        title=f"{tool.upper()} Binding Affinity Distribution",
        xaxis_title="Binding Affinity (kcal/mol)",
        yaxis_title="Count",
        height=350,
        font={"family": "-apple-system, BlinkMacSystemFont, sans-serif", "size": 11},
        paper_bgcolor="white",
        plot_bgcolor="white",
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_docking_top_molecules(
    data: list[dict[str, Any]],
    top_n: int = 10,
    score_col: str = "affinity",
) -> str:
    """Create bar chart of top molecules by binding affinity.

    Args:
        data: List of dicts with molecule_id/smiles and affinity score
        top_n: Number of top molecules to show
        score_col: Column name for the score

    Returns:
        HTML string of the plotly figure
    """
    if not data:
        return _empty_plot("No docking data available")

    # Sort by score (more negative = better)
    valid_data = [d for d in data if d.get(score_col) is not None]
    sorted_data = sorted(valid_data, key=lambda x: x.get(score_col, 0))[:top_n]

    if not sorted_data:
        return _empty_plot("No docking data available")

    # Get molecule identifiers
    labels = []
    for d in sorted_data:
        mol_id = d.get("molecule_id") or d.get("mol_id") or d.get("smiles", "Unknown")
        if len(str(mol_id)) > 20:
            mol_id = str(mol_id)[:17] + "..."
        labels.append(str(mol_id))

    scores = [d.get(score_col, 0) for d in sorted_data]

    # Gradient colors
    n = len(scores)
    colors = [f"rgba(139, 92, 246, {0.95 - i * 0.05})" for i in range(n)]

    fig = go.Figure(
        go.Bar(
            x=scores,
            y=labels,
            orientation="h",
            marker_color=colors,
            text=[f"{s:.2f}" for s in scores],
            textposition="outside",
        )
    )

    fig.update_layout(
        title=f"Top {top_n} Molecules by Binding Affinity",
        xaxis_title="Binding Affinity (kcal/mol)",
        height=max(300, 35 * len(labels)),
        yaxis={"categoryorder": "total descending"},
        font={"family": "-apple-system, BlinkMacSystemFont, sans-serif", "size": 11},
        paper_bgcolor="white",
        plot_bgcolor="white",
        margin={"l": 150},
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_docking_affinity_box(
    data: list[dict[str, Any]], score_col: str = "affinity"
) -> str:
    """Create box plot of binding affinity by model.

    Args:
        data: List of dicts with model_name and affinity score
        score_col: Column name for the score

    Returns:
        HTML string of the plotly figure
    """
    if not data:
        return _empty_plot("No docking data available")

    models = sorted(set(row.get("model_name", "Unknown") for row in data))

    fig = go.Figure()

    for i, model in enumerate(models):
        scores = [
            row.get(score_col)
            for row in data
            if row.get("model_name") == model and row.get(score_col) is not None
        ]
        if scores:
            fig.add_trace(
                go.Box(
                    y=scores,
                    name=model,
                    marker_color=PURPLE_PALETTE[i % len(PURPLE_PALETTE)],
                )
            )

    fig.update_layout(
        title="Binding Affinity by Model",
        yaxis_title="Binding Affinity (kcal/mol)",
        height=400,
        showlegend=False,
        font={"family": "-apple-system, BlinkMacSystemFont, sans-serif", "size": 11},
        paper_bgcolor="white",
        plot_bgcolor="white",
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


# =============================================================================
# Retrosynthesis (AiZynthFinder) Plots
# =============================================================================


def plot_retrosynthesis_route_score_histogram(scores: list[float]) -> str:
    """Create histogram of retrosynthesis route scores.

    Args:
        scores: List of route score values (0-1, higher is better)

    Returns:
        HTML string of the plotly figure
    """
    if not scores:
        return _empty_plot("No route score data available")

    import statistics

    mean_val = statistics.mean(scores)

    fig = go.Figure()

    fig.add_trace(
        go.Histogram(
            x=scores,
            nbinsx=20,
            marker_color="rgba(139, 92, 246, 0.75)",
            marker_line_color="rgba(109, 40, 217, 1)",
            marker_line_width=1,
        )
    )

    fig.add_vline(
        x=mean_val,
        line_dash="dash",
        line_color="rgba(109, 40, 217, 1)",
        annotation_text=f"Mean: {mean_val:.2f}",
        annotation_position="top",
    )

    fig.update_layout(
        xaxis_title="Route Score (higher = better)",
        yaxis_title="Count",
        height=300,
        margin={"l": 50, "r": 30, "t": 30, "b": 50},
        paper_bgcolor="white",
        plot_bgcolor="white",
        font={"family": "-apple-system, BlinkMacSystemFont, sans-serif", "size": 11},
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_retrosynthesis_steps_histogram(steps: list[int]) -> str:
    """Create bar chart of synthesis steps count.

    Args:
        steps: List of number of synthesis steps

    Returns:
        HTML string of the plotly figure
    """
    if not steps:
        return _empty_plot("No synthesis steps data available")

    from collections import Counter

    import statistics

    step_counts = Counter(steps)
    step_values = sorted(step_counts.keys())
    counts = [step_counts[s] for s in step_values]
    mean_val = statistics.mean(steps)

    fig = go.Figure()

    fig.add_trace(
        go.Bar(
            x=step_values,
            y=counts,
            marker_color="rgba(167, 139, 250, 0.75)",
            marker_line_color="rgba(124, 58, 237, 1)",
            marker_line_width=1,
            text=counts,
            textposition="outside",
        )
    )

    fig.add_vline(
        x=mean_val,
        line_dash="dash",
        line_color="rgba(109, 40, 217, 1)",
        annotation_text=f"Avg: {mean_val:.1f}",
        annotation_position="top",
    )

    fig.update_layout(
        xaxis_title="Number of Synthesis Steps",
        yaxis_title="Count",
        height=300,
        margin={"l": 50, "r": 30, "t": 30, "b": 50},
        paper_bgcolor="white",
        plot_bgcolor="white",
        font={"family": "-apple-system, BlinkMacSystemFont, sans-serif", "size": 11},
        xaxis={"dtick": 1},
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)
