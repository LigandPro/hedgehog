"""Plot generation functions for HEDGEHOG reports using Plotly."""

from typing import Any

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Default color palette for models
MODEL_COLORS = px.colors.qualitative.Set2


def plot_funnel(funnel_data: list[dict[str, Any]]) -> str:
    """Create a pipeline funnel chart showing molecule counts at each stage.

    Args:
        funnel_data: List of dicts with 'stage' and 'count' keys

    Returns:
        HTML string of the plotly figure
    """
    if not funnel_data:
        return _empty_plot("No funnel data available")

    stages = [d["stage"] for d in funnel_data]
    counts = [d["count"] for d in funnel_data]

    fig = go.Figure(
        go.Funnel(
            y=stages,
            x=counts,
            textposition="inside",
            textinfo="value+percent initial",
            marker={"color": px.colors.sequential.Teal},
            connector={"line": {"color": "royalblue", "dash": "dot", "width": 2}},
        )
    )

    fig.update_layout(
        title="Pipeline Molecule Funnel",
        font={"size": 14},
        height=400,
        margin={"l": 150, "r": 50, "t": 80, "b": 50},
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


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
        legend={"orientation": "h", "yanchor": "bottom", "y": 1.02, "xanchor": "right", "x": 1},
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
        legend={"orientation": "h", "yanchor": "bottom", "y": 1.02, "xanchor": "right", "x": 1},
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
        legend={"orientation": "h", "yanchor": "bottom", "y": 1.02, "xanchor": "right", "x": 1},
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
