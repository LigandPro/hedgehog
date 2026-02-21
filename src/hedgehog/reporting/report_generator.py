"""Report generator for HEDGEHOG pipeline results."""

import html as html_lib
import json
import logging
from collections.abc import Callable
from datetime import datetime
from pathlib import Path
from typing import Any

from jinja2 import Environment, PackageLoader

from hedgehog.reporting import data_collector, js_data, plots

logger = logging.getLogger(__name__)

_PLOT_REGISTRY: dict[str, Callable[[dict[str, Any]], str]] = {
    "funnel": lambda report_data: plots.plot_funnel(report_data.get("funnel", [])),
    "sankey": lambda report_data: plots.plot_sankey(report_data.get("funnel", [])),
    "stage_summary": lambda report_data: plots.plot_stage_summary(
        report_data.get("stages", [])
    ),
    "model_comparison": lambda report_data: (
        plots.plot_model_comparison(report_data.get("models", []))
        if report_data.get("models")
        else ""
    ),
    "model_losses": lambda report_data: (
        plots.plot_model_stacked_losses(report_data.get("models", []))
        if report_data.get("models")
        else ""
    ),
    "descriptors": lambda report_data: (
        plots.plot_descriptor_distributions(
            report_data["descriptors"]["distributions"],
        )
        if report_data.get("descriptors", {}).get("distributions")
        else ""
    ),
    "descriptors_violin": lambda report_data: (
        plots.plot_descriptors_violin_by_model(
            report_data["descriptors_detailed"]["raw_data"],
            report_data["descriptors_detailed"].get(
                "key_descriptors",
                ["MolWt", "LogP", "TPSA", "QED"],
            ),
        )
        if report_data.get("descriptors_detailed", {}).get("raw_data")
        else ""
    ),
    "descriptors_table": lambda report_data: (
        plots.plot_descriptors_summary_table(
            report_data["descriptors_detailed"]["summary_by_model"],
        )
        if report_data.get("descriptors_detailed", {}).get("summary_by_model")
        else ""
    ),
    "filter_heatmap": lambda report_data: (
        plots.plot_filter_heatmap(report_data["filters"]["by_filter"])
        if report_data.get("filters", {}).get("by_filter")
        else ""
    ),
    "filter_failures": lambda report_data: (
        plots.plot_top_filter_failures(report_data["filters"]["totals"])
        if report_data.get("filters", {}).get("totals")
        else ""
    ),
    "synthesis_dist": lambda report_data: (
        plots.plot_synthesis_distributions(report_data["synthesis"]["distributions"])
        if report_data.get("synthesis", {}).get("distributions")
        else ""
    ),
    "synthesis_scatter": lambda report_data: (
        plots.plot_synthesis_scatter(
            report_data["synthesis"]["scatter_data"]["sa_scores"],
            report_data["synthesis"]["scatter_data"]["syba_scores"],
            report_data["synthesis"]["scatter_data"].get("model_names", []),
        )
        if report_data.get("synthesis", {}).get("scatter_data", {}).get("sa_scores")
        and report_data.get("synthesis", {}).get("scatter_data", {}).get("syba_scores")
        else ""
    ),
    "synthesis_sa_hist": lambda report_data: (
        plots.plot_synthesis_sa_histogram(
            report_data["synthesis_detailed"]["sa_scores"]
        )
        if report_data.get("synthesis_detailed", {}).get("sa_scores")
        else ""
    ),
    "synthesis_syba_hist": lambda report_data: (
        plots.plot_synthesis_syba_histogram(
            report_data["synthesis_detailed"]["syba_scores"]
        )
        if report_data.get("synthesis_detailed", {}).get("syba_scores")
        else ""
    ),
    "synthesis_ra_hist": lambda report_data: (
        plots.plot_synthesis_ra_histogram(
            report_data["synthesis_detailed"]["ra_scores"]
        )
        if report_data.get("synthesis_detailed", {}).get("ra_scores")
        else ""
    ),
    "synthesis_solved_pie": lambda report_data: (
        plots.plot_synthesis_solved_pie(
            report_data["synthesis_detailed"].get("solved_count", 0),
            report_data["synthesis_detailed"].get("unsolved_count", 0),
        )
        if report_data.get("synthesis_detailed", {}).get("solved_count", 0) > 0
        or report_data.get("synthesis_detailed", {}).get("unsolved_count", 0) > 0
        else ""
    ),
    "synthesis_time_box": lambda report_data: (
        plots.plot_synthesis_time_box(report_data["synthesis_detailed"]["raw_data"])
        if report_data.get("synthesis_detailed", {}).get("raw_data")
        else ""
    ),
    "retrosynthesis_route_score_hist": lambda report_data: (
        plots.plot_retrosynthesis_route_score_histogram(
            report_data["retrosynthesis"]["route_scores"]
        )
        if report_data.get("retrosynthesis", {}).get("route_scores")
        else ""
    ),
    "retrosynthesis_steps_hist": lambda report_data: (
        plots.plot_retrosynthesis_steps_histogram(
            report_data["retrosynthesis"]["steps"]
        )
        if report_data.get("retrosynthesis", {}).get("steps")
        else ""
    ),
    "docking_gnina": lambda report_data: (
        plots.plot_docking_distribution(
            report_data["docking"]["gnina"]["scores"], "gnina"
        )
        if report_data.get("docking", {}).get("gnina", {}).get("scores")
        else ""
    ),
    "docking_smina": lambda report_data: (
        plots.plot_docking_distribution(
            report_data["docking"]["smina"]["scores"], "smina"
        )
        if report_data.get("docking", {}).get("smina", {}).get("scores")
        else ""
    ),
}


class ReportGenerator:
    """Generates comprehensive HTML reports for HEDGEHOG pipeline runs."""

    def __init__(
        self,
        base_path: Path,
        stages: list[Any],
        config: dict[str, Any],
        initial_count: int,
        final_count: int,
    ):
        """Initialize the report generator.

        Args:
            base_path: Base path of the pipeline run output
            stages: List of PipelineStage objects
            config: Pipeline configuration dictionary
            initial_count: Number of molecules at pipeline start
            final_count: Number of molecules at pipeline end
        """
        self.base_path = Path(base_path)
        self.stages = stages
        self.config = config
        self.initial_count = initial_count
        self.final_count = final_count
        self.output_dir = self.base_path
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def generate(self) -> Path:
        """Generate the full HTML report.

        Returns:
            Path to the generated report file
        """
        logger.info("Generating pipeline report...")

        # Collect all data
        data = self._collect_data()

        # Generate plots
        plot_htmls = self._generate_plots(data)

        # Render HTML
        html_content = self._render_template(data, plot_htmls)

        # Save report
        report_path = self._save_report(html_content)

        # Save JSON data
        self._save_json_data(data)

        # Generate/update RUN_INFO.md with MolEval metrics table
        self._generate_run_info(data)

        logger.info("Report generated: %s", report_path)
        return report_path

    def _collect_data(self) -> dict[str, Any]:
        """Collect all metrics from stage outputs.

        Returns:
            Dictionary with all collected data
        """
        return data_collector.collect_all_data(
            self.base_path,
            self.config,
            self.stages,
            self.initial_count,
            self.final_count,
            self.output_dir,
        )

    def _collect_stage_smiles(self) -> dict[str, list[str]]:
        """Collect SMILES from key pipeline stages for MolEval analysis."""
        return data_collector.collect_stage_smiles(
            self.base_path,
            self.config,
            self.stages,
            self.initial_count,
            self.final_count,
            self.output_dir,
        )

    def _read_stage_smiles(self, path) -> list[str]:
        """Read SMILES column from a stage CSV file."""
        return data_collector.read_stage_smiles(
            self.base_path,
            self.config,
            self.stages,
            self.initial_count,
            self.final_count,
            self.output_dir,
            path,
        )

    def _get_moleval_metrics(self) -> dict[str, Any]:
        """Compute MolEval generative metrics across pipeline stages."""
        return data_collector.get_moleval_metrics(
            self.base_path,
            self.config,
            self.stages,
            self.initial_count,
            self.final_count,
            self.output_dir,
        )

    def _generate_run_info(self, data: dict[str, Any]) -> None:
        """Append MolEval metrics table to RUN_INFO.md if metrics are available.

        Reads the existing RUN_INFO.md (generated by pipeline.py) and appends
        a formatted MolEval metrics table.  If RUN_INFO.md does not exist yet,
        creates a minimal version with summary + metrics.

        Args:
            data: Collected report data containing 'moleval' key.
        """
        try:
            run_info_path = self.base_path / "RUN_INFO.md"
            moleval = data.get("moleval", {})
            by_stage = moleval.get("by_stage", {})
            stages = moleval.get("stages", [])
            metrics = moleval.get("metrics", [])

            if not by_stage or not stages or not metrics:
                return

            # Build markdown table
            header = "| Metric | " + " | ".join(stages) + " |"
            separator = "|--------|" + "|".join("--------" for _ in stages) + "|"
            rows = []
            for metric in metrics:
                cells = []
                for stage in stages:
                    val = by_stage.get(stage, {}).get(metric)
                    cells.append(f"{val:.4f}" if val is not None else "—")
                rows.append(f"| {metric} | " + " | ".join(cells) + " |")

            table_md = "\n".join([header, separator, *rows])

            section = f"\n\n## MolEval Metrics\n\n{table_md}\n"

            if run_info_path.exists():
                existing = run_info_path.read_text()
                # Replace existing section if present, otherwise append
                marker = "## MolEval Metrics"
                if marker in existing:
                    # Remove old section (everything from marker to next ## or EOF)
                    before = existing[: existing.index(marker)]
                    rest = existing[existing.index(marker) + len(marker) :]
                    # Find next section header
                    next_section = rest.find("\n## ")
                    if next_section != -1:
                        after = rest[next_section:]
                    else:
                        after = ""
                    content = before.rstrip() + section + after
                else:
                    content = existing.rstrip() + section
            else:
                # Create minimal RUN_INFO.md
                retention = (
                    f"{100 * self.final_count / self.initial_count:.2f}%"
                    if self.initial_count > 0
                    else "N/A"
                )
                content = (
                    f"# HEDGEHOG Run Info\n\n"
                    f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"
                    f"## Summary\n\n"
                    f"- Initial molecules: {self.initial_count}\n"
                    f"- Final molecules: {self.final_count}\n"
                    f"- Retention rate: {retention}\n"
                    f"{section}"
                )

            with open(run_info_path, "w") as f:
                f.write(content)
            logger.info("Updated run info with MolEval metrics: %s", run_info_path)
        except (OSError, ValueError, ZeroDivisionError) as e:
            logger.debug("Failed to update RUN_INFO.md with MolEval metrics: %s", e)

    def _build_descriptors_js_data(
        self, desc_detailed: dict[str, Any], available_models: list[str]
    ) -> dict[str, Any]:
        """Build JSON data for JavaScript descriptors visualization."""
        return js_data.build_descriptors_js_data(desc_detailed, available_models)

    def _build_filters_js_data(self, available_models: list[str]) -> dict[str, Any]:
        """Build JSON data for JavaScript filters visualization."""
        return js_data.build_filters_js_data(self.base_path, available_models)

    def _build_descriptor_comparison_data(
        self,
        initial_detailed: dict[str, Any],
        final_detailed: dict[str, Any],
    ) -> dict[str, Any]:
        """Build comparison data between initial and final descriptors."""
        return js_data.build_descriptor_comparison_data(
            initial_detailed, final_detailed
        )

    def _build_sankey_data(self, data: dict[str, Any]) -> dict[str, Any]:
        """Build Sankey JSON data for all/models/compare modes."""
        funnel_by_model = data.get("funnel_by_model", {})
        available_models = data.get("available_models", [])

        sankey_by_model = {"all": plots.plot_sankey_json(data.get("funnel", []))}
        for model, model_funnel in funnel_by_model.items():
            sankey_by_model[model] = plots.plot_sankey_json(model_funnel)

        if available_models and funnel_by_model:
            sankey_by_model["__compare__"] = plots.plot_sankey_compare_json(
                funnel_by_model, available_models
            )
        return sankey_by_model

    def _build_synthesis_js_data(self, data: dict[str, Any]) -> dict[str, Any]:
        """Build synthesis JSON data for JavaScript model filtering."""
        return js_data.build_synthesis_js_data(data.get("synthesis_detailed", {}))

    def _build_docking_js_data(
        self, data: dict[str, Any], docking_detailed: dict[str, Any]
    ) -> dict[str, Any]:
        """Build docking JSON data for JavaScript model filtering."""
        return js_data.build_docking_js_data(data, docking_detailed)

    def _build_docking_filters_plots(self, data: dict[str, Any]) -> dict[str, Any]:
        """Build docking filter plots and corresponding JavaScript data."""
        docking_filter_plots = {}
        df_detailed = data.get("docking_filters_detailed", {})

        if df_detailed.get("per_filter"):
            docking_filter_plots["docking_filters_pass_fail"] = (
                plots.plot_docking_filters_pass_fail_bar(df_detailed["per_filter"])
            )
        if df_detailed.get("numeric_metrics"):
            docking_filter_plots["docking_filters_metric_hists"] = (
                plots.plot_docking_filters_metric_histograms(
                    df_detailed["numeric_metrics"],
                    df_detailed.get("thresholds", {}),
                )
            )
        if df_detailed.get("by_model"):
            docking_filter_plots["docking_filters_by_model"] = (
                plots.plot_docking_filters_by_model_bar(df_detailed["by_model"])
            )
            df_by_model_js = {"all": df_detailed}
            for model, model_data in df_detailed["by_model"].items():
                df_by_model_js[model] = model_data
            docking_filter_plots["docking_filters_data"] = df_by_model_js

        return docking_filter_plots

    def _build_synthesis_thresholds(self, data: dict[str, Any]) -> dict[str, Any]:
        """Build synthesis threshold configuration for JavaScript overlays."""
        synth_config = data_collector.load_stage_config(
            self.base_path,
            self.config,
            self.stages,
            self.initial_count,
            self.final_count,
            self.output_dir,
            "config_synthesis",
        )
        if not synth_config:
            return {}

        score_keys = [
            ("sa_scores", "sa_score"),
            ("syba_scores", "syba_score"),
            ("ra_scores", "ra_score"),
        ]
        thresholds = {}
        for score_key, prefix in score_keys:
            entry = {}
            for bound in ["min", "max"]:
                val = synth_config.get(f"{prefix}_{bound}")
                if val is not None:
                    entry[bound] = val
            if entry:
                thresholds[score_key] = entry
        return thresholds

    def _generate_plots(self, data: dict[str, Any]) -> dict[str, Any]:
        """Generate all visualization plots.

        Args:
            data: Collected pipeline data

        Returns:
            Dictionary mapping plot names to HTML strings
        """
        plot_htmls = {}
        for plot_name, plot_builder in _PLOT_REGISTRY.items():
            result = plot_builder(data)
            if result:
                plot_htmls[plot_name] = result

        plot_htmls["sankey_data"] = self._build_sankey_data(data)

        docking_detailed = data.get("docking_detailed", {})
        for tool in ["gnina", "smina"]:
            tool_data = docking_detailed.get(tool, {})
            raw_data = tool_data.get("raw_data", [])
            scores = [
                d.get("affinity") for d in raw_data if d.get("affinity") is not None
            ]
            docking_plot_registry: dict[str, Callable[[], str]] = {
                f"docking_{tool}_affinity_hist": (
                    lambda tool=tool, scores=scores: (
                        plots.plot_docking_affinity_histogram(scores, tool)
                        if scores
                        else ""
                    )
                ),
                f"docking_{tool}_affinity_box": (
                    lambda raw_data=raw_data: (
                        plots.plot_docking_affinity_box(raw_data) if raw_data else ""
                    )
                ),
                f"docking_{tool}_top_molecules": (
                    lambda top_molecules=tool_data.get("top_molecules", []): (
                        plots.plot_docking_top_molecules(top_molecules)
                        if top_molecules
                        else ""
                    )
                ),
            }
            for plot_name, plot_builder in docking_plot_registry.items():
                result = plot_builder()
                if result:
                    plot_htmls[plot_name] = result

        synthesis_js_data = self._build_synthesis_js_data(data)
        if synthesis_js_data:
            plot_htmls["synthesis_data"] = synthesis_js_data

        desc_detailed = data.get("descriptors_detailed", {})
        available_models = data.get("available_models", [])
        if desc_detailed.get("raw_data"):
            plot_htmls["descriptors_data"] = self._build_descriptors_js_data(
                desc_detailed, available_models
            )
        if available_models:
            filters_js_data = self._build_filters_js_data(available_models)
            if filters_js_data.get("filter_data"):
                plot_htmls["filters_data"] = filters_js_data

        plot_htmls.update(self._build_docking_js_data(data, docking_detailed))
        plot_htmls.update(self._build_docking_filters_plots(data))

        initial_desc = data.get("descriptors_detailed", {})
        final_desc = data.get("descriptors_final_detailed", {})
        if initial_desc.get("raw_data") and final_desc.get("raw_data"):
            comparison = self._build_descriptor_comparison_data(
                initial_desc, final_desc
            )
            if comparison:
                plot_htmls["descriptors_comparison_data"] = comparison

        synthesis_thresholds = self._build_synthesis_thresholds(data)
        if synthesis_thresholds:
            plot_htmls["synthesis_thresholds"] = synthesis_thresholds

        return plot_htmls

    def _render_template(self, data: dict[str, Any], plot_htmls: dict[str, Any]) -> str:
        """Render the HTML template with data and plots.

        Args:
            data: Collected pipeline data
            plot_htmls: Dictionary of plot HTML strings

        Returns:
            Rendered HTML string
        """
        try:
            env = Environment(
                loader=PackageLoader("hedgehog.reporting", "templates"),
                autoescape=True,
            )
            template = env.get_template("report.html")
        except (OSError, ImportError, ValueError) as e:
            logger.warning("Could not load template, using inline template: %s", e)
            return self._render_inline_template(data, plot_htmls)

        return template.render(
            data=data,
            plots=plot_htmls,
            generated_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        )

    def _build_model_options(self, models: list[str]) -> str:
        """Build HTML option elements for model dropdown.

        Args:
            models: List of model names

        Returns:
            HTML string with option elements
        """
        return js_data.build_model_options(models)

    def _render_inline_template(
        self, data: dict[str, Any], plot_htmls: dict[str, Any]
    ) -> str:
        """Render report using inline template (fallback).

        Args:
            data: Collected pipeline data
            plot_htmls: Dictionary of plot HTML strings

        Returns:
            Rendered HTML string
        """
        summary = data.get("summary", {})
        metadata = data.get("metadata", {})

        def _json_for_script(obj: Any) -> str:
            # Prevent breaking out of <script> context.
            return json.dumps(obj).replace("</", "<\\/")

        sankey_data_json = _json_for_script(plot_htmls.get("sankey_data") or {})

        # Build stage status HTML
        stage_rows = ""
        for status in summary.get("stage_statuses", []):
            icon = "✓" if status["completed"] else ("✗" if status["enabled"] else "−")
            color = (
                "#2ecc71"
                if status["completed"]
                else ("#e74c3c" if status["enabled"] else "#95a5a6")
            )
            stage_name = html_lib.escape(str(status.get("name", "")))
            stage_status = html_lib.escape(str(status.get("status", "")).upper())
            stage_rows += f"""
            <tr>
                <td>{stage_name}</td>
                <td style="color: {color}; font-weight: bold;">{icon} {stage_status}</td>
            </tr>
            """

        # Build descriptor summary table
        desc_summary = data.get("descriptors", {}).get("summary", {})
        desc_rows = ""
        for desc, stats in desc_summary.items():
            desc_name = html_lib.escape(str(desc))
            desc_rows += f"""
            <tr>
                <td>{desc_name}</td>
                <td>{stats["mean"]:.2f}</td>
                <td>{stats["std"]:.2f}</td>
                <td>{stats["min"]:.2f}</td>
                <td>{stats["max"]:.2f}</td>
            </tr>
            """

        # Build model table
        model_rows = ""
        for model in data.get("models", []):
            retention = (
                100 * model["final"] / model["initial"] if model["initial"] > 0 else 0
            )
            model_rows += f"""
            <tr>
                <td>{model["model_name"]}</td>
                <td>{model["initial"]}</td>
                <td>{model["final"]}</td>
                <td>{retention:.1f}%</td>
            </tr>
            """

        html = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>HEDGEHOG Pipeline Report</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        :root {{
            --primary: #3498db;
            --success: #2ecc71;
            --danger: #e74c3c;
            --warning: #f39c12;
            --dark: #2c3e50;
            --light: #ecf0f1;
        }}
        * {{ box-sizing: border-box; margin: 0; padding: 0; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
            line-height: 1.6;
            color: var(--dark);
            background: #f5f6fa;
        }}
        .container {{ max-width: 1400px; margin: 0 auto; padding: 20px; }}
        header {{
            background: linear-gradient(135deg, var(--primary), #2980b9);
            color: white;
            padding: 40px 20px;
            margin-bottom: 30px;
            border-radius: 10px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}
        header h1 {{ font-size: 2.5em; margin-bottom: 10px; }}
        header .meta {{ opacity: 0.9; font-size: 0.95em; }}
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        .summary-card {{
            background: white;
            padding: 25px;
            border-radius: 10px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
            text-align: center;
        }}
        .summary-card .value {{
            font-size: 2.5em;
            font-weight: bold;
            color: var(--primary);
        }}
        .summary-card .label {{
            color: #7f8c8d;
            font-size: 0.9em;
            text-transform: uppercase;
            letter-spacing: 1px;
        }}
        section {{
            background: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
        }}
        section h2 {{
            color: var(--dark);
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 2px solid var(--light);
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 15px;
        }}
        th, td {{
            padding: 12px 15px;
            text-align: left;
            border-bottom: 1px solid var(--light);
        }}
        th {{
            background: var(--light);
            font-weight: 600;
            text-transform: uppercase;
            font-size: 0.85em;
            letter-spacing: 0.5px;
        }}
        tr:hover {{ background: #f8f9fa; }}
        .plot-container {{ margin: 20px 0; }}
        .two-col {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
        @media (max-width: 768px) {{
            .two-col {{ grid-template-columns: 1fr; }}
            .summary-grid {{ grid-template-columns: repeat(2, 1fr); }}
        }}
        .retention-badge {{
            display: inline-block;
            padding: 5px 15px;
            border-radius: 20px;
            font-weight: bold;
        }}
        .retention-high {{ background: #d5f4e6; color: #27ae60; }}
        .retention-medium {{ background: #fef9e7; color: #f39c12; }}
        .retention-low {{ background: #fadbd8; color: #e74c3c; }}
        .model-filter {{
            display: flex;
            align-items: center;
            gap: 10px;
            margin-bottom: 15px;
        }}
        .model-filter label {{
            font-weight: 600;
            color: var(--dark);
        }}
        .model-filter select {{
            padding: 8px 12px;
            border: 1px solid var(--light);
            border-radius: 5px;
            font-size: 14px;
            background: white;
            cursor: pointer;
            min-width: 200px;
        }}
        .model-filter select:hover {{
            border-color: var(--primary);
        }}
        .view-toggle {{
            display: flex;
            gap: 10px;
            margin-bottom: 15px;
        }}
        .view-toggle button {{
            padding: 8px 16px;
            border: 1px solid var(--light);
            border-radius: 5px;
            background: white;
            cursor: pointer;
            font-size: 14px;
            transition: all 0.2s;
        }}
        .view-toggle button.active {{
            background: var(--primary);
            color: white;
            border-color: var(--primary);
        }}
        .view-toggle button:hover:not(.active) {{
            border-color: var(--primary);
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>HEDGEHOG Pipeline Report</h1>
            <div class="meta">
                <p>Generated: {metadata.get("generated_at", "N/A")}</p>
                <p>Run path: {metadata.get("run_path", "N/A")}</p>
            </div>
        </header>

        <div class="summary-grid">
            <div class="summary-card">
                <div class="value">{summary.get("initial_molecules", 0)}</div>
                <div class="label">Initial Molecules</div>
            </div>
            <div class="summary-card">
                <div class="value">{summary.get("final_molecules", 0)}</div>
                <div class="label">Final Molecules</div>
            </div>
            <div class="summary-card">
                <div class="value">{summary.get("retention_percent", "N/A")}</div>
                <div class="label">Retention Rate</div>
            </div>
            <div class="summary-card">
                <div class="value">{summary.get("stages_completed", 0)}/{summary.get("stages_enabled", 0)}</div>
                <div class="label">Stages Completed</div>
            </div>
        </div>

        <section>
            <h2>Pipeline Flow</h2>
            <div class="model-filter">
                <label for="model-select">Filter by Model:</label>
                <select id="model-select" onchange="updateSankeyDiagram()">
                    <option value="all">All Models</option>
                    {self._build_model_options(data.get("available_models", []))}
                </select>
            </div>
            <div class="view-toggle">
                <button id="btn-sankey" class="active" onclick="showView('sankey')">Sankey Diagram</button>
                <button id="btn-funnel" onclick="showView('funnel')">Funnel Chart</button>
            </div>
            <div id="sankey-container" class="plot-container">{plot_htmls.get("sankey", "")}</div>
            <div id="funnel-container" class="plot-container" style="display: none;">{plot_htmls.get("funnel", "")}</div>
            <script>
                var sankeyDataByModel = {sankey_data_json};

                function showView(view) {{
                    document.getElementById('sankey-container').style.display = view === 'sankey' ? 'block' : 'none';
                    document.getElementById('funnel-container').style.display = view === 'funnel' ? 'block' : 'none';
                    document.getElementById('btn-sankey').className = view === 'sankey' ? 'active' : '';
                    document.getElementById('btn-funnel').className = view === 'funnel' ? 'active' : '';
                }}

                function updateSankeyDiagram() {{
                    var model = document.getElementById('model-select').value;
                    var data = sankeyDataByModel[model];
                    var container = document.getElementById('sankey-container');
                    if (!data || !data.labels || data.labels.length < 2) {{
                        container.textContent = 'No data available for this model';
                        container.style.textAlign = 'center';
                        container.style.color = 'gray';
                        container.style.padding = '40px';
                        return;
                    }}

                    var trace = {{
                        type: 'sankey',
                        arrangement: 'snap',
                        node: {{
                            pad: 20,
                            thickness: 20,
                            line: {{ color: 'black', width: 0.5 }},
                            label: data.labels,
                            color: data.node_colors,
                            x: data.x_positions,
                            y: data.y_positions
                        }},
                        link: {{
                            source: data.sources,
                            target: data.targets,
                            value: data.values,
                            color: data.link_colors
                        }}
                    }};

                    var layout = {{
                        font: {{ size: 12 }},
                        height: 450,
                        margin: {{ l: 20, r: 20, t: 20, b: 20 }}
                    }};

                    Plotly.react('sankey-container', [trace], layout);
                }}
            </script>
        </section>

        <section>
            <h2>Stage Summary</h2>
            <div class="two-col">
                <div>
                    <table>
                        <thead>
                            <tr><th>Stage</th><th>Status</th></tr>
                        </thead>
                        <tbody>{stage_rows}</tbody>
                    </table>
                </div>
                <div class="plot-container">{plot_htmls.get("stage_summary", "")}</div>
            </div>
        </section>

        {"<section><h2>Model Comparison</h2><div class='plot-container'>" + plot_htmls.get("model_comparison", "") + "</div><table><thead><tr><th>Model</th><th>Initial</th><th>Final</th><th>Retention</th></tr></thead><tbody>" + model_rows + "</tbody></table><div class='plot-container'>" + plot_htmls.get("model_losses", "") + "</div></section>" if model_rows else ""}

        {"<section><h2>Descriptor Analysis</h2><div class='plot-container'>" + plot_htmls.get("descriptors", "") + "</div><table><thead><tr><th>Descriptor</th><th>Mean</th><th>Std</th><th>Min</th><th>Max</th></tr></thead><tbody>" + desc_rows + "</tbody></table></section>" if desc_rows else ""}

        {"<section><h2>Filter Analysis</h2><div class='two-col'><div class='plot-container'>" + plot_htmls.get("filter_heatmap", "") + "</div><div class='plot-container'>" + plot_htmls.get("filter_failures", "") + "</div></div></section>" if plot_htmls.get("filter_heatmap") or plot_htmls.get("filter_failures") else ""}

        {"<section><h2>Synthesis Scores</h2><div class='plot-container'>" + plot_htmls.get("synthesis_dist", "") + "</div><div class='plot-container'>" + plot_htmls.get("synthesis_scatter", "") + "</div></section>" if plot_htmls.get("synthesis_dist") else ""}

        {"<section><h2>Docking Results</h2><div class='two-col'><div class='plot-container'>" + plot_htmls.get("docking_gnina", "") + "</div><div class='plot-container'>" + plot_htmls.get("docking_smina", "") + "</div></div></section>" if plot_htmls.get("docking_gnina") or plot_htmls.get("docking_smina") else ""}

        <footer style="text-align: center; padding: 20px; color: #7f8c8d; font-size: 0.9em;">
            <p>HEDGEHOG - Hierarchical Evaluation of Drug GEnerators tHrOugh riGorous filtration</p>
        </footer>
    </div>
</body>
</html>
        """
        return html

    def _save_report(self, html_content: str) -> Path:
        """Save the HTML report to file.

        Args:
            html_content: Rendered HTML content

        Returns:
            Path to the saved report
        """
        report_path = self.output_dir / "report.html"
        with open(report_path, "w", encoding="utf-8") as f:
            f.write(html_content)
        return report_path

    def _save_json_data(self, data: dict[str, Any]) -> Path:
        """Save report data as JSON for programmatic access.

        Args:
            data: Collected pipeline data

        Returns:
            Path to the saved JSON file
        """
        json_path = self.output_dir / "report_data.json"

        # Clean data for JSON serialization
        clean_data = self._make_json_serializable(data)

        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(clean_data, f, indent=2)

        return json_path

    def _make_json_serializable(self, obj: Any) -> Any:
        """Convert object to JSON-serializable format.

        Args:
            obj: Object to convert

        Returns:
            JSON-serializable version of the object
        """
        if isinstance(obj, dict):
            return {k: self._make_json_serializable(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [self._make_json_serializable(item) for item in obj]
        elif isinstance(obj, (int, float, str, bool, type(None))):
            return obj
        elif hasattr(obj, "tolist"):  # numpy arrays
            return obj.tolist()
        else:
            return str(obj)
