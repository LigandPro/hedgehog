"""Report generator for HEDGEHOG pipeline results."""

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd
from jinja2 import Environment, PackageLoader

from hedgehog.reporting import plots

logger = logging.getLogger(__name__)

# Stage directory names
STAGE_DIRS = {
    "descriptors_initial": "stages/01_descriptors_initial",
    "struct_filters_pre": "stages/02_structural_filters_pre",
    "struct_filters_post": "stages/03_structural_filters_post",
    "synthesis": "stages/04_synthesis",
    "docking": "stages/05_docking",
    "descriptors_final": "stages/06_descriptors_final",
}

# Stage display names
STAGE_DISPLAY_NAMES = {
    "descriptors_initial": "Initial Descriptors",
    "struct_filters_pre": "Pre-Descriptors Filters",
    "struct_filters_post": "Post-Descriptors Filters",
    "synthesis": "Synthesis Analysis",
    "docking": "Molecular Docking",
    "descriptors_final": "Final Descriptors",
}

# Key descriptors to show in report
KEY_DESCRIPTORS = ["MolWt", "LogP", "TPSA", "NumHDonors", "NumHAcceptors", "NumRotatableBonds"]


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
        self.output_dir = self.base_path / "output"
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

        logger.info("Report generated: %s", report_path)
        return report_path

    def _collect_data(self) -> dict[str, Any]:
        """Collect all metrics from stage outputs.

        Returns:
            Dictionary with all collected data
        """
        return {
            "metadata": self._get_metadata(),
            "summary": self._get_summary(),
            "funnel": self._get_funnel_data(),
            "stages": self._get_stage_stats(),
            "models": self._get_model_stats(),
            "descriptors": self._get_descriptor_stats(),
            "filters": self._get_filter_stats(),
            "synthesis": self._get_synthesis_stats(),
            "docking": self._get_docking_stats(),
            "config": self._get_config_summary(),
        }

    def _get_metadata(self) -> dict[str, Any]:
        """Get report metadata."""
        return {
            "generated_at": datetime.now().isoformat(),
            "hedgehog_version": "1.0.0",
            "run_path": str(self.base_path),
        }

    def _get_summary(self) -> dict[str, Any]:
        """Get executive summary statistics."""
        retention_rate = (
            self.final_count / self.initial_count if self.initial_count > 0 else 0
        )

        enabled_stages = [s for s in self.stages if s.enabled]
        completed_stages = [s for s in self.stages if s.completed]

        return {
            "initial_molecules": self.initial_count,
            "final_molecules": self.final_count,
            "retention_rate": retention_rate,
            "retention_percent": f"{retention_rate * 100:.2f}%",
            "stages_enabled": len(enabled_stages),
            "stages_completed": len(completed_stages),
            "stage_statuses": [
                {
                    "name": s.name,
                    "enabled": s.enabled,
                    "completed": s.completed,
                    "status": "completed" if s.completed else ("failed" if s.enabled else "disabled"),
                }
                for s in self.stages
            ],
        }

    def _get_funnel_data(self) -> list[dict[str, Any]]:
        """Get molecule funnel data through pipeline stages."""
        funnel = [{"stage": "Initial", "count": self.initial_count}]

        stage_order = [
            ("struct_filters_pre", "Pre-Filters"),
            ("descriptors_initial", "Descriptors"),
            ("struct_filters_post", "Post-Filters"),
            ("synthesis", "Synthesis"),
            ("docking", "Docking"),
        ]

        for stage_key, display_name in stage_order:
            count = self._get_stage_output_count(stage_key)
            if count is not None:
                funnel.append({"stage": display_name, "count": count})

        # Add final count if different from last stage
        if funnel and funnel[-1]["count"] != self.final_count:
            funnel.append({"stage": "Final", "count": self.final_count})

        return funnel

    def _get_stage_output_count(self, stage_key: str) -> int | None:
        """Get molecule count from stage output file."""
        stage_dir = STAGE_DIRS.get(stage_key)
        if not stage_dir:
            return None

        # Try common output file locations
        paths_to_try = [
            self.base_path / stage_dir / "filtered_molecules.csv",
            self.base_path / stage_dir / "filtered" / "filtered_molecules.csv",
        ]

        for path in paths_to_try:
            if path.exists():
                try:
                    df = pd.read_csv(path)
                    return len(df)
                except Exception:
                    continue
        return None

    def _get_stage_stats(self) -> list[dict[str, Any]]:
        """Get pass/fail statistics for each stage."""
        stats = []

        for stage_key, display_name in STAGE_DISPLAY_NAMES.items():
            stage_dir = self.base_path / STAGE_DIRS.get(stage_key, "")
            if not stage_dir.exists():
                continue

            passed = 0
            failed = 0

            # Try to read filtered and failed molecules
            passed_path = stage_dir / "filtered_molecules.csv"
            if not passed_path.exists():
                passed_path = stage_dir / "filtered" / "filtered_molecules.csv"

            failed_path = stage_dir / "failed_molecules.csv"
            if not failed_path.exists():
                failed_path = stage_dir / "filtered" / "failed_molecules.csv"

            if passed_path.exists():
                try:
                    passed = len(pd.read_csv(passed_path))
                except Exception:
                    pass

            if failed_path.exists():
                try:
                    failed = len(pd.read_csv(failed_path))
                except Exception:
                    pass

            if passed > 0 or failed > 0:
                stats.append({
                    "stage": display_name,
                    "passed": passed,
                    "failed": failed,
                    "total": passed + failed,
                })

        return stats

    def _get_model_stats(self) -> list[dict[str, Any]]:
        """Get per-model statistics."""
        # Try to load final molecules to get model breakdown
        final_path = self.output_dir / "final_molecules.csv"
        if not final_path.exists():
            return []

        try:
            df = pd.read_csv(final_path)
        except Exception:
            return []

        if "model_name" not in df.columns:
            return []

        model_stats = []
        for model in df["model_name"].unique():
            final_count = len(df[df["model_name"] == model])

            # Try to get initial count from input
            initial_count = self._get_initial_model_count(model)

            model_stats.append({
                "model_name": model,
                "initial": initial_count or final_count,
                "final": final_count,
                "losses": self._get_model_losses(model),
            })

        return model_stats

    def _get_initial_model_count(self, model: str) -> int | None:
        """Get initial molecule count for a model."""
        input_path = self.base_path / "input" / "sampled_molecules.csv"
        if not input_path.exists():
            return None

        try:
            df = pd.read_csv(input_path)
            if "model_name" in df.columns:
                return len(df[df["model_name"] == model])
        except Exception:
            pass
        return None

    def _get_model_losses(self, model: str) -> dict[str, int]:
        """Get molecule losses per stage for a model."""
        losses = {}

        stage_pairs = [
            ("descriptors", "descriptors_initial"),
            ("struct_filters", "struct_filters_post"),
            ("synthesis", "synthesis"),
            ("docking", "docking"),
        ]

        for loss_key, stage_key in stage_pairs:
            failed_path = self.base_path / STAGE_DIRS.get(stage_key, "") / "failed_molecules.csv"
            if not failed_path.exists():
                failed_path = self.base_path / STAGE_DIRS.get(stage_key, "") / "filtered" / "failed_molecules.csv"

            if failed_path.exists():
                try:
                    df = pd.read_csv(failed_path)
                    if "model_name" in df.columns:
                        losses[loss_key] = len(df[df["model_name"] == model])
                    else:
                        losses[loss_key] = 0
                except Exception:
                    losses[loss_key] = 0
            else:
                losses[loss_key] = 0

        return losses

    def _get_descriptor_stats(self) -> dict[str, Any]:
        """Get descriptor statistics."""
        desc_dir = self.base_path / STAGE_DIRS["descriptors_initial"]

        # Try multiple possible locations
        paths_to_try = [
            desc_dir / "metrics" / "descriptors_all.csv",
            desc_dir / "filtered" / "descriptors_passed.csv",
        ]

        df = None
        for path in paths_to_try:
            if path.exists():
                try:
                    df = pd.read_csv(path)
                    break
                except Exception:
                    continue

        if df is None:
            return {}

        stats = {"distributions": {}, "summary": {}}

        for desc in KEY_DESCRIPTORS:
            if desc in df.columns:
                values = df[desc].dropna().tolist()
                if values:
                    stats["distributions"][desc] = values
                    stats["summary"][desc] = {
                        "mean": float(df[desc].mean()),
                        "std": float(df[desc].std()),
                        "min": float(df[desc].min()),
                        "max": float(df[desc].max()),
                    }

        return stats

    def _get_filter_stats(self) -> dict[str, Any]:
        """Get filter statistics."""
        filter_stats = {"by_filter": {}, "totals": {}, "by_model": {}}

        # Check both pre and post filter directories
        for stage_key in ["struct_filters_pre", "struct_filters_post"]:
            stage_dir = self.base_path / STAGE_DIRS.get(stage_key, "")
            if not stage_dir.exists():
                continue

            # Look for per-filter subdirectories
            for subdir in stage_dir.iterdir():
                if not subdir.is_dir():
                    continue
                if subdir.name in ["plots", "logs"]:
                    continue

                metrics_path = subdir / "metrics.csv"
                if metrics_path.exists():
                    try:
                        df = pd.read_csv(metrics_path)
                        if "failed_count" in df.columns:
                            total_failed = df["failed_count"].sum()
                            filter_stats["totals"][subdir.name] = int(total_failed)

                            # Get per-model breakdown if available
                            if "model_name" in df.columns:
                                model_counts = df.groupby("model_name")["failed_count"].sum().to_dict()
                                filter_stats["by_filter"][subdir.name] = {
                                    str(k): int(v) for k, v in model_counts.items()
                                }
                    except Exception:
                        continue

        return filter_stats

    def _get_synthesis_stats(self) -> dict[str, Any]:
        """Get synthesis analysis statistics."""
        synth_dir = self.base_path / STAGE_DIRS["synthesis"]
        scores_path = synth_dir / "synthesis_scores.csv"

        if not scores_path.exists():
            return {}

        try:
            df = pd.read_csv(scores_path)
        except Exception:
            return {}

        stats = {"distributions": {}, "scatter_data": {}}

        # Get score distributions
        score_columns = ["sa_score", "syba_score", "ra_score", "sc_score"]
        for col in score_columns:
            if col in df.columns:
                values = df[col].dropna().tolist()
                if values:
                    stats["distributions"][col] = values

        # Get scatter plot data
        if "sa_score" in df.columns and "syba_score" in df.columns:
            stats["scatter_data"] = {
                "sa_scores": df["sa_score"].dropna().tolist(),
                "syba_scores": df["syba_score"].dropna().tolist(),
                "model_names": df["model_name"].tolist() if "model_name" in df.columns else [],
            }

        return stats

    def _get_docking_stats(self) -> dict[str, Any]:
        """Get docking statistics."""
        docking_dir = self.base_path / STAGE_DIRS["docking"]

        stats = {"gnina": {}, "smina": {}}

        # Check GNINA results
        gnina_scores = docking_dir / "gnina" / "scores.csv"
        if gnina_scores.exists():
            try:
                df = pd.read_csv(gnina_scores)
                if "score" in df.columns or "affinity" in df.columns:
                    score_col = "score" if "score" in df.columns else "affinity"
                    stats["gnina"] = {
                        "scores": df[score_col].dropna().tolist(),
                        "count": len(df),
                    }
            except Exception:
                pass

        # Check SMINA results
        smina_scores = docking_dir / "smina" / "scores.csv"
        if smina_scores.exists():
            try:
                df = pd.read_csv(smina_scores)
                if "score" in df.columns or "affinity" in df.columns:
                    score_col = "score" if "score" in df.columns else "affinity"
                    stats["smina"] = {
                        "scores": df[score_col].dropna().tolist(),
                        "count": len(df),
                    }
            except Exception:
                pass

        return stats

    def _get_config_summary(self) -> dict[str, Any]:
        """Get configuration summary."""
        return {
            "folder_to_save": self.config.get("folder_to_save", ""),
            "stages_enabled": [s.name for s in self.stages if s.enabled],
        }

    def _generate_plots(self, data: dict[str, Any]) -> dict[str, str]:
        """Generate all visualization plots.

        Args:
            data: Collected pipeline data

        Returns:
            Dictionary mapping plot names to HTML strings
        """
        plot_htmls = {}

        # Funnel chart
        plot_htmls["funnel"] = plots.plot_funnel(data.get("funnel", []))

        # Stage summary
        plot_htmls["stage_summary"] = plots.plot_stage_summary(data.get("stages", []))

        # Model comparison
        model_stats = data.get("models", [])
        if model_stats:
            plot_htmls["model_comparison"] = plots.plot_model_comparison(model_stats)
            plot_htmls["model_losses"] = plots.plot_model_stacked_losses(model_stats)

        # Descriptor distributions
        desc_data = data.get("descriptors", {})
        if desc_data.get("distributions"):
            plot_htmls["descriptors"] = plots.plot_descriptor_distributions(
                desc_data["distributions"]
            )

        # Filter analysis
        filter_data = data.get("filters", {})
        if filter_data.get("by_filter"):
            plot_htmls["filter_heatmap"] = plots.plot_filter_heatmap(
                filter_data["by_filter"]
            )
        if filter_data.get("totals"):
            plot_htmls["filter_failures"] = plots.plot_top_filter_failures(
                filter_data["totals"]
            )

        # Synthesis scores
        synth_data = data.get("synthesis", {})
        if synth_data.get("distributions"):
            plot_htmls["synthesis_dist"] = plots.plot_synthesis_distributions(
                synth_data["distributions"]
            )
        scatter = synth_data.get("scatter_data", {})
        if scatter.get("sa_scores") and scatter.get("syba_scores"):
            plot_htmls["synthesis_scatter"] = plots.plot_synthesis_scatter(
                scatter["sa_scores"],
                scatter["syba_scores"],
                scatter.get("model_names", []),
            )

        # Docking scores
        docking_data = data.get("docking", {})
        if docking_data.get("gnina", {}).get("scores"):
            plot_htmls["docking_gnina"] = plots.plot_docking_distribution(
                docking_data["gnina"]["scores"], "gnina"
            )
        if docking_data.get("smina", {}).get("scores"):
            plot_htmls["docking_smina"] = plots.plot_docking_distribution(
                docking_data["smina"]["scores"], "smina"
            )

        return plot_htmls

    def _render_template(
        self, data: dict[str, Any], plot_htmls: dict[str, str]
    ) -> str:
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
        except Exception as e:
            logger.warning("Could not load template, using inline template: %s", e)
            return self._render_inline_template(data, plot_htmls)

        return template.render(
            data=data,
            plots=plot_htmls,
            generated_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        )

    def _render_inline_template(
        self, data: dict[str, Any], plot_htmls: dict[str, str]
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

        # Build stage status HTML
        stage_rows = ""
        for status in summary.get("stage_statuses", []):
            icon = "✓" if status["completed"] else ("✗" if status["enabled"] else "−")
            color = "#2ecc71" if status["completed"] else ("#e74c3c" if status["enabled"] else "#95a5a6")
            stage_rows += f"""
            <tr>
                <td>{status['name']}</td>
                <td style="color: {color}; font-weight: bold;">{icon} {status['status'].upper()}</td>
            </tr>
            """

        # Build descriptor summary table
        desc_summary = data.get("descriptors", {}).get("summary", {})
        desc_rows = ""
        for desc, stats in desc_summary.items():
            desc_rows += f"""
            <tr>
                <td>{desc}</td>
                <td>{stats['mean']:.2f}</td>
                <td>{stats['std']:.2f}</td>
                <td>{stats['min']:.2f}</td>
                <td>{stats['max']:.2f}</td>
            </tr>
            """

        # Build model table
        model_rows = ""
        for model in data.get("models", []):
            retention = 100 * model["final"] / model["initial"] if model["initial"] > 0 else 0
            model_rows += f"""
            <tr>
                <td>{model['model_name']}</td>
                <td>{model['initial']}</td>
                <td>{model['final']}</td>
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
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>HEDGEHOG Pipeline Report</h1>
            <div class="meta">
                <p>Generated: {metadata.get('generated_at', 'N/A')}</p>
                <p>Run path: {metadata.get('run_path', 'N/A')}</p>
            </div>
        </header>

        <div class="summary-grid">
            <div class="summary-card">
                <div class="value">{summary.get('initial_molecules', 0)}</div>
                <div class="label">Initial Molecules</div>
            </div>
            <div class="summary-card">
                <div class="value">{summary.get('final_molecules', 0)}</div>
                <div class="label">Final Molecules</div>
            </div>
            <div class="summary-card">
                <div class="value">{summary.get('retention_percent', 'N/A')}</div>
                <div class="label">Retention Rate</div>
            </div>
            <div class="summary-card">
                <div class="value">{summary.get('stages_completed', 0)}/{summary.get('stages_enabled', 0)}</div>
                <div class="label">Stages Completed</div>
            </div>
        </div>

        <section>
            <h2>Pipeline Funnel</h2>
            <div class="plot-container">{plot_htmls.get('funnel', '')}</div>
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
                <div class="plot-container">{plot_htmls.get('stage_summary', '')}</div>
            </div>
        </section>

        {"<section><h2>Model Comparison</h2><div class='plot-container'>" + plot_htmls.get('model_comparison', '') + "</div><table><thead><tr><th>Model</th><th>Initial</th><th>Final</th><th>Retention</th></tr></thead><tbody>" + model_rows + "</tbody></table><div class='plot-container'>" + plot_htmls.get('model_losses', '') + "</div></section>" if model_rows else ""}

        {"<section><h2>Descriptor Analysis</h2><div class='plot-container'>" + plot_htmls.get('descriptors', '') + "</div><table><thead><tr><th>Descriptor</th><th>Mean</th><th>Std</th><th>Min</th><th>Max</th></tr></thead><tbody>" + desc_rows + "</tbody></table></section>" if desc_rows else ""}

        {"<section><h2>Filter Analysis</h2><div class='two-col'><div class='plot-container'>" + plot_htmls.get('filter_heatmap', '') + "</div><div class='plot-container'>" + plot_htmls.get('filter_failures', '') + "</div></div></section>" if plot_htmls.get('filter_heatmap') or plot_htmls.get('filter_failures') else ""}

        {"<section><h2>Synthesis Scores</h2><div class='plot-container'>" + plot_htmls.get('synthesis_dist', '') + "</div><div class='plot-container'>" + plot_htmls.get('synthesis_scatter', '') + "</div></section>" if plot_htmls.get('synthesis_dist') else ""}

        {"<section><h2>Docking Results</h2><div class='two-col'><div class='plot-container'>" + plot_htmls.get('docking_gnina', '') + "</div><div class='plot-container'>" + plot_htmls.get('docking_smina', '') + "</div></div></section>" if plot_htmls.get('docking_gnina') or plot_htmls.get('docking_smina') else ""}

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
