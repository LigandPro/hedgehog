"""Tests for report generator docking filters and final descriptors support."""

import json

import pandas as pd
import pytest

from hedgehog.reporting.report_generator import ReportGenerator


@pytest.fixture
def base_path(tmp_path):
    """Create base directory structure for pipeline output."""
    # Create stage directories
    stages_dir = tmp_path / "stages"
    for d in [
        "01_descriptors_initial",
        "06_docking_filters",
        "07_descriptors_final",
    ]:
        (stages_dir / d).mkdir(parents=True)

    # Create input directory with sampled_molecules
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    pd.DataFrame(
        {
            "smiles": ["CCO", "CCN", "CCC"],
            "model_name": ["ModelA", "ModelA", "ModelB"],
        }
    ).to_csv(input_dir / "sampled_molecules.csv", index=False)

    return tmp_path


@pytest.fixture
def report_gen(base_path):
    """Create ReportGenerator with minimal config."""
    return ReportGenerator(
        base_path=base_path,
        stages=[],
        config={
            "docking_filters": {
                "aggregation": {"mode": "all"},
                "pose_quality": {
                    "max_clashes": 2,
                    "max_strain_energy": 50.0,
                },
                "conformer_deviation": {"max_rmsd_to_conformer": 3.0},
                "search_box": {"max_outside_fraction": 0.0},
            }
        },
        initial_count=3,
        final_count=1,
    )


class TestGetDockingFiltersDetailed:
    """Tests for _get_docking_filters_detailed method."""

    def test_returns_empty_when_no_files(self, report_gen):
        """Should return empty dict when metrics.csv doesn't exist."""
        result = report_gen._get_docking_filters_detailed()
        assert result == {}

    def test_returns_correct_structure(self, report_gen, base_path):
        """Should return correct structure from mock metrics.csv."""
        df_dir = base_path / "stages" / "06_docking_filters"

        # Create metrics.csv with realistic data
        metrics_df = pd.DataFrame(
            {
                "mol_idx": [0, 1, 2, 3, 4],
                "model_name": ["ModelA", "ModelA", "ModelB", "ModelB", "ModelA"],
                "clashes": [0, 1, 3, 0, 2],
                "strain_energy": [5.0, 12.0, 60.0, 3.0, 25.0],
                "min_conformer_rmsd": [1.0, 2.5, 4.0, 0.8, 1.5],
                "frac_atoms_outside_box": [0.0, 0.0, 0.1, 0.0, 0.0],
                "pass_search_box": [True, True, False, True, True],
                "pass_pose_quality": [True, True, False, True, False],
                "pass_conformer_deviation": [True, True, False, True, True],
                "pass": [True, True, False, True, False],
            }
        )
        metrics_df.to_csv(df_dir / "metrics.csv", index=False)

        # Create filtered_molecules.csv
        filtered_df = pd.DataFrame(
            {
                "smiles": ["CCO", "CCN"],
                "model_name": ["ModelA", "ModelB"],
                "mol_idx": [0, 3],
            }
        )
        filtered_df.to_csv(df_dir / "filtered_molecules.csv", index=False)

        result = report_gen._get_docking_filters_detailed()

        assert result["total_poses"] == 5
        assert result["passed_poses"] == 3
        assert result["pass_rate"] == 60.0
        assert result["unique_molecules_passed"] == 2
        assert result["aggregation_mode"] == "all"

        # Per-filter stats
        assert "search_box" in result["per_filter"]
        assert result["per_filter"]["search_box"]["passed"] == 4
        assert result["per_filter"]["search_box"]["total"] == 5
        assert result["per_filter"]["search_box"]["pass_rate"] == 80.0

        assert "pose_quality" in result["per_filter"]
        assert result["per_filter"]["pose_quality"]["passed"] == 3
        assert result["per_filter"]["pose_quality"]["total"] == 5

        # Numeric metrics
        assert "clashes" in result["numeric_metrics"]
        assert len(result["numeric_metrics"]["clashes"]) == 5
        assert "strain_energy" in result["numeric_metrics"]

        # By model
        assert "ModelA" in result["by_model"]
        assert result["by_model"]["ModelA"]["total"] == 3
        assert result["by_model"]["ModelA"]["passed"] == 2
        assert "ModelB" in result["by_model"]
        assert result["by_model"]["ModelB"]["total"] == 2

        # Thresholds from config
        assert result["thresholds"]["clashes"]["max"] == 2
        assert result["thresholds"]["strain_energy"]["max"] == 50.0
        assert result["thresholds"]["min_conformer_rmsd"]["max"] == 3.0

    def test_handles_empty_csv(self, report_gen, base_path):
        """Should return empty dict for empty CSV."""
        df_dir = base_path / "stages" / "06_docking_filters"
        pd.DataFrame().to_csv(df_dir / "metrics.csv", index=False)

        result = report_gen._get_docking_filters_detailed()
        assert result == {}

    def test_no_filtered_molecules_file(self, report_gen, base_path):
        """Should set unique_molecules_passed=0 when file is missing."""
        df_dir = base_path / "stages" / "06_docking_filters"

        metrics_df = pd.DataFrame(
            {
                "mol_idx": [0, 1],
                "model_name": ["ModelA", "ModelA"],
                "pass_search_box": [True, False],
                "pass": [True, False],
            }
        )
        metrics_df.to_csv(df_dir / "metrics.csv", index=False)

        result = report_gen._get_docking_filters_detailed()
        assert result["unique_molecules_passed"] == 0
        assert result["total_poses"] == 2

    def test_includes_interaction_summary_when_artifacts_present(
        self, report_gen, base_path
    ):
        """Should include normalized interaction summary when artifact files exist."""
        df_dir = base_path / "stages" / "06_docking_filters"

        pd.DataFrame(
            {
                "mol_idx": [0, 1, 2],
                "model_name": ["ModelA", "ModelA", "ModelB"],
                "pass_search_box": [True, True, True],
                "pass_interactions": [True, True, False],
                "pass": [True, True, False],
            }
        ).to_csv(df_dir / "metrics.csv", index=False)

        pd.DataFrame(
            {
                "smiles": ["CCO", "CCN"],
                "model_name": ["ModelA", "ModelA"],
                "mol_idx": [0, 1],
            }
        ).to_csv(df_dir / "filtered_molecules.csv", index=False)

        pd.DataFrame(
            {
                "mol_idx": [0, 0, 1],
                "source_mol_idx": ["10", "10", "11"],
                "model_name": ["ModelA", "ModelA", "ModelA"],
                "pass": [True, True, True],
                "pass_interactions": [True, True, True],
                "residue": ["ASP:123", "ASP:123", "SER:45"],
                "interaction_type": ["HBDonor", "Hydrophobic", "HBAcceptor"],
                "label": [
                    "ASP:123|HBDonor",
                    "ASP:123|Hydrophobic",
                    "SER:45|HBAcceptor",
                ],
            }
        ).to_csv(df_dir / "interaction_events.csv", index=False)

        pd.DataFrame(
            {
                "residue": ["ASP:123", "SER:45"],
                "total_events": ["2", "1"],
                "unique_poses": ["1", "1"],
                "unique_passed_poses": ["1", "1"],
                "interaction_types": ["2", "1"],
            }
        ).to_csv(df_dir / "interaction_residue_summary.csv", index=False)

        pd.DataFrame(
            {
                "interaction_type": ["HBDonor", "Hydrophobic", "HBAcceptor"],
                "total_events": [1, 1, 1],
                "unique_poses": [1, 1, 1],
                "unique_residues": [1, 1, 1],
            }
        ).to_csv(df_dir / "interaction_type_summary.csv", index=False)

        pd.DataFrame(
            {
                "residue": ["ASP:123", "ASP:123", "SER:45"],
                "interaction_type": ["HBDonor", "Hydrophobic", "HBAcceptor"],
                "event_count": [1, 1, 1],
                "unique_poses": [1, 1, 1],
            }
        ).to_csv(df_dir / "interaction_matrix.csv", index=False)

        (df_dir / "interaction_report_meta.json").write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "total_events": "3",
                    "poses_with_interactions": "2",
                    "unique_residues": "2",
                    "unique_interaction_types": "3",
                }
            )
        )

        result = report_gen._get_docking_filters_detailed()
        summary = result["interaction_summary"]

        assert summary["total_events"] == 3
        assert summary["poses_with_interactions"] == 2
        assert summary["unique_residues"] == 2
        assert summary["unique_interaction_types"] == 3
        assert summary["top_residues"][0]["residue"] == "ASP:123"
        assert summary["top_residues"][0]["total_events"] == 2
        assert len(summary["type_distribution"]) == 3
        assert len(summary["matrix"]) == 3

    def test_interaction_artifact_parsing_is_robust(self, report_gen, base_path):
        """Should tolerate missing/malformed interaction artifact files."""
        df_dir = base_path / "stages" / "06_docking_filters"

        pd.DataFrame(
            {
                "mol_idx": [0],
                "model_name": ["ModelA"],
                "pass_search_box": [True],
                "pass": [True],
            }
        ).to_csv(df_dir / "metrics.csv", index=False)

        pd.DataFrame({"residue": ["ASP:123"]}).to_csv(
            df_dir / "interaction_residue_summary.csv", index=False
        )
        pd.DataFrame(
            {
                "interaction_type": ["HBDonor"],
                "total_events": ["not-a-number"],
            }
        ).to_csv(df_dir / "interaction_type_summary.csv", index=False)
        (df_dir / "interaction_matrix.csv").write_text(
            'residue,interaction_type,event_count\n"ASP:123,HBDonor,2\n'
        )
        (df_dir / "interaction_report_meta.json").write_text("{invalid-json}")

        result = report_gen._get_docking_filters_detailed()
        summary = result.get("interaction_summary", {})

        assert isinstance(summary, dict)
        assert summary["top_residues"][0]["residue"] == "ASP:123"
        assert summary["top_residues"][0]["total_events"] == 0
        assert summary["type_distribution"][0]["interaction_type"] == "HBDonor"
        assert summary["type_distribution"][0]["total_events"] == 0

    def test_interaction_fallback_unique_passed_poses_uses_pass_flag(
        self, report_gen, base_path
    ):
        """Fallback aggregation should compute unique_passed_poses from pass==True only."""
        df_dir = base_path / "stages" / "06_docking_filters"

        pd.DataFrame(
            {
                "mol_idx": [0, 1, 2],
                "model_name": ["ModelA", "ModelA", "ModelA"],
                "pass_search_box": [True, True, True],
                "pass": [True, False, True],
            }
        ).to_csv(df_dir / "metrics.csv", index=False)

        pd.DataFrame(
            {
                "mol_idx": [0, 1, 2],
                "source_mol_idx": ["10", "11", "12"],
                "model_name": ["ModelA", "ModelA", "ModelA"],
                "pass": [True, False, True],
                "pass_interactions": [True, True, True],
                "residue": ["ASP:123", "ASP:123", "SER:45"],
                "interaction_type": ["HBDonor", "Hydrophobic", "HBAcceptor"],
                "label": ["L1", "L2", "L3"],
            }
        ).to_csv(df_dir / "interaction_events.csv", index=False)

        # Intentionally skip residue summary to trigger fallback aggregation.
        result = report_gen._get_docking_filters_detailed()
        summary = result["interaction_summary"]
        residue_rows = {row["residue"]: row for row in summary["top_residues"]}

        assert residue_rows["ASP:123"]["unique_poses"] == 2
        assert residue_rows["ASP:123"]["unique_passed_poses"] == 1
        assert residue_rows["SER:45"]["unique_poses"] == 1
        assert residue_rows["SER:45"]["unique_passed_poses"] == 1

    def test_generate_plots_includes_interaction_plots(self, report_gen):
        """_generate_plots should include interaction plot keys when data exists."""
        data = {
            "docking_filters_detailed": {
                "interaction_summary": {
                    "top_residues": [
                        {"residue": "ASP:123", "total_events": 3, "unique_poses": 2}
                    ],
                    "type_distribution": [
                        {"interaction_type": "HBDonor", "total_events": 3}
                    ],
                    "matrix": [
                        {
                            "residue": "ASP:123",
                            "interaction_type": "HBDonor",
                            "event_count": 3,
                        }
                    ],
                }
            }
        }

        plot_htmls = report_gen._generate_plots(data)
        assert "docking_filters_interactions_top_residues" in plot_htmls
        assert "docking_filters_interactions_type_distribution" in plot_htmls
        assert "docking_filters_interactions_matrix" in plot_htmls


class TestParametrizedDescriptors:
    """Tests for parametrized descriptor methods."""

    def test_get_descriptor_stats_initial(self, report_gen, base_path):
        """Should read from initial descriptors directory."""
        desc_dir = base_path / "stages" / "01_descriptors_initial" / "metrics"
        desc_dir.mkdir(parents=True, exist_ok=True)

        df = pd.DataFrame(
            {
                "smiles": ["CCO", "CCN"],
                "MolWt": [46.07, 45.08],
                "LogP": [-0.31, -0.64],
            }
        )
        df.to_csv(desc_dir / "descriptors_all.csv", index=False)

        result = report_gen._get_descriptor_stats("descriptors_initial")
        assert "MolWt" in result.get("distributions", {})
        assert "LogP" in result.get("distributions", {})

    def test_get_descriptor_stats_final(self, report_gen, base_path):
        """Should read from final descriptors directory."""
        desc_dir = base_path / "stages" / "07_descriptors_final" / "metrics"
        desc_dir.mkdir(parents=True, exist_ok=True)

        df = pd.DataFrame(
            {
                "smiles": ["CCO"],
                "MolWt": [46.07],
                "LogP": [-0.31],
            }
        )
        df.to_csv(desc_dir / "descriptors_all.csv", index=False)

        result = report_gen._get_descriptor_stats("descriptors_final")
        assert "MolWt" in result.get("distributions", {})

    def test_get_descriptors_detailed_final(self, report_gen, base_path):
        """Should read from final descriptors directory."""
        desc_dir = base_path / "stages" / "07_descriptors_final"

        df = pd.DataFrame(
            {
                "smiles": ["CCO", "CCN"],
                "model_name": ["ModelA", "ModelB"],
                "MolWt": [46.07, 45.08],
                "LogP": [-0.31, -0.64],
                "QED": [0.35, 0.40],
            }
        )
        df.to_csv(desc_dir / "filtered_molecules.csv", index=False)

        result = report_gen._get_descriptors_detailed("descriptors_final")
        assert len(result.get("raw_data", [])) == 2
        assert "ModelA" in result.get("summary_by_model", {})

    def test_get_descriptors_detailed_missing_dir(self, report_gen):
        """Should return empty for nonexistent stage directory."""
        result = report_gen._get_descriptors_detailed("descriptors_final")
        # If no CSV found, returns {}
        assert result == {} or result.get("raw_data") == []


class TestBuildDescriptorComparisonData:
    """Tests for _build_descriptor_comparison_data."""

    def test_comparison_with_data(self, report_gen):
        """Should produce comparison data for common descriptors."""
        initial = {
            "raw_data": [
                {"MolWt": 300, "LogP": 2.0, "model_name": "A"},
                {"MolWt": 350, "LogP": 2.5, "model_name": "A"},
            ]
        }
        final = {
            "raw_data": [
                {"MolWt": 310, "LogP": 1.8, "model_name": "A"},
            ]
        }

        result = report_gen._build_descriptor_comparison_data(initial, final)
        assert "MolWt" in result["data"]
        assert "LogP" in result["data"]
        assert len(result["data"]["MolWt"]["initial_values"]) == 2
        assert len(result["data"]["MolWt"]["final_values"]) == 1

    def test_comparison_empty_data(self, report_gen):
        """Should return empty for missing raw_data."""
        result = report_gen._build_descriptor_comparison_data(
            {"raw_data": []}, {"raw_data": []}
        )
        assert result == {}


class TestNewPlotFunctions:
    """Tests for new docking filters plot functions."""

    def test_pass_fail_bar_returns_html(self):
        """plot_docking_filters_pass_fail_bar should return non-empty HTML."""
        from hedgehog.reporting.plots import plot_docking_filters_pass_fail_bar

        per_filter = {
            "search_box": {"passed": 90, "total": 100, "pass_rate": 90.0},
            "pose_quality": {"passed": 70, "total": 100, "pass_rate": 70.0},
        }
        html = plot_docking_filters_pass_fail_bar(per_filter)
        assert html
        assert "plotly" in html.lower() or "div" in html.lower()

    def test_pass_fail_bar_empty_data(self):
        """Should handle empty data gracefully."""
        from hedgehog.reporting.plots import plot_docking_filters_pass_fail_bar

        html = plot_docking_filters_pass_fail_bar({})
        assert html  # Returns empty plot message

    def test_metric_histograms_returns_html(self):
        """plot_docking_filters_metric_histograms should return non-empty HTML."""
        from hedgehog.reporting.plots import plot_docking_filters_metric_histograms

        metrics = {
            "clashes": [0, 1, 2, 0, 3, 1],
            "strain_energy": [5.0, 12.0, 25.0, 3.0, 45.0],
        }
        thresholds = {"clashes": {"max": 2}, "strain_energy": {"max": 50.0}}
        html = plot_docking_filters_metric_histograms(metrics, thresholds)
        assert html
        assert "plotly" in html.lower() or "div" in html.lower()

    def test_metric_histograms_empty_data(self):
        """Should handle empty data gracefully."""
        from hedgehog.reporting.plots import plot_docking_filters_metric_histograms

        html = plot_docking_filters_metric_histograms({})
        assert html

    def test_by_model_bar_returns_html(self):
        """plot_docking_filters_by_model_bar should return non-empty HTML."""
        from hedgehog.reporting.plots import plot_docking_filters_by_model_bar

        by_model = {
            "ModelA": {"total": 50, "passed": 40, "pass_rate": 80.0},
            "ModelB": {"total": 50, "passed": 30, "pass_rate": 60.0},
        }
        html = plot_docking_filters_by_model_bar(by_model)
        assert html
        assert "plotly" in html.lower() or "div" in html.lower()

    def test_by_model_bar_empty_data(self):
        """Should handle empty data gracefully."""
        from hedgehog.reporting.plots import plot_docking_filters_by_model_bar

        html = plot_docking_filters_by_model_bar({})
        assert html

    def test_interaction_top_residues_returns_html(self):
        """Top-residues interaction plot should return non-empty HTML."""
        from hedgehog.reporting.plots import (
            plot_docking_filters_interaction_top_residues_bar,
        )

        html = plot_docking_filters_interaction_top_residues_bar(
            [
                {"residue": "ASP:123", "total_events": 5, "unique_poses": 3},
                {"residue": "SER:45", "total_events": 2, "unique_poses": 2},
            ]
        )
        assert html
        assert "plotly" in html.lower() or "div" in html.lower()

    def test_interaction_top_residues_empty(self):
        """Top-residues interaction plot should handle empty input."""
        from hedgehog.reporting.plots import (
            plot_docking_filters_interaction_top_residues_bar,
        )

        html = plot_docking_filters_interaction_top_residues_bar([])
        assert html

    def test_interaction_type_distribution_returns_html(self):
        """Interaction-type distribution plot should return non-empty HTML."""
        from hedgehog.reporting.plots import (
            plot_docking_filters_interaction_type_distribution_bar,
        )

        html = plot_docking_filters_interaction_type_distribution_bar(
            [
                {
                    "interaction_type": "HBDonor",
                    "total_events": 4,
                    "unique_residues": 2,
                },
                {
                    "interaction_type": "Hydrophobic",
                    "total_events": 3,
                    "unique_residues": 2,
                },
            ]
        )
        assert html
        assert "plotly" in html.lower() or "div" in html.lower()

    def test_interaction_type_distribution_empty(self):
        """Interaction-type distribution plot should handle empty input."""
        from hedgehog.reporting.plots import (
            plot_docking_filters_interaction_type_distribution_bar,
        )

        html = plot_docking_filters_interaction_type_distribution_bar([])
        assert html

    def test_interaction_matrix_heatmap_returns_html(self):
        """Interaction matrix heatmap should return non-empty HTML."""
        from hedgehog.reporting.plots import (
            plot_docking_filters_interaction_matrix_heatmap,
        )

        html = plot_docking_filters_interaction_matrix_heatmap(
            [
                {
                    "residue": "ASP:123",
                    "interaction_type": "HBDonor",
                    "event_count": 2,
                },
                {
                    "residue": "SER:45",
                    "interaction_type": "Hydrophobic",
                    "event_count": 1,
                },
            ]
        )
        assert html
        assert "plotly" in html.lower() or "div" in html.lower()

    def test_interaction_matrix_heatmap_empty(self):
        """Interaction matrix heatmap should handle empty input."""
        from hedgehog.reporting.plots import (
            plot_docking_filters_interaction_matrix_heatmap,
        )

        html = plot_docking_filters_interaction_matrix_heatmap([])
        assert html
