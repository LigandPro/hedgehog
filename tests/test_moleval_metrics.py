"""Tests for MolEval generative metrics integration."""

from unittest.mock import patch

import pandas as pd
import pytest

from hedgehog.reporting import moleval_metrics, plots
from hedgehog.reporting.report_generator import ReportGenerator

# =====================================================================
# moleval_metrics module tests
# =====================================================================


class TestIsAvailable:
    def test_returns_bool(self):
        result = moleval_metrics.is_available()
        assert isinstance(result, bool)

    def test_returns_true_when_moleval_installed(self):
        assert moleval_metrics.is_available() is True


class TestFilterMetrics:
    def test_filters_to_configured_metrics(self):
        raw = {
            "#": 100,
            "# valid": 99,
            "# valid & unique": 95,
            "Validity": 0.99,
            "Uniqueness": 0.95,
            "IntDiv1": 0.85,
            "IntDiv2": 0.72,
            "SEDiv": 0.8,
            "ScaffDiv": 0.6,
            "ScaffUniqueness": 0.45,
            "FG": 0.3,
            "RS": 0.4,
            "Filters": 0.7,
        }
        config = {
            "validity": True,
            "uniqueness": True,
            "internal_diversity": True,
            "se_diversity": True,
            "scaffold_diversity": True,
            "functional_groups": True,
            "ring_systems": True,
            "filters": True,
        }

        result = moleval_metrics._filter_metrics(raw, config)

        # All intrinsic metrics should be present
        assert "Validity" in result
        assert "Uniqueness" in result
        assert "IntDiv1" in result
        assert "IntDiv2" in result
        assert "SEDiv" in result
        assert "ScaffDiv" in result
        assert "ScaffUniqueness" in result
        assert "FG" in result
        assert "RS" in result
        assert "Filters" in result

        # Count-based keys should NOT be present
        assert "#" not in result
        assert "# valid" not in result

    def test_respects_disabled_flags(self):
        raw = {
            "Validity": 0.99,
            "Uniqueness": 0.95,
            "IntDiv1": 0.85,
            "IntDiv2": 0.72,
            "FG": 0.3,
            "Filters": 0.7,
        }
        config = {
            "validity": True,
            "uniqueness": False,
            "internal_diversity": False,
            "functional_groups": True,
            "filters": False,
        }

        result = moleval_metrics._filter_metrics(raw, config)

        assert "Validity" in result
        assert "FG" in result
        assert "Uniqueness" not in result
        assert "IntDiv1" not in result
        assert "IntDiv2" not in result
        assert "Filters" not in result

    def test_rounds_to_four_decimals(self):
        raw = {"Validity": 0.999999999, "IntDiv1": 0.123456789}
        config = {"validity": True, "internal_diversity": True}

        result = moleval_metrics._filter_metrics(raw, config)

        assert result["Validity"] == 1.0
        assert result["IntDiv1"] == 0.1235

    def test_empty_raw_returns_empty(self):
        result = moleval_metrics._filter_metrics({}, {"validity": True})
        assert result == {}


class TestComputeStageMetrics:
    def test_run_false_returns_empty(self):
        result = moleval_metrics.compute_stage_metrics(
            {"Input": ["CCO", "c1ccccc1"]},
            {"run": False},
        )
        assert result == {}

    def test_empty_smiles_returns_empty(self):
        result = moleval_metrics.compute_stage_metrics(
            {},
            {"run": True, "n_jobs": 1, "device": "cpu"},
        )
        assert result == {}

    def test_none_smiles_lists_skipped(self):
        result = moleval_metrics.compute_stage_metrics(
            {"Input": [], "Stage2": []},
            {"run": True, "n_jobs": 1, "device": "cpu"},
        )
        assert result == {}

    def test_with_valid_smiles(self):
        smiles = {
            "Input": [
                "CCO",
                "c1ccccc1",
                "CC(=O)O",
                "CCCC",
                "CCCCC",
                "c1ccc(O)cc1",
                "CC(C)O",
                "CCN",
                "CCCO",
                "CCCN",
            ],
        }
        config = {"run": True, "n_jobs": 1, "device": "cpu", "max_molecules": 2000}

        result = moleval_metrics.compute_stage_metrics(smiles, config)

        assert "by_stage" in result
        assert "stages" in result
        assert "metrics" in result
        assert "Input" in result["by_stage"]
        assert "IntDiv1" in result["by_stage"]["Input"]
        assert 0 <= result["by_stage"]["Input"]["IntDiv1"] <= 1
        assert 0 <= result["by_stage"]["Input"]["Validity"] <= 1

    def test_unavailable_moleval_returns_empty(self):
        with patch.object(moleval_metrics, "is_available", return_value=False):
            result = moleval_metrics.compute_stage_metrics(
                {"Input": ["CCO"]},
                {"run": True, "n_jobs": 1},
            )
        assert result == {}

    def test_compute_stage_metrics_deterministic(self):
        """Verify that seeded RNG produces identical subsamples across calls."""
        # Generate a large pool so subsampling actually happens
        base_smiles = [
            "CCO",
            "c1ccccc1",
            "CC(=O)O",
            "CCCC",
            "CCCCC",
            "c1ccc(O)cc1",
            "CC(C)O",
            "CCN",
            "CCCO",
            "CCCN",
        ]
        # Create many unique SMILES by simple variation
        large_pool = [
            f"{s}" for s in base_smiles
        ] * 50  # 500 total (not unique, but list is large)

        # Capture which subsample is passed to GetMetrics.calculate by mocking it
        captured_inputs = []

        def fake_calculate(gen, **kwargs):
            captured_inputs.append(list(gen))
            return {"IntDiv1": 0.85, "Validity": 1.0}

        with patch("hedgehog.vendor.moleval.metrics.metrics.GetMetrics") as MockGM:
            instance = MockGM.return_value
            instance.calculate = fake_calculate

            config = {
                "run": True,
                "n_jobs": 1,
                "device": "cpu",
                "max_molecules": 5,  # Force subsampling
            }

            moleval_metrics.compute_stage_metrics(
                {"Input": large_pool},
                config,
                seed=42,
            )
            moleval_metrics.compute_stage_metrics(
                {"Input": large_pool},
                config,
                seed=42,
            )

        assert len(captured_inputs) == 2
        assert captured_inputs[0] == captured_inputs[1]


# =====================================================================
# ReportGenerator integration tests
# =====================================================================


@pytest.fixture
def moleval_base_path(tmp_path):
    """Create base directory structure with stage CSVs for moleval testing."""
    # Input
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    pd.DataFrame(
        {
            "smiles": [
                "CCO",
                "c1ccccc1",
                "CC(=O)O",
                "CCCC",
                "CCCCC",
                "c1ccc(O)cc1",
                "CC(C)O",
                "CCN",
                "CCCO",
                "CCCN",
            ],
            "model_name": ["M1"] * 5 + ["M2"] * 5,
        }
    ).to_csv(input_dir / "sampled_molecules.csv", index=False)

    # Descriptors stage
    desc_dir = tmp_path / "stages" / "01_descriptors_initial" / "filtered"
    desc_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "smiles": [
                "CCO",
                "c1ccccc1",
                "CC(=O)O",
                "CCCC",
                "CCCCC",
                "c1ccc(O)cc1",
                "CC(C)O",
                "CCN",
            ],
        }
    ).to_csv(desc_dir / "filtered_molecules.csv", index=False)

    # Struct filters
    sf_dir = tmp_path / "stages" / "03_structural_filters_post"
    sf_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "smiles": ["CCO", "c1ccccc1", "CC(=O)O", "CCCC", "CCCCC", "CCN"],
        }
    ).to_csv(sf_dir / "filtered_molecules.csv", index=False)

    return tmp_path


@pytest.fixture
def moleval_report_gen(moleval_base_path):
    """Create ReportGenerator with moleval config."""
    return ReportGenerator(
        base_path=moleval_base_path,
        stages=[],
        config={"config_moleval": ""},
        initial_count=10,
        final_count=6,
    )


class TestCollectStageSmiles:
    def test_reads_input_smiles(self, moleval_report_gen):
        result = moleval_report_gen._collect_stage_smiles()
        assert "Input" in result
        assert len(result["Input"]) == 10

    def test_reads_descriptor_smiles(self, moleval_report_gen):
        result = moleval_report_gen._collect_stage_smiles()
        assert "Descriptors" in result
        assert len(result["Descriptors"]) == 8

    def test_reads_struct_filters_smiles(self, moleval_report_gen):
        result = moleval_report_gen._collect_stage_smiles()
        assert "StructFilters" in result
        assert len(result["StructFilters"]) == 6

    def test_missing_stages_skipped(self, moleval_report_gen):
        result = moleval_report_gen._collect_stage_smiles()
        # Synthesis and DockingFilters directories don't exist
        assert "Synthesis" not in result
        assert "DockingFilters" not in result

    def test_empty_csv_returns_no_smiles(self, moleval_base_path):
        # Create empty CSV in synthesis stage
        synth_dir = moleval_base_path / "stages" / "04_synthesis"
        synth_dir.mkdir(parents=True)
        pd.DataFrame({"smiles": []}).to_csv(
            synth_dir / "filtered_molecules.csv", index=False
        )

        gen = ReportGenerator(
            base_path=moleval_base_path,
            stages=[],
            config={},
            initial_count=10,
            final_count=6,
        )
        result = gen._collect_stage_smiles()
        assert "Synthesis" not in result


class TestReadStageSmiles:
    def test_reads_smiles_column(self, moleval_report_gen, moleval_base_path):
        path = moleval_base_path / "input" / "sampled_molecules.csv"
        result = moleval_report_gen._read_stage_smiles(path)
        assert len(result) == 10
        assert "CCO" in result

    def test_missing_file_returns_empty(self, moleval_report_gen, moleval_base_path):
        path = moleval_base_path / "nonexistent.csv"
        result = moleval_report_gen._read_stage_smiles(path)
        assert result == []

    def test_no_smiles_column_returns_empty(self, moleval_report_gen, tmp_path):
        csv_path = tmp_path / "no_smiles.csv"
        pd.DataFrame({"name": ["A", "B"], "value": [1, 2]}).to_csv(
            csv_path, index=False
        )
        result = moleval_report_gen._read_stage_smiles(csv_path)
        assert result == []


class TestGetMolevalMetrics:
    def test_no_config_returns_empty(self, moleval_report_gen):
        # config_moleval is empty string
        result = moleval_report_gen._get_moleval_metrics()
        assert result == {}

    def test_run_false_returns_empty(self, moleval_base_path, tmp_path):
        import yaml

        config_file = tmp_path / "moleval_off.yml"
        with open(config_file, "w") as f:
            yaml.dump({"run": False}, f)

        gen = ReportGenerator(
            base_path=moleval_base_path,
            stages=[],
            config={"config_moleval": str(config_file)},
            initial_count=10,
            final_count=6,
        )
        result = gen._get_moleval_metrics()
        assert result == {}


class TestCollectDataIncludesMoleval:
    def test_collect_data_has_moleval_key(self, moleval_report_gen):
        data = moleval_report_gen._collect_data()
        assert "moleval" in data


# =====================================================================
# Plot function tests
# =====================================================================


class TestMolevalPlotFunctions:
    def test_heatmap_empty_data(self):
        result = plots.plot_moleval_heatmap({}, [], [])
        assert "No MolEval metrics" in result

    def test_heatmap_with_data(self):
        by_stage = {
            "Input": {"IntDiv1": 0.85, "Validity": 1.0},
            "Descriptors": {"IntDiv1": 0.82, "Validity": 1.0},
        }
        stages = ["Input", "Descriptors"]
        metrics = ["IntDiv1", "Validity"]

        result = plots.plot_moleval_heatmap(by_stage, stages, metrics)
        assert "<div" in result  # Contains plotly HTML

    def test_lines_empty_data(self):
        result = plots.plot_moleval_stage_lines({}, [])
        assert "No MolEval metrics" in result

    def test_lines_with_data(self):
        by_stage = {
            "Input": {"IntDiv1": 0.85, "ScaffDiv": 0.6, "Filters": 0.9},
            "Descriptors": {"IntDiv1": 0.82, "ScaffDiv": 0.55, "Filters": 0.88},
        }
        stages = ["Input", "Descriptors"]

        result = plots.plot_moleval_stage_lines(by_stage, stages)
        assert "<div" in result

    def test_lines_shows_all_metrics_from_data(self):
        by_stage = {
            "Input": {
                "IntDiv1": 0.85,
                "IntDiv2": 0.72,
                "SEDiv": 0.8,
                "ScaffDiv": 0.6,
                "ScaffUniqueness": 0.45,
                "FG": 0.3,
                "RS": 0.4,
                "Filters": 0.7,
            },
            "Descriptors": {
                "IntDiv1": 0.82,
                "IntDiv2": 0.70,
                "SEDiv": 0.78,
                "ScaffDiv": 0.58,
                "ScaffUniqueness": 0.42,
                "FG": 0.28,
                "RS": 0.38,
                "Filters": 0.68,
            },
        }
        stages = ["Input", "Descriptors"]

        result = plots.plot_moleval_stage_lines(by_stage, stages)
        # All 8 metrics should appear as trace names in the plotly HTML
        for metric in [
            "IntDiv1",
            "IntDiv2",
            "SEDiv",
            "ScaffDiv",
            "ScaffUniqueness",
            "FG",
            "RS",
            "Filters",
        ]:
            assert metric in result, f"{metric} not found in line chart HTML"


# =====================================================================
# RUN_INFO.md generation tests
# =====================================================================


class TestGenerateRunInfo:
    def test_creates_run_info_with_moleval(self, tmp_path):
        gen = ReportGenerator(
            base_path=tmp_path,
            stages=[],
            config={},
            initial_count=100,
            final_count=39,
        )
        data = {
            "moleval": {
                "by_stage": {
                    "Input": {"IntDiv1": 0.85, "ScaffDiv": 0.60},
                    "Descriptors": {"IntDiv1": 0.82, "ScaffDiv": 0.55},
                },
                "stages": ["Input", "Descriptors"],
                "metrics": ["IntDiv1", "ScaffDiv"],
            }
        }

        gen._generate_run_info(data)

        run_info = tmp_path / "RUN_INFO.md"
        assert run_info.exists()
        content = run_info.read_text()
        assert "## MolEval Metrics" in content
        assert "IntDiv1" in content
        assert "0.8500" in content

    def test_appends_to_existing_run_info(self, tmp_path):
        run_info = tmp_path / "RUN_INFO.md"
        run_info.write_text("# HEDGEHOG Run Info\n\n## Summary\n\n- Initial: 100\n")

        gen = ReportGenerator(
            base_path=tmp_path,
            stages=[],
            config={},
            initial_count=100,
            final_count=39,
        )
        data = {
            "moleval": {
                "by_stage": {"Input": {"IntDiv1": 0.85}},
                "stages": ["Input"],
                "metrics": ["IntDiv1"],
            }
        }

        gen._generate_run_info(data)

        content = run_info.read_text()
        assert "# HEDGEHOG Run Info" in content
        assert "## Summary" in content
        assert "## MolEval Metrics" in content
        assert "IntDiv1" in content

    def test_no_moleval_data_skips(self, tmp_path):
        gen = ReportGenerator(
            base_path=tmp_path,
            stages=[],
            config={},
            initial_count=100,
            final_count=39,
        )

        gen._generate_run_info({"moleval": {}})

        assert not (tmp_path / "RUN_INFO.md").exists()

    def test_replaces_existing_moleval_section(self, tmp_path):
        run_info = tmp_path / "RUN_INFO.md"
        run_info.write_text(
            "# HEDGEHOG Run Info\n\n"
            "## MolEval Metrics\n\n"
            "| Metric | Input |\n|--------|--------|\n| Old | 0.5 |\n"
        )

        gen = ReportGenerator(
            base_path=tmp_path,
            stages=[],
            config={},
            initial_count=100,
            final_count=39,
        )
        data = {
            "moleval": {
                "by_stage": {"Input": {"IntDiv1": 0.85}},
                "stages": ["Input"],
                "metrics": ["IntDiv1"],
            }
        }

        gen._generate_run_info(data)

        content = run_info.read_text()
        assert content.count("## MolEval Metrics") == 1
        assert "Old" not in content
        assert "IntDiv1" in content


# =====================================================================
# Rules compliance tests
# =====================================================================


class TestComputeRulesCompliance:
    def test_disabled_returns_empty(self):
        result = moleval_metrics.compute_rules_compliance(
            {"Input": ["CCO", "c1ccccc1"]},
            {"rules_compliance": False},
        )
        assert result == {}

    def test_empty_smiles_returns_empty(self):
        result = moleval_metrics.compute_rules_compliance(
            {},
            {"rules_compliance": True},
        )
        assert result == {}

    def test_with_valid_smiles(self):
        smiles = {
            "Input": [
                "CCO",
                "c1ccccc1",
                "CC(=O)O",
                "CCCC",
                "CCCCC",
                "c1ccc(O)cc1",
                "CC(C)O",
                "CCN",
                "CCCO",
                "CCCN",
            ],
        }
        config = {
            "rules_compliance": True,
            "rules_sets": ["rule_of_five", "rule_of_veber"],
            "max_molecules": 2000,
        }
        result = moleval_metrics.compute_rules_compliance(smiles, config)
        assert "by_stage" in result
        assert "stages" in result
        assert "rules" in result
        assert "Input" in result["by_stage"]
        # All small molecules should pass Lipinski
        assert result["by_stage"]["Input"]["Lipinski"] > 0.5

    def test_multi_stage(self):
        smiles = {
            "Input": ["CCO", "c1ccccc1", "CC(=O)O"],
            "Descriptors": ["CCO", "c1ccccc1"],
        }
        config = {
            "rules_compliance": True,
            "rules_sets": ["rule_of_five"],
        }
        result = moleval_metrics.compute_rules_compliance(smiles, config)
        assert "Input" in result["stages"]
        assert "Descriptors" in result["stages"]

    def test_default_rule_sets(self):
        smiles = {"Input": ["CCO", "c1ccccc1"]}
        config = {"rules_compliance": True}
        result = moleval_metrics.compute_rules_compliance(smiles, config)
        assert "by_stage" in result
        # Default rules include Lipinski, Veber, GenDesign
        rules = result["rules"]
        assert "Lipinski" in rules
        assert "Veber" in rules


class TestComputeGroupMatchRates:
    def test_disabled_returns_empty(self):
        result = moleval_metrics.compute_group_match_rates(
            {"Input": ["CCO"]},
            {"chemical_groups": False},
        )
        assert result == {}

    def test_empty_smiles_returns_empty(self):
        result = moleval_metrics.compute_group_match_rates(
            {},
            {"chemical_groups": True},
        )
        assert result == {}

    def test_with_valid_smiles(self):
        smiles = {
            "Input": [
                "c1ccccc1",
                "c1ccc2ccccc2c1",
                "c1ccncc1",
                "c1ccc(O)cc1",
                "CC(C)O",
                "CCN",
            ],
        }
        config = {
            "chemical_groups": True,
            "chemical_group_names": ["rings_in_drugs"],
            "max_molecules": 2000,
        }
        result = moleval_metrics.compute_group_match_rates(smiles, config)
        assert "by_stage" in result
        assert "Input" in result["by_stage"]
        assert "rings_in_drugs" in result["by_stage"]["Input"]
        # Match rate should be between 0 and 1
        rate = result["by_stage"]["Input"]["rings_in_drugs"]
        assert 0.0 <= rate <= 1.0


# =====================================================================
# Plot function tests for rules compliance
# =====================================================================


class TestRulesCompliancePlot:
    def test_empty_data(self):
        result = plots.plot_rules_compliance_lines({}, [], [])
        assert "No compliance data" in result

    def test_with_data(self):
        by_stage = {
            "Input": {"Lipinski": 0.87, "Veber": 0.92},
            "Descriptors": {"Lipinski": 0.90, "Veber": 0.95},
        }
        stages = ["Input", "Descriptors"]
        rules = ["Lipinski", "Veber"]
        result = plots.plot_rules_compliance_lines(by_stage, stages, rules)
        assert "<div" in result
        assert "Lipinski" in result
        assert "Veber" in result
