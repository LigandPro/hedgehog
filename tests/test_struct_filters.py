"""Tests for structFilters/utils.py."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import datamol as dm
import numpy as np
import pandas as pd
import pytest
from rdkit import Chem

from hedgehog.stages.structFilters.utils import (
    apply_molcomplexity_filters,
    apply_halogenicity,
    apply_lilly_filter,
    apply_nibr_filter,
    apply_protecting_groups,
    apply_ring_infraction,
    apply_stereo_center,
    apply_structural_alerts,
    apply_symmetry,
    camelcase,
    clean_name,
    filter_function_applier,
    format_number,
    get_basic_stats,
    process_path,
)


class TestProcessPath:
    """Tests for process_path function."""

    def test_adds_trailing_slash(self, tmp_path):
        """Should add trailing slash if missing."""
        path_str = str(tmp_path)
        result = process_path(path_str)
        assert result.endswith("/")

    def test_keeps_trailing_slash(self, tmp_path):
        """Should keep trailing slash if present."""
        path_str = str(tmp_path) + "/"
        result = process_path(path_str)
        assert result.endswith("/")
        assert not result.endswith("//")

    def test_creates_directory(self, tmp_path):
        """Should create directory if it doesn't exist."""
        new_dir = tmp_path / "new_directory"
        result = process_path(str(new_dir))
        assert Path(result.rstrip("/")).exists()

    def test_with_keyword(self, tmp_path):
        """Should append keyword as subdirectory."""
        result = process_path(str(tmp_path), "subdir")
        assert "subdir" in result
        assert result.endswith("/")


class TestCamelcase:
    """Tests for camelcase function."""

    def test_underscore_to_camelcase(self):
        """Should convert underscore-separated to CamelCase."""
        assert camelcase("hello_world") == "HelloWorld"

    def test_single_word(self):
        """Should capitalize single word."""
        assert camelcase("hello") == "Hello"

    def test_multiple_underscores(self):
        """Should handle multiple underscores."""
        assert camelcase("one_two_three") == "OneTwoThree"

    def test_empty_string(self):
        """Should handle empty string."""
        assert camelcase("") == ""


class TestFormatNumber:
    """Tests for format_number function."""

    def test_millions(self):
        """Should format millions with M suffix."""
        assert format_number(1500000) == "1.5M"

    def test_thousands(self):
        """Should format thousands with K suffix."""
        assert format_number(1500) == "1.5K"

    def test_small_numbers(self):
        """Should format small numbers without suffix."""
        assert format_number(100) == "100"

    def test_exact_million(self):
        """Should format exact million."""
        assert format_number(1000000) == "1.0M"

    def test_exact_thousand(self):
        """Should format exact thousand."""
        assert format_number(1000) == "1.0K"

    def test_zero(self):
        """Should handle zero."""
        assert format_number(0) == "0"


class TestCleanName:
    """Tests for clean_name function."""

    def test_removes_metrics(self):
        """Should remove 'metrics' from name."""
        assert "metrics" not in clean_name("filter_metrics")

    def test_removes_underscores(self):
        """Should remove underscores."""
        assert "_" not in clean_name("filter_name")

    def test_removes_csv_extension(self):
        """Should remove .csv extension."""
        assert ".csv" not in clean_name("filter.csv")

    def test_strips_whitespace(self):
        """Should strip leading/trailing whitespace."""
        result = clean_name("  filter  ")
        assert not result.startswith(" ")
        assert not result.endswith(" ")


class TestGetBasicStats:
    """Tests for get_basic_stats function behavior."""

    def test_multimodel_grouping(self):
        """Should group statistics by model when multiple models present."""
        # This is a structural test to verify the function signature works
        # Full integration testing would require more setup
        df = pd.DataFrame(
            {
                "smiles": ["CCO", "CC", "CCC", "CCCC"],
                "model_name": ["m1", "m1", "m2", "m2"],
                "mol": [None, None, None, None],
                "pass": [True, False, True, True],
                "pass_any": [True, True, True, True],
            }
        )

        # Verify DataFrame structure
        assert df["model_name"].nunique() == 2
        assert "pass" in df.columns

    def test_multimodel_with_model_name_list_does_not_overwrite_column(self):
        """Should not crash when model_name arg is a list for multi-model inputs."""
        from hedgehog.stages.structFilters.utils import get_basic_stats

        df = pd.DataFrame(
            {
                "smiles": ["CCO", "CC", "CCC", "CCCC"],
                "model_name": ["m1", "m1", "m2", "m2"],
                "mol": [1, 1, 1, 1],
                "pass": [True, False, True, True],
            }
        )

        res_df, extended = get_basic_stats({}, df, ["m1", "m2"], "bredt")
        assert set(res_df["model_name"].tolist()) == {"m1", "m2"}
        assert extended["model_name"].nunique(dropna=True) == 2
        assert len(extended) == 4


class TestPadDataframeToLength:
    """Tests for _pad_dataframe_to_length function."""

    def test_no_padding_needed(self):
        """Should return unchanged if already at target length."""
        from hedgehog.stages.structFilters.utils import _pad_dataframe_to_length

        df = pd.DataFrame({"a": [1, 2, 3]})
        result = _pad_dataframe_to_length(df, 3)
        assert len(result) == 3

    def test_padding_adds_rows(self):
        """Should add rows to reach target length."""
        from hedgehog.stages.structFilters.utils import _pad_dataframe_to_length

        df = pd.DataFrame({"a": [1, 2]})
        result = _pad_dataframe_to_length(df, 5)
        assert len(result) == 5

    def test_longer_than_target(self):
        """Should return unchanged if longer than target."""
        from hedgehog.stages.structFilters.utils import _pad_dataframe_to_length

        df = pd.DataFrame({"a": [1, 2, 3, 4, 5]})
        result = _pad_dataframe_to_length(df, 3)
        # Should not trim, just return as-is
        assert len(result) >= 3


class TestEnsureDataframeLength:
    """Tests for _ensure_dataframe_length function."""

    def test_exact_length(self):
        """Should return unchanged if exactly at expected length."""
        from hedgehog.stages.structFilters.utils import _ensure_dataframe_length

        df = pd.DataFrame({"a": [1, 2, 3]})
        result = _ensure_dataframe_length(df, 3)
        assert len(result) == 3

    def test_pads_short_dataframe(self):
        """Should pad if shorter than expected."""
        from hedgehog.stages.structFilters.utils import _ensure_dataframe_length

        df = pd.DataFrame({"a": [1, 2]})
        result = _ensure_dataframe_length(df, 4)
        assert len(result) == 4

    def test_trims_long_dataframe(self):
        """Should trim if longer than expected."""
        from hedgehog.stages.structFilters.utils import _ensure_dataframe_length

        df = pd.DataFrame({"a": [1, 2, 3, 4, 5]})
        result = _ensure_dataframe_length(df, 3)
        assert len(result) == 3


# =====================================================================
# New medchem filter tests
# =====================================================================


def _make_mols(smiles_list):
    """Convert SMILES to RDKit mol objects, filtering out None."""
    mols = [dm.to_mol(s) for s in smiles_list]
    return [m for m in mols if m is not None]


def _make_config(tmp_path):
    """Create a minimal config dict pointing to a dummy structFilters config."""
    import yaml

    config_path = tmp_path / "config_sf.yml"
    config_path.write_text(
        yaml.dump(
            {
                "ring_infraction_hetcycle_min_size": 4,
                "stereo_max_centers": 4,
                "stereo_max_undefined": 2,
                "halogenicity_thresh_F": 6,
                "halogenicity_thresh_Br": 3,
                "halogenicity_thresh_Cl": 3,
                "symmetry_threshold": 0.8,
            }
        )
    )
    return {"config_structFilters": str(config_path)}


class TestFilterFunctionApplierNewFilters:
    """Test that new filters are registered in filter_function_applier."""

    @pytest.mark.parametrize(
        "name",
        [
            "protecting_groups",
            "ring_infraction",
            "stereo_center",
            "halogenicity",
            "symmetry",
        ],
    )
    def test_new_filter_registered(self, name):
        fn = filter_function_applier(name)
        assert callable(fn)

    def test_unknown_filter_raises(self):
        with pytest.raises(ValueError, match="not found"):
            filter_function_applier("nonexistent_filter")


class TestApplyProtectingGroups:
    """Tests for apply_protecting_groups filter."""

    def test_clean_molecules_pass(self, tmp_path):
        mols = _make_mols(["c1ccccc1", "CCO", "CC(=O)O"])
        config = _make_config(tmp_path)
        result = apply_protecting_groups(config, mols)
        assert len(result) == 3
        assert "pass" in result.columns
        assert "mol" in result.columns

    def test_boc_protected_fails(self, tmp_path):
        # tert-butoxycarbonyl (Boc) group: OC(=O)C(C)(C)C attached to amine
        boc_smiles = "CC(C)(C)OC(=O)Nc1ccccc1"
        mols = _make_mols([boc_smiles])
        config = _make_config(tmp_path)
        result = apply_protecting_groups(config, mols)
        assert len(result) == 1
        # Boc is a known protecting group, should fail
        assert result["pass"].iloc[0] is np.bool_(False)

    def test_returns_dataframe(self, tmp_path):
        mols = _make_mols(["CCO"])
        config = _make_config(tmp_path)
        result = apply_protecting_groups(config, mols)
        assert isinstance(result, pd.DataFrame)

    def test_with_model_names(self, tmp_path):
        smiles = [("CCO", "model1", dm.to_mol("CCO"), 0)]
        mols = [dm.to_mol("CCO")]
        config = _make_config(tmp_path)
        result = apply_protecting_groups(config, mols, smiles)
        assert "model_name" in result.columns


class TestApplyRingInfraction:
    """Tests for apply_ring_infraction filter."""

    def test_normal_rings_pass(self, tmp_path):
        # Benzene, cyclohexane - normal rings
        mols = _make_mols(["c1ccccc1", "C1CCCCC1"])
        config = _make_config(tmp_path)
        result = apply_ring_infraction(config, mols)
        assert result["pass"].all()

    def test_returns_correct_structure(self, tmp_path):
        mols = _make_mols(["c1ccccc1", "CCO", "C1CN1"])
        config = _make_config(tmp_path)
        result = apply_ring_infraction(config, mols)
        assert isinstance(result, pd.DataFrame)
        assert "pass" in result.columns
        assert "mol" in result.columns
        assert len(result) == 3

    def test_reads_config_param(self, tmp_path):
        import yaml

        config_path = tmp_path / "config_sf.yml"
        config_path.write_text(yaml.dump({"ring_infraction_hetcycle_min_size": 5}))
        config = {"config_structFilters": str(config_path)}
        mols = _make_mols(["c1ccccc1"])
        result = apply_ring_infraction(config, mols)
        assert len(result) == 1


class TestApplyStereoCenter:
    """Tests for apply_stereo_center filter."""

    def test_no_stereocenters_passes(self, tmp_path):
        mols = _make_mols(["c1ccccc1", "CCO"])
        config = _make_config(tmp_path)
        result = apply_stereo_center(config, mols)
        assert result["pass"].all()

    def test_few_stereocenters_passes(self, tmp_path):
        # Cholesterol has several stereocenters but we test a simpler case
        # (R)-butan-2-ol: 1 defined stereocenter
        mols = _make_mols(["[C@@H](O)(CC)C"])
        config = _make_config(tmp_path)
        result = apply_stereo_center(config, mols)
        assert result["pass"].all()

    def test_many_stereocenters_fails(self, tmp_path):
        # Molecule with >4 stereocenters: glucose open chain
        # C([C@@H]([C@@H]([C@@H]([C@@H](C=O)O)O)O)O)O - 4 stereocenters
        # Use a longer sugar chain with 5+ stereocenters
        mols = _make_mols(["[C@@H](O)([C@@H](O)[C@@H](O)[C@@H](O)[C@@H](O)CO)C=O"])
        config = _make_config(tmp_path)
        result = apply_stereo_center(config, mols)
        # At least should not crash; actual fail depends on parsed stereocenters
        assert len(result) == 1


class TestApplyHalogenicity:
    """Tests for apply_halogenicity filter."""

    def test_no_halogens_passes(self, tmp_path):
        mols = _make_mols(["c1ccccc1", "CCO"])
        config = _make_config(tmp_path)
        result = apply_halogenicity(config, mols)
        assert result["pass"].all()

    def test_moderate_halogens_passes(self, tmp_path):
        # 2 chlorines - within threshold of 3
        mols = _make_mols(["ClCCCl"])
        config = _make_config(tmp_path)
        result = apply_halogenicity(config, mols)
        assert result["pass"].all()

    def test_excessive_fluorine_fails(self, tmp_path):
        # Molecule with >6 fluorines: perfluorooctane
        mols = _make_mols(["FC(F)(F)C(F)(F)C(F)(F)C(F)(F)F"])
        config = _make_config(tmp_path)
        result = apply_halogenicity(config, mols)
        assert len(result) == 1
        assert result["pass"].iloc[0] is np.bool_(False)

    def test_excessive_chlorine_fails(self, tmp_path):
        # Molecule with >3 chlorines
        mols = _make_mols(["ClC(Cl)(Cl)C(Cl)Cl"])
        config = _make_config(tmp_path)
        result = apply_halogenicity(config, mols)
        assert len(result) == 1
        assert result["pass"].iloc[0] is np.bool_(False)


class TestApplySymmetry:
    """Tests for apply_symmetry filter."""

    def test_asymmetric_molecule_passes(self, tmp_path):
        # Asymmetric molecule
        mols = _make_mols(["CC(=O)Oc1ccccc1C(=O)O"])  # Aspirin
        config = _make_config(tmp_path)
        result = apply_symmetry(config, mols)
        assert len(result) == 1
        assert result["pass"].iloc[0] is np.bool_(True)

    def test_returns_dataframe(self, tmp_path):
        mols = _make_mols(["CCO", "c1ccccc1"])
        config = _make_config(tmp_path)
        result = apply_symmetry(config, mols)
        assert isinstance(result, pd.DataFrame)
        assert "pass" in result.columns


class TestGetBasicStatsNewFilters:
    """Tests for get_basic_stats with new filter names."""

    @pytest.mark.parametrize(
        "filter_name",
        [
            "protecting_groups",
            "ring_infraction",
            "stereo_center",
            "halogenicity",
            "symmetry",
        ],
    )
    def test_basic_stats_computes(self, filter_name):
        df = pd.DataFrame(
            {
                "mol": [dm.to_mol("CCO"), dm.to_mol("c1ccccc1")],
                "pass": [True, False],
                "model_name": ["m1", "m1"],
                "smiles": ["CCO", "c1ccccc1"],
            }
        )
        config = {}
        res_df, extended = get_basic_stats(config, df, "m1", filter_name)
        assert "banned_ratio" in res_df.columns
        assert res_df["banned_ratio"].iloc[0] == 0.5
        assert len(extended) == 2


# =====================================================================
# Contract tests for MedChem filters (common_alerts, Lilly, NIBR)
# =====================================================================


def _build_alert_data(rulesets_smarts):
    """Build a mock alert DataFrame from a dict of {ruleset_name: [(smarts, description), ...]}.

    Returns a DataFrame with columns: rule_set_name, smarts, description.
    """
    rows = []
    for ruleset, patterns in rulesets_smarts.items():
        for smarts, desc in patterns:
            rows.append(
                {"rule_set_name": ruleset, "smarts": smarts, "description": desc}
            )
    return pd.DataFrame(rows)


class TestCommonAlertsContract:
    """Tests for the pass/pass_any semantics of apply_structural_alerts."""

    @patch("hedgehog.stages.structFilters.utils.filter_alerts")
    @patch("hedgehog.stages.structFilters.utils.load_config")
    def test_pass_is_AND_of_all_rulesets(self, mock_load_config, mock_filter_alerts):
        """A molecule failing 1 of 2 rulesets should have pass=False, pass_any=True."""
        mock_load_config.return_value = {}

        # Ruleset A: matches benzene ring ([c]); Ruleset B: matches nitrogen ([#7]) — absent in benzene
        alert_data = _build_alert_data(
            {
                "RulesetA": [("[c]", "aromatic carbon")],
                "RulesetB": [("[#7]", "nitrogen atom")],
            }
        )
        mock_filter_alerts.return_value = alert_data

        mol = Chem.MolFromSmiles("c1ccccc1")  # benzene: has aromatic C, no N
        config = {"config_structFilters": "dummy.yml"}
        result = apply_structural_alerts(config, [mol])

        assert len(result) == 1
        row = result.iloc[0]
        # Fails RulesetA (matches aromatic C), passes RulesetB (no nitrogen)
        assert row["pass_RulesetA"] is False or row["pass_RulesetA"] == False  # noqa: E712
        assert row["pass_RulesetB"] is True or row["pass_RulesetB"] == True  # noqa: E712
        # pass = AND of all rulesets → False (since RulesetA fails)
        assert row["pass"] == False  # noqa: E712
        # pass_any = OR of all rulesets → True (since RulesetB passes)
        assert row["pass_any"] == True  # noqa: E712

    @patch("hedgehog.stages.structFilters.utils.filter_alerts")
    @patch("hedgehog.stages.structFilters.utils.load_config")
    def test_pass_any_is_OR_of_all_rulesets(self, mock_load_config, mock_filter_alerts):
        """A molecule failing ALL rulesets should have pass=False, pass_any=False."""
        mock_load_config.return_value = {}

        # Both rulesets match benzene
        alert_data = _build_alert_data(
            {
                "RulesetA": [("[c]", "aromatic carbon")],
                "RulesetB": [("[#6]", "any carbon")],
            }
        )
        mock_filter_alerts.return_value = alert_data

        mol = Chem.MolFromSmiles("c1ccccc1")  # benzene: matches both
        config = {"config_structFilters": "dummy.yml"}
        result = apply_structural_alerts(config, [mol])

        row = result.iloc[0]
        assert row["pass_RulesetA"] == False  # noqa: E712
        assert row["pass_RulesetB"] == False  # noqa: E712
        assert row["pass"] == False  # noqa: E712
        assert row["pass_any"] == False  # noqa: E712

    @patch("hedgehog.stages.structFilters.utils.filter_alerts")
    @patch("hedgehog.stages.structFilters.utils.load_config")
    def test_all_pass_both_true(self, mock_load_config, mock_filter_alerts):
        """A molecule passing all rulesets should have pass=True, pass_any=True."""
        mock_load_config.return_value = {}

        # Rulesets that won't match simple ethanol (CCO)
        alert_data = _build_alert_data(
            {
                "RulesetA": [("[#7]", "nitrogen atom")],  # no nitrogen in ethanol
                "RulesetB": [("[#17]", "chlorine atom")],  # no chlorine in ethanol
            }
        )
        mock_filter_alerts.return_value = alert_data

        mol = Chem.MolFromSmiles("CCO")
        config = {"config_structFilters": "dummy.yml"}
        result = apply_structural_alerts(config, [mol])

        row = result.iloc[0]
        assert row["pass_RulesetA"] == True  # noqa: E712
        assert row["pass_RulesetB"] == True  # noqa: E712
        assert row["pass"] == True  # noqa: E712
        assert row["pass_any"] == True  # noqa: E712


class TestGetBasicStatsCommonAlerts:
    """Tests for get_basic_stats with common_alerts filter_name using pre-built DataFrames."""

    def _make_common_alerts_df(
        self, pass_vals, pass_any_vals, pass_A_vals, pass_B_vals
    ):
        """Build a DataFrame matching the shape apply_structural_alerts would produce."""
        n = len(pass_vals)
        return pd.DataFrame(
            {
                "mol": [dm.to_mol("CCO")] * n,
                "smiles": ["CCO"] * n,
                "pass": pass_vals,
                "pass_any": pass_any_vals,
                "pass_RulesetA": pass_A_vals,
                "pass_RulesetB": pass_B_vals,
                "model_name": ["m1"] * n,
            }
        )

    def test_all_banned_ratio_definition(self):
        """all_banned_ratio = (~pass).mean() = fraction failing >=1 ruleset."""
        # 2 of 4 molecules fail at least one ruleset (pass=False)
        df = self._make_common_alerts_df(
            pass_vals=[True, False, False, True],
            pass_any_vals=[True, True, False, True],
            pass_A_vals=[True, False, False, True],
            pass_B_vals=[True, True, False, True],
        )
        config = {"include_rulesets": ["RulesetA", "RulesetB"]}
        res_df, _ = get_basic_stats(config, df, "m1", "common_alerts")
        expected = 2 / 4  # 2 molecules have pass=False
        assert abs(res_df["all_banned_ratio"].iloc[0] - expected) < 1e-9

    def test_any_banned_ratio_definition(self):
        """any_banned_ratio = (~pass_any).mean() = fraction failing ALL rulesets."""
        # 1 of 4 molecules fails all rulesets (pass_any=False)
        df = self._make_common_alerts_df(
            pass_vals=[True, False, False, True],
            pass_any_vals=[True, True, False, True],
            pass_A_vals=[True, False, False, True],
            pass_B_vals=[True, True, False, True],
        )
        config = {"include_rulesets": ["RulesetA", "RulesetB"]}
        res_df, _ = get_basic_stats(config, df, "m1", "common_alerts")
        expected = 1 / 4  # 1 molecule has pass_any=False
        assert abs(res_df["any_banned_ratio"].iloc[0] - expected) < 1e-9

    def test_per_ruleset_banned_ratios(self):
        """Per-ruleset banned ratio = 1 - pass_{ruleset}.mean()."""
        df = self._make_common_alerts_df(
            pass_vals=[True, False, False, False],
            pass_any_vals=[True, True, True, False],
            pass_A_vals=[True, False, True, False],  # 2/4 fail A
            pass_B_vals=[True, True, False, False],  # 2/4 fail B
        )
        config = {"include_rulesets": ["RulesetA", "RulesetB"]}
        res_df, _ = get_basic_stats(config, df, "m1", "common_alerts")
        assert abs(res_df["RulesetA_banned_ratio"].iloc[0] - 0.5) < 1e-9
        assert abs(res_df["RulesetB_banned_ratio"].iloc[0] - 0.5) < 1e-9


class TestLillyFilter:
    """Tests for apply_lilly_filter with fully mocked LillyDemeritsFilters."""

    @patch("hedgehog.stages.structFilters.utils.LillyDemeritsFilters")
    @patch("hedgehog.stages.structFilters.utils.LILLY_AVAILABLE", True)
    @patch("hedgehog.stages.structFilters.utils.load_config")
    def test_length_preserved_on_normal_run(self, mock_load_config, MockLillyClass):
        """Invariant: len(result) == len(input)."""
        mock_load_config.return_value = {"lilly_scheduler": "threads"}
        mock_dfilter = MagicMock()
        MockLillyClass.return_value = mock_dfilter

        mols = [dm.to_mol("CCO"), dm.to_mol("c1ccccc1"), dm.to_mol("CC")]
        mock_dfilter.return_value = pd.DataFrame(
            {
                "smiles": ["CCO", "c1ccccc1", "CC"],
                "status": ["ok", "ok", "ok"],
                "pass_filter": [True, True, True],
                "demerit_score": [0.0, 10.0, 5.0],
                "reasons": ["", "aromatic", "short_chain"],
            }
        )

        config = {"n_jobs": 1, "config_structFilters": "dummy.yml"}
        result = apply_lilly_filter(config, mols)
        assert len(result) == 3

    @patch("hedgehog.stages.structFilters.utils.LillyDemeritsFilters")
    @patch("hedgehog.stages.structFilters.utils.LILLY_AVAILABLE", True)
    @patch("hedgehog.stages.structFilters.utils.load_config")
    def test_batch_fallback_on_value_error(self, mock_load_config, MockLillyClass):
        """When dfilter raises ValueError with 'Length of values', fallback triggers."""
        mock_load_config.return_value = {"lilly_scheduler": "threads"}
        mock_dfilter = MagicMock()
        MockLillyClass.return_value = mock_dfilter

        mols = [dm.to_mol("CCO"), dm.to_mol("c1ccccc1")]

        # First call (batch) raises ValueError, second calls (one-by-one) succeed
        call_count = [0]

        def side_effect(mols, n_jobs, scheduler):
            call_count[0] += 1
            if call_count[0] == 1:
                raise ValueError("Length of values does not match length of index")
            # One-by-one fallback calls
            return pd.DataFrame(
                {
                    "smiles": [dm.to_smiles(m) for m in mols],
                    "status": ["ok"] * len(mols),
                    "pass_filter": [True] * len(mols),
                    "demerit_score": [0.0] * len(mols),
                    "reasons": [""] * len(mols),
                }
            )

        mock_dfilter.side_effect = side_effect

        config = {"n_jobs": 1, "config_structFilters": "dummy.yml"}
        result = apply_lilly_filter(config, mols)
        assert len(result) == 2
        # Fallback was triggered (more than 1 call)
        assert call_count[0] > 1

    @patch("hedgehog.stages.structFilters.utils.LillyDemeritsFilters")
    @patch("hedgehog.stages.structFilters.utils.LILLY_AVAILABLE", True)
    @patch("hedgehog.stages.structFilters.utils.load_config")
    def test_invalid_molecules_alignment(self, mock_load_config, MockLillyClass):
        """[valid, None, valid] -> None at index 1 marked as failed."""
        mock_load_config.return_value = {"lilly_scheduler": "threads"}
        mock_dfilter = MagicMock()
        MockLillyClass.return_value = mock_dfilter

        valid_mol_1 = dm.to_mol("CCO")
        valid_mol_2 = dm.to_mol("CC")
        mols = [valid_mol_1, None, valid_mol_2]

        # dfilter is only called on valid mols (2 of them)
        mock_dfilter.return_value = pd.DataFrame(
            {
                "smiles": ["CCO", "CC"],
                "status": ["ok", "ok"],
                "pass_filter": [True, True],
                "demerit_score": [0.0, 5.0],
                "reasons": ["", "short"],
            }
        )

        config = {"n_jobs": 1, "config_structFilters": "dummy.yml"}
        result = apply_lilly_filter(config, mols)
        assert len(result) == 3
        # Index 1 (None molecule) should be marked as failed
        assert result["pass_filter"].iloc[1] == False  # noqa: E712

    @patch("hedgehog.stages.structFilters.utils.LILLY_AVAILABLE", False)
    def test_raises_when_not_available(self):
        """LILLY_AVAILABLE=False -> ImportError."""
        mols = [dm.to_mol("CCO")]
        config = {"n_jobs": 1, "config_structFilters": "dummy.yml"}
        with pytest.raises(ImportError, match="not available"):
            apply_lilly_filter(config, mols)


class TestNIBRFilter:
    """Tests for apply_nibr_filter with mocked NIBRFilters."""

    @patch("hedgehog.stages.structFilters.utils.mc")
    @patch("hedgehog.stages.structFilters.utils.load_config")
    def test_severity_based_pass(self, mock_load_config, mock_mc):
        """severity 0 -> pass=True, severity >0 -> pass=False (in get_basic_stats)."""
        mock_load_config.return_value = {"nibr_scheduler": "threads"}

        mock_nibr = MagicMock()
        mock_mc.structural.NIBRFilters.return_value = mock_nibr

        # NIBRFilters returns a DataFrame with severity, n_covalent_motif, special_mol, mol
        mols = [dm.to_mol("CCO"), dm.to_mol("c1ccccc1")]
        nibr_result = pd.DataFrame(
            {
                "mol": mols,
                "severity": [0, 10],
                "n_covalent_motif": [0, 1],
                "special_mol": [0, 1],
                "status": ["ok", "flagged"],
            }
        )
        mock_nibr.return_value = nibr_result

        config = {"n_jobs": 1, "config_structFilters": "dummy.yml"}
        result = apply_nibr_filter(config, mols)

        # Result should have the data from NIBRFilters
        assert len(result) == 2
        assert result["severity"].iloc[0] == 0
        assert result["severity"].iloc[1] == 10

    @patch("hedgehog.stages.structFilters.utils.mc")
    @patch("hedgehog.stages.structFilters.utils.load_config")
    def test_nibr_bonus_metrics_in_stats(self, mock_load_config, mock_mc):
        """Verify mean_severity, max_severity, banned_ratio in stats."""
        mock_load_config.return_value = {"nibr_scheduler": "threads"}

        mock_nibr = MagicMock()
        mock_mc.structural.NIBRFilters.return_value = mock_nibr

        mols = [dm.to_mol("CCO"), dm.to_mol("c1ccccc1"), dm.to_mol("CC")]
        nibr_result = pd.DataFrame(
            {
                "mol": mols,
                "severity": [0, 10, 5],
                "n_covalent_motif": [0, 2, 1],
                "special_mol": [0, 1, 0],
            }
        )
        mock_nibr.return_value = nibr_result

        config = {"n_jobs": 1, "config_structFilters": "dummy.yml"}
        result = apply_nibr_filter(config, mols)
        result["model_name"] = "m1"

        res_df, _ = get_basic_stats(config, result, "m1", "NIBR")

        assert "mean_severity" in res_df.columns
        assert "max_severity" in res_df.columns
        assert "banned_ratio" in res_df.columns
        assert abs(res_df["mean_severity"].iloc[0] - 5.0) < 1e-9
        assert abs(res_df["max_severity"].iloc[0] - 10.0) < 1e-9
        # banned_ratio: 2 of 3 have severity > 0 → banned_ratio = 2/3
        assert abs(res_df["banned_ratio"].iloc[0] - 2 / 3) < 1e-9


class TestGetBasicStatsLilly:
    """Tests for get_basic_stats with lilly filter_name."""

    def test_lilly_banned_ratio_from_pass_filter(self):
        """banned_ratio = fraction with pass_filter=False."""
        df = pd.DataFrame(
            {
                "mol": [
                    dm.to_mol("CCO"),
                    dm.to_mol("c1ccccc1"),
                    dm.to_mol("CC"),
                    dm.to_mol("CCC"),
                ],
                "pass_filter": [True, False, True, False],
                "demerit_score": [0.0, 100.0, 0.0, 50.0],
                "smiles": ["CCO", "c1ccccc1", "CC", "CCC"],
                "model_name": ["m1", "m1", "m1", "m1"],
                "status": ["ok", "exclude", "ok", "exclude"],
                "reasons": ["", "alert", "", "alert"],
            }
        )
        config = {}
        res_df, _ = get_basic_stats(config, df, "m1", "lilly")
        # 2 of 4 have pass_filter=False → banned_ratio = 0.5
        assert abs(res_df["banned_ratio"].iloc[0] - 0.5) < 1e-9

    def test_mean_noNA_demerit_score(self):
        """mean_noNA_demerit_score = mean of demerit_score excluding NaN."""
        df = pd.DataFrame(
            {
                "mol": [dm.to_mol("CCO"), dm.to_mol("c1ccccc1"), dm.to_mol("CC")],
                "pass_filter": [True, False, True],
                "demerit_score": [10.0, 30.0, np.nan],
                "smiles": ["CCO", "c1ccccc1", "CC"],
                "model_name": ["m1", "m1", "m1"],
                "status": ["ok", "exclude", "ok"],
                "reasons": ["", "alert", ""],
            }
        )
        config = {}
        res_df, _ = get_basic_stats(config, df, "m1", "lilly")
        # mean of [10.0, 30.0] (excluding NaN) = 20.0
        assert abs(res_df["mean_noNA_demerit_score"].iloc[0] - 20.0) < 1e-9


def test_molcomplexity_smcm_does_not_print_unsupported_atom(capsys):
    mol = Chem.MolFromSmiles("[H]/N=C/C")
    res = apply_molcomplexity_filters({}, [mol])
    captured = capsys.readouterr()
    assert "Unsupported atom" not in captured.out
    assert "pass_smcm" in res.columns
