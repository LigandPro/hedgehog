"""Tests for shipped Stage 3 hard-gate and diagnostic policies."""

from pathlib import Path

import yaml

from hedgehog.struct_filters.utils import filter_alerts

CONFIG_DIR = Path(__file__).resolve().parents[1] / "src/hedgehog/configs"
ALL_COMMON_ALERT_RULESETS = {
    "Alarm-NMR",
    "AlphaScreen-Hitters",
    "BMS",
    "Chelator",
    "DNABinder",
    "Dundee",
    "Electrophilic",
    "Frequent-Hitter",
    "GST-Hitters",
    "Genotoxic-Carcinogenicity",
    "Glaxo",
    "HIS-Hitters",
    "Inpharmatica",
    "LD50-Oral",
    "LINT",
    "LuciferaseInhibitor",
    "MLSMR",
    "Non-Genotoxic-Carcinogenicity",
    "PAINS",
    "Reactive-Unstable-Toxic",
    "Skin",
    "SureChEMBL",
    "Toxicophore",
}
MAIN_PROFILES = (
    "config_structFilters.yml",
    "config_structFilters_exploration.yml",
    "config_structFilters_strict.yml",
)
FILTER_FLAG_KEYS = {
    "common_alerts",
    "bredt",
    "protecting_groups",
    "ring_infraction",
    "stereo_center",
    "undefined_stereo_center",
    "halogenicity",
    "symmetry",
    "molgraph_stats",
    "molcomplexity",
    "NIBR",
    "lilly",
}
DIAGNOSTIC_FILTERS = {
    "common_alerts",
    "molcomplexity",
    "bredt",
    "lilly",
    "protecting_groups",
    "ring_infraction",
    "stereo_center",
    "halogenicity",
    "symmetry",
}
EXPECTED_HARD_BY_PROFILE = {
    "config_structFilters.yml": {
        "common_alerts",
        "NIBR",
        "lilly",
        "molgraph_stats",
        "protecting_groups",
        "undefined_stereo_center",
    },
    "config_structFilters_exploration.yml": {
        "NIBR",
        "molgraph_stats",
        "protecting_groups",
        "undefined_stereo_center",
    },
    "config_structFilters_strict.yml": {
        "common_alerts",
        "NIBR",
        "lilly",
        "molgraph_stats",
        "protecting_groups",
        "bredt",
        "molcomplexity",
        "stereo_center",
        "undefined_stereo_center",
    },
}
EXPECTED_CA_HARD_RULESETS = {
    "config_structFilters.yml": ["PAINS"],
    "config_structFilters_exploration.yml": [],
    "config_structFilters_strict.yml": [
        "PAINS",
        "LD50-Oral",
        "Toxicophore",
        "Skin",
        "MLSMR",
    ],
}
EXPECTED_LILLY_CUTOFF = {
    "config_structFilters.yml": 160,
    "config_structFilters_exploration.yml": 160,
    "config_structFilters_strict.yml": 100,
}


def load_profile(name: str) -> dict:
    """Load one shipped structural-filter profile."""
    return yaml.safe_load((CONFIG_DIR / name).read_text())


def test_main_profiles_declare_expected_hard_gates():
    for name in MAIN_PROFILES:
        config = load_profile(name)
        flags = {
            key.removeprefix("filter_"): value
            for key, value in config.items()
            if key.startswith("filter_") and key != "filter_data"
        }
        assert set(flags) == FILTER_FLAG_KEYS
        assert {key for key, value in flags.items() if value} == EXPECTED_HARD_BY_PROFILE[
            name
        ]
        assert config["nibr_max_severity"] == 10
        assert config["molgraph_max_severity"] == 5
        assert config["lilly_demerit_cutoff"] == EXPECTED_LILLY_CUTOFF[name]


def test_all_common_alert_rules_are_calculated_in_every_main_profile():
    for name in MAIN_PROFILES:
        config = load_profile(name)
        assert config["include_rulesets"] == "all"
        assert config["exclude_smarts"] == []
        assert "exclude_descriptions" not in config
        alerts = filter_alerts(config)
        assert len(alerts) == 2458
        assert alerts["rule_id"].nunique() == 2458
        assert set(alerts["rule_set_name"].astype(str)) == ALL_COMMON_ALERT_RULESETS
        assert config["filter_common_alerts"] is (
            name != "config_structFilters_exploration.yml"
        )
        assert (
            config["common_alerts_filter_include_rulesets"]
            == EXPECTED_CA_HARD_RULESETS[name]
        )
        assert config["common_alerts_filter_exclude_rulesets"] == []


def test_exclude_smarts_drops_catalog_patterns_in_every_main_profile():
    sample_smarts = "C=[CH]OS(=O)(=O)O[#6]"
    for name in MAIN_PROFILES:
        config = load_profile(name)
        full = filter_alerts(config)
        assert sample_smarts in set(full["smarts"].astype(str))
        alerts = filter_alerts({**config, "exclude_smarts": [sample_smarts]})
        assert len(alerts) == len(full) - 1
        assert sample_smarts not in set(alerts["smarts"].astype(str))


def test_requested_diagnostics_are_enabled_across_profiles():
    for name in MAIN_PROFILES:
        config = load_profile(name)
        enabled = {
            key.removeprefix("calculate_")
            for key, value in config.items()
            if key.startswith("calculate_") and value is True
        }
        assert DIAGNOSTIC_FILTERS.issubset(enabled)
        assert config["write_structural_liability_profile"] is True
