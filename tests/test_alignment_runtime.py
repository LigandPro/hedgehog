from pathlib import Path

import yaml

from hedgehog._constants import CFG_DOCKING, CFG_STRUCT_FILTERS
from hedgehog.alignment_runtime import alignment_stage_thresholds
from hedgehog.struct_filters.waves.registry import get_aligned_enforced_filters


def test_runtime_finds_adjacent_audit_without_master_metadata(tmp_path: Path):
    docking_config = tmp_path / "config_docking.yml"
    structural_config = tmp_path / "config_structFilters.yml"
    docking_config.write_text("run: true\n", encoding="utf-8")
    structural_config.write_text("run: true\n", encoding="utf-8")
    audit = {
        "stages": {
            "docking": {
                "thresholds": {
                    "score_thresholds": {
                        "gnina": {
                            "score_property": "minimizedAffinity",
                            "max": -8.5,
                        }
                    }
                }
            },
            "struct_filters": {
                "thresholds": {
                    "enabled_rules": [
                        "common_alerts:PAINS",
                        "bredt",
                    ]
                }
            },
        }
    }
    (tmp_path / "alignment_thresholds.yml").write_text(
        yaml.safe_dump(audit, sort_keys=False),
        encoding="utf-8",
    )
    master = {
        CFG_DOCKING: str(docking_config),
        CFG_STRUCT_FILTERS: str(structural_config),
        "alignment": {
            "enabled": False,
            "target_coverage_percent": 95,
        },
    }

    docking = alignment_stage_thresholds(master, CFG_DOCKING, "docking")
    assert docking == audit["stages"]["docking"]["thresholds"]
    assert get_aligned_enforced_filters(master) == {"common_alerts", "bredt"}


def test_structural_config_enforced_filters_are_fallback_policy():
    master = {CFG_STRUCT_FILTERS: "missing.yml"}
    structural = {"enforced_filters": ["NIBR", "molgraph_stats"]}

    assert get_aligned_enforced_filters(master, structural) == {
        "NIBR",
        "molgraph_stats",
    }


def test_target_audit_cannot_change_generic_structural_hard_gate(tmp_path: Path):
    structural_config = tmp_path / "config_structFilters.yml"
    structural_config.write_text("run: true\n", encoding="utf-8")
    audit = {
        "stages": {
            "struct_filters": {
                "thresholds": {
                    "enabled_rules": [
                        "common_alerts:PAINS",
                        "lilly",
                    ]
                }
            }
        }
    }
    (tmp_path / "alignment_thresholds.yml").write_text(
        yaml.safe_dump(audit, sort_keys=False), encoding="utf-8"
    )
    master = {CFG_STRUCT_FILTERS: str(structural_config)}
    generic = {
        "filter_NIBR": True,
        "filter_molgraph_stats": True,
        "filter_common_alerts": False,
        "filter_lilly": True,
        "filter_stereo_center": False,
        "filter_undefined_stereo_center": True,
    }

    assert get_aligned_enforced_filters(master, generic) == {
        "NIBR",
        "lilly",
        "molgraph_stats",
        "undefined_stereo_center",
    }

def test_partial_policy_flags_preserve_legacy_hard_filters():
    structural = {
        "calculate_bredt": True,
        "calculate_protecting_groups": True,
        "filter_bredt": False,
    }

    assert get_aligned_enforced_filters({}, structural) == {"protecting_groups"}


def test_undefined_stereo_toggle_makes_total_stereo_diagnostic_in_legacy_config():
    structural = {
        "calculate_stereo_center": True,
        "filter_undefined_stereo_center": True,
    }

    assert get_aligned_enforced_filters({}, structural) == {
        "undefined_stereo_center"
    }
