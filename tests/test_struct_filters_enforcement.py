from hedgehog._constants import CFG_STRUCT_FILTERS
from hedgehog.struct_filters.waves.registry import get_aligned_enforced_filters


def test_structural_config_enforced_filters_are_fallback_policy():
    master = {CFG_STRUCT_FILTERS: "missing.yml"}
    structural = {"enforced_filters": ["NIBR", "molgraph_stats"]}

    assert get_aligned_enforced_filters(master, structural) == {
        "NIBR",
        "molgraph_stats",
    }


def test_explicit_filter_flags_ignore_missing_alignment_audit():
    master = {CFG_STRUCT_FILTERS: "missing.yml"}
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
    assert get_aligned_enforced_filters(master) is None


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

    assert get_aligned_enforced_filters({}, structural) == {"undefined_stereo_center"}
