import json

import pandas as pd
import pytest

from hedgehog.struct_filters.utils import (
    attach_structural_enforcement_pass,
    build_structural_enforcement_pass_mask,
    finalize_structural_liability_profile,
    initialize_structural_liability_profile,
    merge_structural_liability_profile,
)


def _extended(passed, **columns):
    payload = {
        "smiles": ["CCO"],
        "model_name": ["reference"],
        "mol_idx": [7],
        "pass": [passed],
    }
    payload.update({name: [value] for name, value in columns.items()})
    return pd.DataFrame(payload)


def test_liability_profile_separates_hard_gate_from_diagnostics():
    base = pd.DataFrame(
        {"smiles": ["CCO"], "model_name": ["reference"], "mol_idx": [7]}
    )
    profile = initialize_structural_liability_profile(base)
    profile = merge_structural_liability_profile(
        profile, "NIBR", _extended(True, severity=2, reasons="minor flag")
    )
    profile = merge_structural_liability_profile(
        profile,
        "molgraph_stats",
        _extended(True, molgraph_max_severity=0),
    )
    profile = merge_structural_liability_profile(
        profile,
        "lilly",
        _extended(False, demerit_score=180, status="exclude", reasons="rule A"),
    )
    profile = finalize_structural_liability_profile(
        profile,
        {"NIBR", "molgraph_stats"},
        ["NIBR", "molgraph_stats", "lilly"],
    )

    row = profile.iloc[0]
    assert bool(row["stage3_hard_pass"]) is True
    assert row["hard_failed_filters"] == ""
    assert row["stage3_hard_filters"] == "NIBR;molgraph_stats"
    assert bool(row["NIBR__filter_enabled"]) is True
    assert bool(row["molgraph_stats__filter_enabled"]) is True
    assert bool(row["lilly__filter_enabled"]) is False
    assert row["structural_warning_count"] == 1
    assert row["diagnostic_failed_filters"] == "lilly"
    assert row["lilly__demerit_score"] == 180
    assert bool(row["nibr_molgraph_policy_pass"]) is True
    assert bool(row["lilly160_molgraph_policy_pass"]) is False


def test_liability_profile_summarizes_but_keeps_individual_common_alert_hits():
    hits = [
        {"ruleset": "PAINS", "rule_id": 10, "smarts": "[N]"},
        {"ruleset": "PAINS", "rule_id": 10, "smarts": "[N]"},
        {"ruleset": "Dundee", "rule_id": 20, "smarts": "[O]"},
    ]
    base = pd.DataFrame(
        {"smiles": ["CCO"], "model_name": ["reference"], "mol_idx": [7]}
    )
    profile = initialize_structural_liability_profile(base)
    profile = merge_structural_liability_profile(
        profile,
        "common_alerts",
        _extended(
            False,
            pass_PAINS=False,
            pass_Dundee=False,
            pass_any=False,
            alert_hits_json=json.dumps(hits),
        ),
    )
    profile = finalize_structural_liability_profile(profile, set(), ["common_alerts"])

    row = profile.iloc[0]
    assert row["common_alert_ruleset_hit_count"] == 2
    assert row["common_alert_match_count"] == 3
    assert row["common_alert_unique_rule_count"] == 2
    assert row["common_alert_rule_ids"] == "10;20"
    assert row["common_alert_ruleset__PAINS__match_count"] == 2
    assert row["common_alert_ruleset__PAINS__unique_rule_count"] == 1
    assert row["common_alert_ruleset__Dundee__match_count"] == 1
    assert row["common_alert_ruleset__Dundee__unique_rule_count"] == 1
    assert json.loads(row["common_alerts__alert_hits_json"]) == hits


def test_common_alert_enforcement_selects_rulesets():
    hits = [
        [
            {"ruleset": "PAINS", "rule_id": 10, "smarts": "[N]"},
            {"ruleset": "Dundee", "rule_id": 20, "smarts": "[O]"},
        ],
        [{"ruleset": "PAINS", "rule_id": 11, "smarts": "[S]"}],
        [],
    ]
    extended = pd.DataFrame(
        {
            "smiles": ["CCN", "CCS", "CCC"],
            "model_name": ["m", "m", "m"],
            "mol_idx": [1, 2, 3],
            "pass_PAINS": [False, False, True],
            "pass_Dundee": [False, True, True],
            "pass": [False, False, True],
            "alert_hits_json": [json.dumps(value) for value in hits],
        }
    )

    pains_only = build_structural_enforcement_pass_mask(
        {"common_alerts_filter_include_rulesets": ["PAINS"]},
        "common_alerts",
        extended,
    )
    assert pains_only["pass"].tolist() == [False, False, True]

    dundee_only = build_structural_enforcement_pass_mask(
        {"common_alerts_filter_include_rulesets": ["Dundee"]},
        "common_alerts",
        extended,
    )
    assert dundee_only["pass"].tolist() == [False, True, True]

    attached, _mask = attach_structural_enforcement_pass(
        {"common_alerts_filter_include_rulesets": ["Dundee"]},
        "common_alerts",
        extended,
    )
    assert attached["pass"].tolist() == [False, False, True]
    assert attached["enforcement_pass"].tolist() == [False, True, True]


def test_common_alert_enforcement_rejects_unknown_ruleset():
    extended = pd.DataFrame(
        {
            "smiles": ["CCC"],
            "model_name": ["m"],
            "mol_idx": [1],
            "pass_PAINS": [True],
            "pass": [True],
            "alert_hits_json": ["[]"],
        }
    )
    with pytest.raises(ValueError, match="Unknown Common Alerts filter ruleset"):
        build_structural_enforcement_pass_mask(
            {"common_alerts_filter_include_rulesets": ["NOT_A_RULESET"]},
            "common_alerts",
            extended,
        )
