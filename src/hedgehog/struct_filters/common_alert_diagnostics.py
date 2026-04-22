"""Diagnostic artifacts for Common Alerts structural filtering."""

from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem

HITS_JSON_COLUMN = "alert_hits_json"


HITS_LONG_COLUMNS = [
    "mol_idx",
    "model_name",
    "smiles",
    "match_smiles",
    "ruleset",
    "rule_id",
    "description",
    "smarts",
    "priority",
    "source",
    "catalog_description",
    "match_id",
    "matched_atom_ids",
    "matched_bond_ids",
    "n_matches_for_rule",
    "total_alerts_for_molecule",
    "total_rulesets_failed",
    "is_unique_ruleset_failure",
    "is_unique_description_failure",
    "would_pass_if_disable_ruleset",
    "would_pass_if_disable_description",
]


def bond_ids_between_atoms(mol: Chem.Mol, atom_ids: list[int]) -> list[int]:
    """Return bond ids whose endpoints are both inside the matched atom set."""
    atom_set = set(atom_ids)
    bond_ids: list[int] = []
    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() in atom_set and bond.GetEndAtomIdx() in atom_set:
            bond_ids.append(int(bond.GetIdx()))
    return bond_ids


def hit_records_to_json(records: list[dict[str, Any]]) -> str:
    """Serialize per-molecule alert hit records for extended.csv."""
    return json.dumps(records, separators=(",", ":"))


def parse_hit_records(value: Any) -> list[dict[str, Any]]:
    """Parse per-molecule alert hit records from extended.csv."""
    if value is None or pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    try:
        parsed = json.loads(text)
    except json.JSONDecodeError:
        return []
    return parsed if isinstance(parsed, list) else []


def _json_list(value: Any) -> str:
    if isinstance(value, str):
        return value
    if value is None:
        value = []
    return json.dumps(list(value), separators=(",", ":"))


def make_compiled_alert_rule(
    row: Any,
    *,
    ruleset: str,
    patt: Chem.Mol,
    smarts: str,
    description: str,
) -> dict[str, Any]:
    """Normalize one alert rule payload for fast worker-side matching."""
    return {
        "ruleset": ruleset,
        "patt": patt,
        "rule_id": row.get("rule_id"),
        "description": description,
        "smarts": smarts,
        "priority": row.get("priority"),
        "source": row.get("source"),
        "catalog_description": row.get("catalog_description"),
    }


def make_hit_record(
    *,
    mol: Chem.Mol,
    rule: dict[str, Any],
    atom_ids: tuple[int, ...] | list[int],
    match_id: int,
    n_matches_for_rule: int,
) -> dict[str, Any]:
    """Build one atom-level hit record for extended/hits_long diagnostics."""
    atom_id_list = [int(atom_id) for atom_id in atom_ids]
    return {
        "ruleset": rule["ruleset"],
        "rule_id": rule.get("rule_id"),
        "description": rule["description"],
        "smarts": rule.get("smarts", ""),
        "priority": rule.get("priority"),
        "source": rule.get("source"),
        "catalog_description": rule.get("catalog_description"),
        "match_smiles": Chem.MolToSmiles(mol, canonical=False),
        "match_id": match_id,
        "matched_atom_ids": atom_id_list,
        "matched_bond_ids": bond_ids_between_atoms(mol, atom_id_list),
        "n_matches_for_rule": n_matches_for_rule,
    }


def _pass_columns(df: pd.DataFrame) -> list[str]:
    return [
        c for c in df.columns if c.startswith("pass_") and c not in {"pass", "pass_any"}
    ]


def _description_key(hit: dict[str, Any]) -> str:
    return f"{hit.get('ruleset', '')}\t{hit.get('description', '')}"


def build_hits_long(extended_df: pd.DataFrame) -> pd.DataFrame:
    """Build one-row-per-alert-match diagnostics from common_alerts extended data."""
    if HITS_JSON_COLUMN not in extended_df.columns:
        return pd.DataFrame(columns=HITS_LONG_COLUMNS)

    pass_cols = _pass_columns(extended_df)
    rows: list[dict[str, Any]] = []

    for _, mol_row in extended_df.iterrows():
        hits = parse_hit_records(mol_row.get(HITS_JSON_COLUMN))
        if not hits:
            continue

        failed_rulesets = {
            col.removeprefix("pass_") for col in pass_cols if not bool(mol_row.get(col))
        }
        description_keys = {_description_key(hit) for hit in hits}
        total_rulesets_failed = len(failed_rulesets)
        total_alerts_for_molecule = len(hits)

        for hit in hits:
            ruleset = str(hit.get("ruleset", ""))
            description_key = _description_key(hit)
            unique_ruleset = total_rulesets_failed == 1 and ruleset in failed_rulesets
            unique_description = (
                len(description_keys) == 1 and description_key in description_keys
            )
            rows.append(
                {
                    "mol_idx": mol_row.get("mol_idx"),
                    "model_name": mol_row.get("model_name"),
                    "smiles": mol_row.get("smiles"),
                    "match_smiles": hit.get("match_smiles", mol_row.get("smiles")),
                    "ruleset": ruleset,
                    "rule_id": hit.get("rule_id"),
                    "description": hit.get("description", ""),
                    "smarts": hit.get("smarts", ""),
                    "priority": hit.get("priority"),
                    "source": hit.get("source"),
                    "catalog_description": hit.get("catalog_description"),
                    "match_id": hit.get("match_id", 0),
                    "matched_atom_ids": _json_list(hit.get("matched_atom_ids")),
                    "matched_bond_ids": _json_list(hit.get("matched_bond_ids")),
                    "n_matches_for_rule": hit.get("n_matches_for_rule", 1),
                    "total_alerts_for_molecule": total_alerts_for_molecule,
                    "total_rulesets_failed": total_rulesets_failed,
                    "is_unique_ruleset_failure": unique_ruleset,
                    "is_unique_description_failure": unique_description,
                    "would_pass_if_disable_ruleset": unique_ruleset,
                    "would_pass_if_disable_description": unique_description,
                }
            )

    return pd.DataFrame(rows, columns=HITS_LONG_COLUMNS)


def _nunique(series: pd.Series) -> int:
    return int(series.nunique(dropna=False))


def _molecule_identity_values(df: pd.DataFrame) -> pd.Series:
    """Build stable molecule identities for analytics across multiple models."""
    values = []
    for idx, row in df.iterrows():
        model_name = row.get("model_name", "")
        if pd.isna(model_name):
            model_name = ""
        mol_idx = row.get("mol_idx")
        if pd.notna(mol_idx):
            values.append((model_name, mol_idx))
            continue
        values.append((model_name, row.get("smiles", idx)))
    return pd.Series(values, index=df.index)


def _with_molecule_identity(hits_long: pd.DataFrame) -> pd.DataFrame:
    hits = hits_long.copy()
    hits["_molecule_identity"] = _molecule_identity_values(hits)
    return hits


def _unique_rescue_counts(
    hits_long: pd.DataFrame, group_cols: str | list[str], flag_col: str
) -> pd.Series:
    """Count rescued molecules once per group, even with multiple alert matches."""
    rescue_hits = hits_long[hits_long[flag_col].astype(bool)]
    if rescue_hits.empty:
        return pd.Series(dtype=int)
    rescue_hits = _with_molecule_identity(rescue_hits)
    return rescue_hits.groupby(group_cols, dropna=False)["_molecule_identity"].nunique()


def build_ruleset_summary(hits_long: pd.DataFrame) -> pd.DataFrame:
    """Summarize total and unique rescue impact by ruleset."""
    if hits_long.empty:
        return pd.DataFrame(
            columns=["ruleset", "total_hits", "unique_rescue", "overlap_hits"]
        )
    hits_with_id = _with_molecule_identity(hits_long)
    grouped = hits_with_id.groupby("ruleset", dropna=False)
    summary = grouped.agg(total_hits=("_molecule_identity", _nunique))
    unique_rescue = _unique_rescue_counts(
        hits_long, "ruleset", "would_pass_if_disable_ruleset"
    )
    summary["unique_rescue"] = unique_rescue
    summary["unique_rescue"] = summary["unique_rescue"].fillna(0).astype(int)
    summary = summary.reset_index()
    summary["overlap_hits"] = summary["total_hits"] - summary["unique_rescue"]
    return summary.sort_values(
        ["unique_rescue", "total_hits", "ruleset"], ascending=[False, False, True]
    )


def build_description_summary(hits_long: pd.DataFrame) -> pd.DataFrame:
    """Summarize total and unique rescue impact by ruleset/description."""
    if hits_long.empty:
        return pd.DataFrame(
            columns=[
                "ruleset",
                "description",
                "total_hits",
                "unique_rescue",
                "overlap_hits",
            ]
        )
    hits_with_id = _with_molecule_identity(hits_long)
    grouped = hits_with_id.groupby(["ruleset", "description"], dropna=False)
    summary = grouped.agg(total_hits=("_molecule_identity", _nunique))
    unique_rescue = _unique_rescue_counts(
        hits_long,
        ["ruleset", "description"],
        "would_pass_if_disable_description",
    )
    summary["unique_rescue"] = unique_rescue
    summary["unique_rescue"] = summary["unique_rescue"].fillna(0).astype(int)
    summary = summary.reset_index()
    summary["overlap_hits"] = summary["total_hits"] - summary["unique_rescue"]
    return summary.sort_values(
        ["unique_rescue", "total_hits", "ruleset", "description"],
        ascending=[False, False, True, True],
    )


def _cooccurrence(hits_long: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    if hits_long.empty:
        return pd.DataFrame(
            columns=[
                "left",
                "right",
                "left_hits",
                "right_hits",
                "intersection",
                "jaccard",
                "containment",
            ]
        )

    hits_with_id = _with_molecule_identity(hits_long)
    mol_sets = {
        key: set(group["_molecule_identity"].tolist())
        for key, group in hits_with_id.groupby(keys, dropna=False)
    }
    normalized_sets = {
        "\t".join(str(part) for part in key)
        if isinstance(key, tuple)
        else str(key): value
        for key, value in mol_sets.items()
    }

    rows: list[dict[str, Any]] = []
    for left, right in combinations(sorted(normalized_sets), 2):
        left_set = normalized_sets[left]
        right_set = normalized_sets[right]
        intersection = len(left_set & right_set)
        union = len(left_set | right_set)
        min_size = min(len(left_set), len(right_set))
        rows.append(
            {
                "left": left,
                "right": right,
                "left_hits": len(left_set),
                "right_hits": len(right_set),
                "intersection": intersection,
                "jaccard": intersection / union if union else 0.0,
                "containment": intersection / min_size if min_size else 0.0,
            }
        )

    return pd.DataFrame(rows).sort_values(
        ["intersection", "jaccard", "left", "right"],
        ascending=[False, False, True, True],
    )


def write_common_alert_diagnostics(output_dir: Path, extended_df: pd.DataFrame) -> None:
    """Write Common Alerts long-format hit data and summary CSV artifacts."""
    output_dir.mkdir(parents=True, exist_ok=True)
    hits_long = build_hits_long(extended_df)
    hits_long.to_csv(output_dir / "hits_long.csv", index=False)
    build_ruleset_summary(hits_long).to_csv(
        output_dir / "ruleset_summary.csv", index=False
    )
    build_description_summary(hits_long).to_csv(
        output_dir / "description_summary.csv", index=False
    )
    _cooccurrence(hits_long, ["ruleset"]).to_csv(
        output_dir / "cooccurrence_ruleset.csv", index=False
    )
    _cooccurrence(hits_long, ["ruleset", "description"]).to_csv(
        output_dir / "cooccurrence_description.csv", index=False
    )
