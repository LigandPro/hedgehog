from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pandas as pd

from hedgehog.molprep.types import MolPrepFailure


def _build_molprep_detail(
    *,
    input_df: pd.DataFrame,
    passed_df: pd.DataFrame,
    failed_df: pd.DataFrame,
    duplicates_df: pd.DataFrame,
    smiles_raw_col: str,
) -> pd.DataFrame:
    """Build a per-molecule MolPrep detail table with step/filter pass flags."""
    detail = pd.DataFrame()
    detail["model_name"] = (
        input_df["model_name"] if "model_name" in input_df.columns else pd.NA
    )
    detail["mol_idx"] = input_df["mol_idx"] if "mol_idx" in input_df.columns else pd.NA
    if smiles_raw_col in input_df.columns:
        detail["smiles_raw"] = input_df[smiles_raw_col]
    elif "smiles" in input_df.columns:
        detail["smiles_raw"] = input_df["smiles"]
    else:
        detail["smiles_raw"] = pd.NA
    detail["smiles"] = pd.NA

    key_cols = ["model_name", "mol_idx", "smiles_raw"]
    for key in key_cols:
        if key not in detail.columns:
            detail[key] = pd.NA

    success_union = pd.concat(
        [passed_df.copy(), duplicates_df.copy()], ignore_index=True
    )
    if not success_union.empty:
        success_cols = {
            "smiles": "smiles",
            "model_name": "model_name",
            "mol_idx": "mol_idx",
            smiles_raw_col: "smiles_raw",
        }
        success_src = pd.DataFrame()
        for src_col, dst_col in success_cols.items():
            if src_col in success_union.columns:
                success_src[dst_col] = success_union[src_col]
        if all(col in success_src.columns for col in key_cols):
            success_src = success_src.drop_duplicates(subset=key_cols, keep="first")
            detail = detail.merge(success_src, on=key_cols, how="left", suffixes=("", "_ok"))
            if "smiles_ok" in detail.columns:
                detail["smiles"] = detail["smiles_ok"].combine_first(detail["smiles"])
                detail = detail.drop(columns=["smiles_ok"])

    fail_src = pd.DataFrame()
    for src_col, dst_col in (
        ("model_name", "model_name"),
        ("mol_idx", "mol_idx"),
        (smiles_raw_col, "smiles_raw"),
        ("reason", "fail_reason"),
        ("step", "fail_step"),
        ("reason_detail", "fail_reason_detail"),
    ):
        if src_col in failed_df.columns:
            fail_src[dst_col] = failed_df[src_col]
    if not fail_src.empty and all(col in fail_src.columns for col in key_cols):
        fail_src = fail_src.drop_duplicates(subset=key_cols, keep="first")
        detail = detail.merge(fail_src, on=key_cols, how="left")
    else:
        detail["fail_reason"] = pd.NA
        detail["fail_step"] = pd.NA
        detail["fail_reason_detail"] = pd.NA

    passed_keys = set()
    dup_keys = set()
    if not passed_df.empty and all(col in passed_df.columns for col in ["model_name", "mol_idx", smiles_raw_col]):
        passed_keys = set(
            zip(
                passed_df["model_name"].astype(str),
                passed_df["mol_idx"].astype(str),
                passed_df[smiles_raw_col].astype(str),
                strict=False,
            )
        )
    if not duplicates_df.empty and all(col in duplicates_df.columns for col in ["model_name", "mol_idx", smiles_raw_col]):
        dup_keys = set(
            zip(
                duplicates_df["model_name"].astype(str),
                duplicates_df["mol_idx"].astype(str),
                duplicates_df[smiles_raw_col].astype(str),
                strict=False,
            )
        )

    key_tuples = list(
        zip(
            detail["model_name"].astype(str),
            detail["mol_idx"].astype(str),
            detail["smiles_raw"].astype(str),
            strict=False,
        )
    )
    detail["pre_dedup_pass"] = [k in passed_keys or k in dup_keys for k in key_tuples]
    detail["dedup_pass"] = pd.NA
    detail.loc[[k in passed_keys for k in key_tuples], "dedup_pass"] = True
    detail.loc[[k in dup_keys for k in key_tuples], "dedup_pass"] = False
    detail["stage_pass"] = [k in passed_keys for k in key_tuples]

    # Step-level pass flags. False when that step is the failure point.
    step_fail_reason = detail["fail_reason"].astype(str)
    detail["to_mol_pass"] = step_fail_reason != "parse_failed"
    detail["fix_mol_pass"] = step_fail_reason != "fix_failed"
    detail["sanitize_mol_pass"] = step_fail_reason != "sanitize_failed"
    detail["remove_salts_solvents_pass"] = step_fail_reason != "remove_salts_failed"
    detail["keep_largest_fragment_pass"] = ~step_fail_reason.isin(
        {"largest_fragment_failed", "multifragment_after_largest"}
    )
    detail["standardize_mol_pass"] = step_fail_reason != "standardize_mol_failed"
    detail["remove_stereochemistry_pass"] = step_fail_reason != "remove_stereo_failed"
    detail["to_smiles_pass"] = step_fail_reason != "to_smiles_failed"
    detail["standardize_smiles_pass"] = ~step_fail_reason.isin(
        {"standardize_smiles_failed", "post_standardize_parse_failed"}
    )

    filters_reached = detail["pre_dedup_pass"] | (detail["fail_step"].astype(str) == "filters")
    for col in (
        "filter_allowed_atoms_pass",
        "filter_radicals_pass",
        "filter_isotopes_pass",
        "filter_single_fragment_pass",
    ):
        detail[col] = pd.NA
        detail.loc[filters_reached, col] = True
    detail.loc[detail["fail_reason"].astype(str) == "disallowed_atoms", "filter_allowed_atoms_pass"] = False
    detail.loc[detail["fail_reason"].astype(str) == "radicals", "filter_radicals_pass"] = False
    detail.loc[detail["fail_reason"].astype(str) == "isotopes", "filter_isotopes_pass"] = False
    detail.loc[detail["fail_reason"].astype(str) == "multifragment", "filter_single_fragment_pass"] = False

    # If reason_detail stores independent filter outcomes (JSON), use them.
    for idx, raw in detail["fail_reason_detail"].items():
        if pd.isna(raw):
            continue
        try:
            payload = json.loads(str(raw))
        except Exception:
            continue
        if not isinstance(payload, dict):
            continue
        for col in (
            "filter_allowed_atoms_pass",
            "filter_radicals_pass",
            "filter_isotopes_pass",
            "filter_single_fragment_pass",
        ):
            if col in payload:
                val = payload[col]
                if isinstance(val, bool):
                    detail.at[idx, col] = val

    # Explicit aliases for analysis tables.
    detail["valid"] = detail["pre_dedup_pass"]
    detail["unique"] = detail["dedup_pass"]
    detail["allowed_atoms_ok"] = detail["filter_allowed_atoms_pass"]
    detail["single_fragment_ok"] = detail["filter_single_fragment_pass"]
    detail["radicals_ok"] = detail["filter_radicals_pass"]
    detail["isotopes_ok"] = detail["filter_isotopes_pass"]

    ordered_cols = [
        "smiles",
        "smiles_raw",
        "model_name",
        "mol_idx",
        "to_mol_pass",
        "fix_mol_pass",
        "sanitize_mol_pass",
        "remove_salts_solvents_pass",
        "keep_largest_fragment_pass",
        "standardize_mol_pass",
        "remove_stereochemistry_pass",
        "to_smiles_pass",
        "standardize_smiles_pass",
        "filter_allowed_atoms_pass",
        "filter_radicals_pass",
        "filter_isotopes_pass",
        "filter_single_fragment_pass",
        "pre_dedup_pass",
        "dedup_pass",
        "stage_pass",
        "valid",
        "unique",
        "allowed_atoms_ok",
        "single_fragment_ok",
        "radicals_ok",
        "isotopes_ok",
        "fail_reason",
        "fail_step",
        "fail_reason_detail",
    ]
    for col in ordered_cols:
        if col not in detail.columns:
            detail[col] = pd.NA
    return detail[ordered_cols]


def write_outputs(
    *,
    out_dir: Path,
    input_df: pd.DataFrame,
    passed_df: pd.DataFrame,
    failed_df: pd.DataFrame,
    failures: list[MolPrepFailure],
    dedup_removed: int,
    duplicates_df: pd.DataFrame,
    smiles_raw_col: str,
    write_duplicates_removed: bool,
) -> tuple[pd.DataFrame, Counter[str], list[str]]:
    """Write MolPrep output files and return prepared data for final return/logging."""
    filtered_path = out_dir / "filtered_molecules.csv"
    failed_path = out_dir / "failed_molecules.csv"
    metrics_path = out_dir / "metrics.csv"
    detail_path = out_dir / "molprep_detail.csv"

    stable_cols = ["smiles", smiles_raw_col, "model_name", "mol_idx"]
    for c in stable_cols:
        if c not in passed_df.columns:
            passed_df[c] = pd.NA
        if c not in failed_df.columns:
            failed_df[c] = pd.NA

    passed_df[stable_cols].to_csv(filtered_path, index=False)
    if failed_df.empty:
        failed_df = pd.DataFrame(
            columns=[*stable_cols, "reason", "step", "reason_detail"]
        )
    failed_df.to_csv(failed_path, index=False)

    reasons = Counter([f.reason for f in failures])
    metrics = {
        "total_in": int(len(input_df)),
        "passed": int(len(passed_df)),
        "failed": int(len(failures)),
        "dedup_removed": int(dedup_removed),
        "failed_by_reason_json": json.dumps(
            dict(sorted(reasons.items())), ensure_ascii=False
        ),
    }
    pd.DataFrame([metrics]).to_csv(metrics_path, index=False)

    if write_duplicates_removed and not duplicates_df.empty:
        duplicates_df.to_csv(out_dir / "duplicates_removed.csv", index=False)

    detail_df = _build_molprep_detail(
        input_df=input_df,
        passed_df=passed_df,
        failed_df=failed_df,
        duplicates_df=duplicates_df,
        smiles_raw_col=smiles_raw_col,
    )
    detail_df.to_csv(detail_path, index=False)

    return passed_df, reasons, stable_cols
