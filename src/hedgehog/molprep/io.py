from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pandas as pd

from hedgehog.molprep.types import MolPrepFailure


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

    stable_cols = ["smiles", smiles_raw_col, "model_name", "mol_idx"]
    for c in stable_cols:
        if c not in passed_df.columns:
            passed_df[c] = pd.NA
        if c not in failed_df.columns:
            failed_df[c] = pd.NA

    passed_df[stable_cols].to_csv(filtered_path, index=False)
    if failed_df.empty:
        failed_df = pd.DataFrame(columns=[*stable_cols, "reason", "step", "reason_detail"])
    failed_df.to_csv(failed_path, index=False)

    reasons = Counter([f.reason for f in failures])
    metrics = {
        "total_in": int(len(input_df)),
        "passed": int(len(passed_df)),
        "failed": int(len(failures)),
        "dedup_removed": int(dedup_removed),
        "failed_by_reason_json": json.dumps(dict(sorted(reasons.items())), ensure_ascii=False),
    }
    pd.DataFrame([metrics]).to_csv(metrics_path, index=False)

    if write_duplicates_removed and not duplicates_df.empty:
        duplicates_df.to_csv(out_dir / "duplicates_removed.csv", index=False)

    return passed_df, reasons, stable_cols
