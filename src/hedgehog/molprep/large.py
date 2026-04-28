from __future__ import annotations

import json
import shutil
import sqlite3
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd

from hedgehog._constants import KEY_FOLDER_TO_SAVE
from hedgehog.configs.logger import logger
from hedgehog.large_dataset import (
    ShardedCsvWriter,
    StreamingMolIdxAssigner,
    count_part_rows,
    get_chunk_rows,
    get_single_csv_limit,
    iter_input_chunks,
    materialize_csv_if_small,
    parts_dir_for_csv,
    write_manifest,
    write_stage_counts,
)
from hedgehog.molprep.orchestrator import _get_cfg, _resolve_molprep_n_jobs
from hedgehog.molprep.types import MolPrepFailure
from hedgehog.utils.parallel import parallel_map
from hedgehog.utils.paths import process_path


def _ensure_dedup_table(conn: sqlite3.Connection) -> None:
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS molprep_seen (
            model_name TEXT NOT NULL,
            smiles TEXT NOT NULL,
            PRIMARY KEY(model_name, smiles)
        )
        """
    )
    conn.commit()


def _deduplicate_global(
    conn: sqlite3.Connection, passed_df: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if passed_df.empty:
        return passed_df, pd.DataFrame(columns=passed_df.columns)

    keep_indices: list[int] = []
    duplicate_indices: list[int] = []
    for idx, row in passed_df.iterrows():
        model_name = str(row.get("model_name") or "")
        smiles = str(row.get("smiles") or "")
        cur = conn.execute(
            "INSERT OR IGNORE INTO molprep_seen(model_name, smiles) VALUES (?, ?)",
            (model_name, smiles),
        )
        if cur.rowcount == 1:
            keep_indices.append(idx)
        else:
            duplicate_indices.append(idx)
    conn.commit()

    kept = passed_df.loc[keep_indices].reset_index(drop=True)
    duplicates = passed_df.loc[duplicate_indices].reset_index(drop=True)
    return kept, duplicates


def _write_metrics(out_dir: Path, rows: list[dict[str, Any]]) -> None:
    totals = Counter()
    reason_counts = Counter()
    for row in rows:
        totals.update(
            {
                "input": int(row["input"]),
                "passed": int(row["passed"]),
                "failed": int(row["failed"]),
                "duplicates_removed": int(row["duplicates_removed"]),
            }
        )
        if row.get("failure_reasons_json"):
            reason_counts.update(dict(json.loads(str(row["failure_reasons_json"]))))

    metrics = pd.DataFrame(
        [
            {
                "input_molecules": totals["input"],
                "passed_molecules": totals["passed"],
                "failed_molecules": totals["failed"],
                "duplicates_removed": totals["duplicates_removed"],
            }
        ]
    )
    metrics.to_csv(out_dir / "metrics.csv", index=False)
    write_stage_counts(out_dir / "summary" / "stage_counts.tsv", rows)
    if reason_counts:
        pd.DataFrame(reason_counts.most_common(), columns=["reason", "count"]).to_csv(
            out_dir / "summary" / "failure_reasons.tsv", sep="\t", index=False
        )


def run_large(config: dict, cfg: dict, out_dir: Path, reporter=None) -> None:
    """Run MolPrep in streaming large-dataset mode."""
    from hedgehog.molprep.workers import (
        _init_molprep_worker,
        _is_missing_smiles_value,
        _process_molprep_item,
    )

    out_dir = Path(out_dir)
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    run_base = Path(process_path(config[KEY_FOLDER_TO_SAVE]))
    work_dir = run_base / "_workdir" / "large_dataset"
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    chunk_rows = get_chunk_rows(config)
    single_csv_limit = get_single_csv_limit(config)
    input_path = config["generated_mols_path"]

    assigner = StreamingMolIdxAssigner(work_dir / "identity.sqlite", run_base)
    dedup_conn = sqlite3.connect(work_dir / "state.sqlite")
    _ensure_dedup_table(dedup_conn)

    passed_writer = ShardedCsvWriter(parts_dir_for_csv(out_dir / "filtered_molecules.csv"), config)
    failed_writer = ShardedCsvWriter(parts_dir_for_csv(out_dir / "failed_molecules.csv"), config)
    duplicates_writer = ShardedCsvWriter(
        parts_dir_for_csv(out_dir / "duplicates_removed.csv"), config
    )

    cols_cfg = _get_cfg(cfg, ["columns"], {}) or {}
    smiles_col = str(cols_cfg.get("smiles", "smiles"))
    model_col = str(cols_cfg.get("model_name", "model_name"))
    idx_col = str(cols_cfg.get("mol_idx", "mol_idx"))
    smiles_raw_col = str(cols_cfg.get("smiles_raw", "smiles_raw"))
    n_jobs = _resolve_molprep_n_jobs(cfg)
    logger.info("MolPrep large mode: chunk_rows=%d, workers=%d", chunk_rows, n_jobs)

    manifest_rows: list[dict[str, Any]] = []
    stage_rows: list[dict[str, Any]] = []
    total_seen = 0

    try:
        for chunk_index, chunk in enumerate(iter_input_chunks(input_path, chunk_rows), start=1):
            chunk = assigner.assign(chunk)
            total_seen += len(chunk)
            if reporter is not None:
                reporter.progress(
                    total_seen,
                    total_seen,
                    message=f"Mol Prep chunk {chunk_index}",
                    molecules_in=total_seen,
                )

            items: list[tuple[str, str | None, str | None]] = []
            for _, row in chunk.iterrows():
                smiles_value = row.get(smiles_col, "")
                smiles_raw = "" if _is_missing_smiles_value(smiles_value) else str(smiles_value).strip()
                model_name = row.get(model_col)
                mol_idx = row.get(idx_col)
                items.append(
                    (
                        smiles_raw,
                        str(model_name) if model_name is not None else None,
                        str(mol_idx) if mol_idx is not None else None,
                    )
                )

            results = parallel_map(
                _process_molprep_item,
                items,
                n_jobs=n_jobs,
                initializer=_init_molprep_worker,
                initargs=(cfg,),
            )

            passed_rows: list[dict[str, Any]] = []
            failures: list[MolPrepFailure] = []
            for passed, failure in results:
                if passed is not None:
                    passed_rows.append(passed)
                elif failure is not None:
                    failures.append(failure)

            passed_df = pd.DataFrame(passed_rows)
            if not passed_df.empty and smiles_raw_col != "smiles_raw":
                passed_df[smiles_raw_col] = passed_df["smiles_raw"]

            failed_df = pd.DataFrame(
                [
                    {
                        "smiles": f.smiles_raw,
                        "smiles_raw": f.smiles_raw,
                        "model_name": f.model_name,
                        "mol_idx": f.mol_idx,
                        "reason": f.reason,
                        "step": f.step,
                        "reason_detail": f.reason_detail,
                    }
                    for f in failures
                ]
            )
            if not failed_df.empty and smiles_raw_col != "smiles_raw":
                failed_df[smiles_raw_col] = failed_df["smiles_raw"]

            passed_df, duplicates_df = _deduplicate_global(dedup_conn, passed_df)
            stable_cols = ["smiles", "model_name", "mol_idx"]
            extra_cols = [c for c in ["smiles_raw", smiles_raw_col] if c in passed_df.columns]
            stable_cols = [*stable_cols, *[c for c in extra_cols if c not in stable_cols]]

            passed_writer.write(passed_df[[c for c in stable_cols if c in passed_df.columns]])
            failed_writer.write(failed_df)
            duplicates_writer.write(duplicates_df)

            failure_reasons = Counter(f.reason for f in failures)
            row = {
                "chunk": chunk_index,
                "input": len(chunk),
                "passed": len(passed_df),
                "failed": len(failed_df),
                "duplicates_removed": len(duplicates_df),
                "failure_reasons_json": json.dumps(dict(failure_reasons), sort_keys=True),
            }
            stage_rows.append(row)
            manifest_rows.append(
                {
                    "table": "mol_prep.filtered_molecules",
                    "chunk": chunk_index,
                    "rows": len(passed_df),
                    "parts_dir": str(passed_writer.parts_dir),
                }
            )
            logger.info(
                "MolPrep large chunk %d: %d in, %d passed, %d failed, %d dedup removed",
                chunk_index,
                len(chunk),
                len(passed_df),
                len(failed_df),
                len(duplicates_df),
            )
    finally:
        assigner.close()
        dedup_conn.close()

    materialize_csv_if_small(
        passed_writer.parts_dir,
        out_dir / "filtered_molecules.csv",
        single_csv_limit,
        columns=["smiles", "model_name", "mol_idx", "smiles_raw"],
    )
    materialize_csv_if_small(
        failed_writer.parts_dir,
        out_dir / "failed_molecules.csv",
        single_csv_limit,
        columns=[
            "smiles",
            "smiles_raw",
            "model_name",
            "mol_idx",
            "reason",
            "step",
            "reason_detail",
        ],
    )
    materialize_csv_if_small(
        duplicates_writer.parts_dir,
        out_dir / "duplicates_removed.csv",
        single_csv_limit,
        columns=["smiles", "model_name", "mol_idx", "smiles_raw"],
    )
    _write_metrics(out_dir, stage_rows)
    write_manifest(out_dir / "manifest.tsv", manifest_rows)

    logger.info(
        "MolPrep large mode complete: %d input, %d passed",
        total_seen,
        count_part_rows(passed_writer.parts_dir) or 0,
    )
