from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import pandas as pd

from hedgehog.configs.logger import logger
from hedgehog.molprep.io import write_outputs
from hedgehog.molprep.types import MolPrepFailure
from hedgehog.utils.parallel import parallel_map, resolve_n_jobs


def _resolve_molprep_n_jobs(cfg: dict[str, Any]) -> int:
    default = -1
    env_n_jobs = os.environ.get("MOLSCORE_NJOBS")
    if env_n_jobs:
        try:
            parsed = int(env_n_jobs)
            if parsed > 0:
                default = parsed
        except ValueError:
            default = -1

    if isinstance(cfg, dict) and "n_jobs" in cfg and default > 0:
        try:
            stage_n_jobs = int(cfg["n_jobs"])
        except (TypeError, ValueError):
            stage_n_jobs = None
        if stage_n_jobs is not None and stage_n_jobs <= 0:
            return default

    return resolve_n_jobs(stage_config=cfg, default=default)


def _get_cfg(cfg: dict[str, Any], path: list[str], default: Any) -> Any:
    cur: Any = cfg
    for key in path:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def run_mol_prep(
    df: pd.DataFrame,
    cfg: dict[str, Any],
    out_dir: Path,
    reporter=None,
) -> pd.DataFrame:
    """Run MolPrep on a dataframe and write stage outputs."""
    from hedgehog.molprep.workers import (
        _init_molprep_worker,
        _is_missing_smiles_value,
        _process_molprep_item,
    )

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cols_cfg = _get_cfg(cfg, ["columns"], {}) or {}
    smiles_col = str(cols_cfg.get("smiles", "smiles"))
    model_col = str(cols_cfg.get("model_name", "model_name"))
    idx_col = str(cols_cfg.get("mol_idx", "mol_idx"))
    smiles_raw_col = str(cols_cfg.get("smiles_raw", "smiles_raw"))

    items: list[tuple[str, str | None, str | None]] = []
    for _, row in df.iterrows():
        smiles_value = row.get(smiles_col, "")
        if _is_missing_smiles_value(smiles_value):
            smiles_raw = ""
        else:
            smiles_raw = str(smiles_value).strip()
        model_name = row.get(model_col)
        mol_idx = row.get(idx_col)
        items.append(
            (
                smiles_raw,
                str(model_name) if model_name is not None else None,
                str(mol_idx) if mol_idx is not None else None,
            )
        )

    n_jobs = _resolve_molprep_n_jobs(cfg)
    logger.info("MolPrep workers: %d", n_jobs)
    progress_cb = None
    if reporter is not None:

        def _progress_cb(done: int, total: int) -> None:
            reporter.progress(
                done,
                total,
                message="Mol Prep",
                molecules_in=len(df),
            )

        progress_cb = _progress_cb

    results = parallel_map(
        _process_molprep_item,
        items,
        n_jobs=n_jobs,
        progress=progress_cb,
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
                "smiles": f.smiles_raw,  # Keep compatibility: 'smiles' column present
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

    # Deduplicate after standardization (within model_name)
    dedup_removed = 0
    duplicates_df = pd.DataFrame()
    if not passed_df.empty:
        before = len(passed_df)
        dup_mask = passed_df.duplicated(subset=["smiles", "model_name"], keep="first")
        duplicates_df = passed_df[dup_mask].copy() if dup_mask.any() else pd.DataFrame()
        passed_df = passed_df[~dup_mask].reset_index(drop=True)
        dedup_removed = before - len(passed_df)

    passed_df, reasons, stable_cols = write_outputs(
        out_dir=out_dir,
        input_df=df,
        passed_df=passed_df,
        failed_df=failed_df,
        failures=failures,
        dedup_removed=dedup_removed,
        duplicates_df=duplicates_df,
        smiles_raw_col=smiles_raw_col,
        write_duplicates_removed=bool(
            _get_cfg(cfg, ["output", "write_duplicates_removed"], True)
        ),
    )

    logger.info(
        "MolPrep: %d in, %d passed, %d failed, %d dedup removed",
        len(df),
        len(passed_df),
        len(failures),
        dedup_removed,
    )

    # Also log top failure reasons (human-friendly)
    if reasons:
        top = reasons.most_common(10)
        logger.info("MolPrep failure reasons (top): %s", json.dumps(top))

    return passed_df[stable_cols]
