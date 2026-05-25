from __future__ import annotations

import shutil
from collections import defaultdict
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
    iter_csv_parts,
    iter_input_chunks,
    materialize_csv_if_small,
    parts_dir_for_csv,
    should_filter_large_outputs,
    stage_output_or_parts,
    write_manifest,
    write_stage_counts,
)
from hedgehog.synthesis.utils import (
    _build_score_filter_mask,
    calculate_synthesis_scores,
)
from hedgehog.utils.paths import process_path

IDENTITY_COLUMNS = ["smiles", "model_name", "mol_idx"]
SCORE_FILTERS = [
    ("sa_score", "sa_score_min", "sa_score_max"),
    ("ra_score", "ra_score_min", "ra_score_max"),
    ("syba_score", "syba_score_min", "syba_score_max"),
]


def _synthesis_input_chunks(config: dict, chunk_rows: int):
    base = Path(process_path(config[KEY_FOLDER_TO_SAVE]))
    candidates = [
        base / "stages" / "03_structural_filters_post" / "filtered_molecules.csv",
        base
        / "stages"
        / "02_descriptors_initial"
        / "filtered"
        / "filtered_molecules.csv",
        base / "stages" / "01_mol_prep" / "filtered_molecules.csv",
    ]
    for candidate in candidates:
        source = stage_output_or_parts(candidate)
        if source is not None:
            yield from iter_csv_parts(source, chunk_rows=chunk_rows)
            return

    assigner = StreamingMolIdxAssigner(
        base / "_workdir" / "large_dataset" / "identity.sqlite", base
    )
    try:
        for chunk in iter_input_chunks(config["generated_mols_path"], chunk_rows):
            yield assigner.assign(chunk)
    finally:
        assigner.close()


def _score_pass_flags(
    scored_df: pd.DataFrame, config_synthesis: dict[str, Any]
) -> tuple[pd.DataFrame, pd.Series]:
    id_cols = [col for col in IDENTITY_COLUMNS if col in scored_df.columns]
    flags = scored_df[id_cols].copy()
    pass_mask = pd.Series(True, index=scored_df.index)

    for column, min_key, max_key in SCORE_FILTERS:
        mask = _build_score_filter_mask(
            scored_df,
            column,
            config_synthesis.get(min_key, 0),
            config_synthesis.get(max_key, "inf"),
        )
        if mask is None:
            column_pass = pd.Series(True, index=scored_df.index)
        else:
            column_pass = mask.fillna(True).astype(bool)
        flags[f"{column}_pass"] = column_pass.to_numpy()
        pass_mask = pass_mask & column_pass

    flags["synthesis_score_pass"] = pass_mask.to_numpy()
    return flags, pass_mask


def _update_pass_counts(
    pass_counts: dict[str, dict[str, int]], flags: pd.DataFrame
) -> None:
    for col in [c for c in flags.columns if c.endswith("_pass")]:
        values = flags[col].fillna(False).astype(bool)
        pass_counts[col]["passed"] += int(values.sum())
        pass_counts[col]["total"] += int(len(values))


def _write_summary(
    output_folder: Path,
    stage_rows: list[dict[str, Any]],
    pass_counts: dict[str, dict[str, int]],
) -> None:
    summary_dir = output_folder / "summary"
    write_stage_counts(summary_dir / "stage_counts.tsv", stage_rows)
    count_rows = []
    for col, counts in sorted(pass_counts.items()):
        total = counts["total"]
        passed = counts["passed"]
        count_rows.append(
            {
                "filter": col,
                "passed": passed,
                "failed": total - passed,
                "total": total,
                "passed_ratio": (passed / total) if total else 0.0,
            }
        )
    if count_rows:
        pd.DataFrame(count_rows).to_csv(
            summary_dir / "synthesis_score_counts.tsv", sep="\t", index=False
        )


def run_large(
    config: dict,
    config_synthesis: dict[str, Any],
    output_folder: Path,
    reporter=None,
) -> None:
    """Run synthesis scores in large-dataset compute-only mode."""
    output_folder = Path(output_folder)
    if output_folder.exists():
        shutil.rmtree(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)
    base = Path(process_path(config[KEY_FOLDER_TO_SAVE]))
    work_dir = base / "_workdir" / "large_dataset" / "synthesis"
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    chunk_rows = get_chunk_rows(config)
    single_csv_limit = get_single_csv_limit(config)
    filter_outputs = should_filter_large_outputs(config)

    scores_writer = ShardedCsvWriter(
        parts_dir_for_csv(output_folder / "synthesis_scores.csv"), config
    )
    flags_writer = ShardedCsvWriter(
        parts_dir_for_csv(output_folder / "pass_flags.csv"), config
    )
    passed_writer = ShardedCsvWriter(
        parts_dir_for_csv(output_folder / "synthesis_passed.csv"), config
    )
    failed_writer = ShardedCsvWriter(
        parts_dir_for_csv(output_folder / "synthesis_failed.csv"), config
    )
    filtered_writer = ShardedCsvWriter(
        parts_dir_for_csv(output_folder / "filtered_molecules.csv"), config
    )

    stage_rows: list[dict[str, Any]] = []
    manifest_rows: list[dict[str, Any]] = []
    pass_counts: dict[str, dict[str, int]] = defaultdict(
        lambda: {"passed": 0, "total": 0}
    )
    processed_total = 0

    logger.info(
        "Synthesis large mode: chunk_rows=%d, retrosynthesis=disabled",
        chunk_rows,
    )

    for chunk_index, chunk in enumerate(
        _synthesis_input_chunks(config, chunk_rows), start=1
    ):
        if chunk.empty:
            continue
        processed_total += len(chunk)
        if reporter is not None:
            reporter.progress(
                processed_total,
                processed_total,
                message=f"Synthesis chunk {chunk_index}",
            )

        scored_df = calculate_synthesis_scores(
            chunk,
            str(base),
            config_synthesis,
            progress_cb=None,
        )
        flags, pass_mask = _score_pass_flags(scored_df, config_synthesis)
        passed_df = scored_df.loc[pass_mask].copy()
        failed_df = scored_df.loc[~pass_mask].copy()
        filtered_df = passed_df if filter_outputs else scored_df

        scores_writer.write(scored_df)
        flags_writer.write(flags)
        passed_writer.write(passed_df)
        failed_writer.write(failed_df)
        filtered_writer.write(
            filtered_df[[c for c in IDENTITY_COLUMNS if c in filtered_df.columns]]
        )
        _update_pass_counts(pass_counts, flags)

        stage_rows.append(
            {
                "chunk": chunk_index,
                "input": len(chunk),
                "scored": len(scored_df),
                "passed": len(passed_df),
                "failed": len(failed_df),
                "filtered": len(filtered_df),
            }
        )
        manifest_rows.append(
            {
                "table": "synthesis.filtered_molecules",
                "chunk": chunk_index,
                "rows": len(filtered_df),
                "parts_dir": str(filtered_writer.parts_dir),
            }
        )
        logger.info(
            "Synthesis large chunk %d: %d input, %d scored, %d output",
            chunk_index,
            len(chunk),
            len(scored_df),
            len(filtered_df),
        )

    for parts_dir, output_csv, columns in [
        (scores_writer.parts_dir, output_folder / "synthesis_scores.csv", None),
        (flags_writer.parts_dir, output_folder / "pass_flags.csv", None),
        (passed_writer.parts_dir, output_folder / "synthesis_passed.csv", None),
        (failed_writer.parts_dir, output_folder / "synthesis_failed.csv", None),
        (
            filtered_writer.parts_dir,
            output_folder / "filtered_molecules.csv",
            IDENTITY_COLUMNS,
        ),
    ]:
        materialize_csv_if_small(
            parts_dir, output_csv, single_csv_limit, columns=columns
        )

    _write_summary(output_folder, stage_rows, pass_counts)
    write_manifest(output_folder / "manifest.tsv", manifest_rows)

    logger.info(
        "Synthesis large mode complete: %d input, %d output",
        processed_total,
        count_part_rows(filtered_writer.parts_dir) or 0,
    )
