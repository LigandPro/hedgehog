from __future__ import annotations

import shutil
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd

from hedgehog._constants import KEY_FOLDER_TO_SAVE
from hedgehog.configs.logger import logger
from hedgehog.descriptors.compute import compute_metrics
from hedgehog.descriptors.filtering import filter_molecules
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
    resolve_large_output,
    should_enable_all_large_filters,
    should_filter_large_outputs,
    write_manifest,
    write_stage_counts,
)
from hedgehog.utils.paths import process_path


def _descriptor_input_chunks(config: dict, chunk_rows: int):
    base = Path(process_path(config[KEY_FOLDER_TO_SAVE]))
    molprep_output = resolve_large_output(
        base, "stages", "01_mol_prep", "filtered_molecules.csv"
    )
    if molprep_output is not None:
        yield from iter_csv_parts(molprep_output, chunk_rows=chunk_rows)
        return

    assigner = StreamingMolIdxAssigner(
        base / "_workdir" / "large_dataset" / "identity.sqlite", base
    )
    try:
        for chunk in iter_input_chunks(config["generated_mols_path"], chunk_rows):
            yield assigner.assign(chunk)
    finally:
        assigner.close()


def _update_numeric_summary(
    summary: dict[str, dict[str, float]], df: pd.DataFrame
) -> None:
    id_cols = {"smiles", "model_name", "mol_idx"}
    for col in df.columns:
        if col in id_cols:
            continue
        values = pd.to_numeric(df[col], errors="coerce").dropna()
        if values.empty:
            continue
        current = summary.setdefault(
            col,
            {
                "count": 0.0,
                "sum": 0.0,
                "min": float(values.min()),
                "max": float(values.max()),
            },
        )
        current["count"] += float(len(values))
        current["sum"] += float(values.sum())
        current["min"] = min(current["min"], float(values.min()))
        current["max"] = max(current["max"], float(values.max()))


def _write_descriptor_summary(
    descriptors_folder: Path,
    stage_rows: list[dict[str, Any]],
    pass_counts: dict[str, dict[str, int]],
    numeric_summary: dict[str, dict[str, float]],
) -> None:
    summary_dir = descriptors_folder / "summary"
    write_stage_counts(summary_dir / "stage_counts.tsv", stage_rows)

    pass_rows = []
    for col, counts in sorted(pass_counts.items()):
        total = counts.get("total", 0)
        passed = counts.get("passed", 0)
        pass_rows.append(
            {
                "filter": col,
                "passed": passed,
                "failed": total - passed,
                "total": total,
                "passed_ratio": (passed / total) if total else 0.0,
            }
        )
    if pass_rows:
        pd.DataFrame(pass_rows).to_csv(
            summary_dir / "descriptor_filter_counts.tsv", sep="\t", index=False
        )

    numeric_rows = []
    for col, stats in sorted(numeric_summary.items()):
        count = stats["count"]
        numeric_rows.append(
            {
                "descriptor": col,
                "count": int(count),
                "mean": (stats["sum"] / count) if count else 0.0,
                "min": stats["min"],
                "max": stats["max"],
            }
        )
    if numeric_rows:
        pd.DataFrame(numeric_rows).to_csv(
            summary_dir / "descriptor_numeric_summary.tsv", sep="\t", index=False
        )


def run_large(
    data,
    config: dict,
    config_descriptors: dict,
    descriptors_folder: Path,
    reporter=None,
) -> pd.DataFrame | None:
    """Compute descriptors in streaming large-dataset mode."""
    del data

    descriptors_folder = Path(descriptors_folder)
    if descriptors_folder.exists():
        shutil.rmtree(descriptors_folder)
    metrics_folder = descriptors_folder / "metrics"
    filtered_folder = descriptors_folder / "filtered"
    metrics_folder.mkdir(parents=True, exist_ok=True)
    filtered_folder.mkdir(parents=True, exist_ok=True)

    chunk_rows = get_chunk_rows(config)
    single_csv_limit = get_single_csv_limit(config)
    work_dir = (
        Path(process_path(config[KEY_FOLDER_TO_SAVE]))
        / "_workdir"
        / "large_dataset"
        / "descriptors"
    )
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    descriptors_writer = ShardedCsvWriter(
        parts_dir_for_csv(metrics_folder / "descriptors_all.csv"), config
    )
    flags_writer = ShardedCsvWriter(
        parts_dir_for_csv(filtered_folder / "pass_flags.csv"), config
    )
    passed_writer = ShardedCsvWriter(
        parts_dir_for_csv(filtered_folder / "descriptors_passed.csv"), config
    )
    filtered_writer = ShardedCsvWriter(
        parts_dir_for_csv(filtered_folder / "filtered_molecules.csv"), config
    )
    failed_writer = ShardedCsvWriter(
        parts_dir_for_csv(filtered_folder / "descriptors_failed.csv"), config
    )
    skipped_writer = ShardedCsvWriter(
        parts_dir_for_csv(metrics_folder / "skipped_molecules.csv"), config
    )

    stage_rows: list[dict[str, Any]] = []
    manifest_rows: list[dict[str, Any]] = []
    pass_counts: dict[str, dict[str, int]] = defaultdict(
        lambda: {"passed": 0, "total": 0}
    )
    numeric_summary: dict[str, dict[str, float]] = {}
    processed_total = 0
    calculate_filters = bool(config_descriptors.get("filter_data", False)) or (
        should_enable_all_large_filters(config)
    )
    filter_outputs = should_filter_large_outputs(config)

    for chunk_index, chunk in enumerate(
        _descriptor_input_chunks(config, chunk_rows), start=1
    ):
        if chunk.empty:
            continue
        processed_total += len(chunk)
        if reporter is not None:
            reporter.progress(
                processed_total,
                processed_total,
                message=f"Descriptors chunk {chunk_index}",
            )

        chunk_dir = work_dir / f"chunk-{chunk_index:06d}"
        chunk_metrics = chunk_dir / "metrics"
        chunk_filtered = chunk_dir / "filtered"
        chunk_metrics.mkdir(parents=True, exist_ok=True)
        chunk_filtered.mkdir(parents=True, exist_ok=True)

        metrics_df = compute_metrics(
            chunk,
            chunk_metrics,
            config=config,
            config_descriptors=config_descriptors,
            reporter=None,
        )
        descriptors_writer.write(metrics_df)
        _update_numeric_summary(numeric_summary, metrics_df)

        skipped_path = chunk_metrics / "skipped_molecules.csv"
        if skipped_path.exists() and skipped_path.stat().st_size > 0:
            skipped_writer.write(pd.read_csv(skipped_path))

        filtered_count = 0
        failed_count = 0
        if calculate_filters:
            filter_molecules(metrics_df, config_descriptors, chunk_filtered)
            outputs = {
                "pass_flags": (chunk_filtered / "pass_flags.csv", flags_writer),
                "descriptors_passed": (
                    chunk_filtered / "descriptors_passed.csv",
                    passed_writer,
                ),
                "filtered_molecules": (
                    chunk_filtered / "filtered_molecules.csv",
                    filtered_writer,
                ),
                "descriptors_failed": (
                    chunk_filtered / "descriptors_failed.csv",
                    failed_writer,
                ),
            }
            for table, (path, writer) in outputs.items():
                if path.exists() and path.stat().st_size > 0:
                    df_out = pd.read_csv(path)
                    if table != "filtered_molecules" or filter_outputs:
                        writer.write(df_out)
                    if table == "filtered_molecules":
                        filtered_count = len(df_out)
                    elif table == "descriptors_failed":
                        failed_count = len(df_out)
                    if table == "pass_flags":
                        for col in [c for c in df_out.columns if c.endswith("_pass")]:
                            values = df_out[col].fillna(False).astype(bool)
                            pass_counts[col]["passed"] += int(values.sum())
                            pass_counts[col]["total"] += int(len(values))
            if not filter_outputs:
                filtered_writer.write(metrics_df[["smiles", "model_name", "mol_idx"]])
                filtered_count = len(metrics_df)
        else:
            filtered_writer.write(metrics_df[["smiles", "model_name", "mol_idx"]])
            filtered_count = len(metrics_df)

        stage_rows.append(
            {
                "chunk": chunk_index,
                "input": len(chunk),
                "descriptors": len(metrics_df),
                "filtered": filtered_count,
                "failed": failed_count,
            }
        )
        manifest_rows.append(
            {
                "table": "descriptors.metrics.descriptors_all",
                "chunk": chunk_index,
                "rows": len(metrics_df),
                "parts_dir": str(descriptors_writer.parts_dir),
            }
        )
        logger.info(
            "Descriptors large chunk %d: %d input, %d descriptors, %d filtered",
            chunk_index,
            len(chunk),
            len(metrics_df),
            filtered_count,
        )

    for parts_dir, output_csv, columns in [
        (descriptors_writer.parts_dir, metrics_folder / "descriptors_all.csv", None),
        (flags_writer.parts_dir, filtered_folder / "pass_flags.csv", None),
        (passed_writer.parts_dir, filtered_folder / "descriptors_passed.csv", None),
        (
            filtered_writer.parts_dir,
            filtered_folder / "filtered_molecules.csv",
            ["smiles", "model_name", "mol_idx"],
        ),
        (failed_writer.parts_dir, filtered_folder / "descriptors_failed.csv", None),
        (
            skipped_writer.parts_dir,
            metrics_folder / "skipped_molecules.csv",
            ["smiles", "model_name", "mol_idx"],
        ),
    ]:
        materialize_csv_if_small(
            parts_dir, output_csv, single_csv_limit, columns=columns
        )

    _write_descriptor_summary(
        descriptors_folder,
        stage_rows,
        pass_counts,
        numeric_summary,
    )
    write_manifest(descriptors_folder / "manifest.tsv", manifest_rows)

    logger.info(
        "Descriptors large mode complete: %d input, %d filtered",
        processed_total,
        count_part_rows(filtered_writer.parts_dir) or 0,
    )
    return None
