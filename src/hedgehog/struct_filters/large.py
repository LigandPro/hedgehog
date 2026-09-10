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
    resolve_large_output,
    should_enable_all_large_filters,
    should_filter_large_outputs,
    stage_output_or_parts,
    write_manifest,
    write_stage_counts,
)
from hedgehog.struct_filters.utils import (
    attach_structural_enforcement_pass,
    build_structural_policy_pass_masks,
    filter_function_applier,
    finalize_structural_liability_profile,
    get_basic_stats,
    initialize_structural_liability_profile,
    merge_structural_liability_profile,
    merge_structural_policy_aliases,
    prepare_structfilters_input,
    process_prepared_payload,
)
from hedgehog.struct_filters.waves.registry import (
    get_aligned_enforced_filters,
    get_calculated_policy_names,
    policy_calculation_filter,
)
from hedgehog.utils.parallel import resolve_n_jobs
from hedgehog.utils.paths import process_path

IDENTITY_COLUMNS = ["smiles", "model_name", "mol_idx"]


def _enabled_filters(
    config_struct_filters: dict[str, Any], config: dict[str, Any] | None = None
) -> dict[str, bool]:
    force_all = should_enable_all_large_filters(config)
    return {
        key.replace("calculate_", ""): value
        for key, value in config_struct_filters.items()
        if key.startswith("calculate_")
        and isinstance(value, bool)
        and (value or force_all)
    }


def _struct_input_chunks(config: dict, stage_dir: str, chunk_rows: int):
    base = Path(process_path(config[KEY_FOLDER_TO_SAVE]))
    is_post_descriptors = (
        "03_structural_filters_post" in stage_dir or stage_dir == "StructFilters"
    )
    if is_post_descriptors:
        descriptors_output = stage_output_or_parts(
            base / "stages" / "02_descriptors_initial" / "filtered_molecules.csv"
        )
        if descriptors_output is None:
            descriptors_output = stage_output_or_parts(
                base
                / "stages"
                / "02_descriptors_initial"
                / "filtered"
                / "filtered_molecules.csv"
            )
        if descriptors_output is not None:
            if descriptors_output.is_file() and descriptors_output.stat().st_size == 0:
                return
            yield from iter_csv_parts(descriptors_output, chunk_rows=chunk_rows)
            return

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


def _merge_named_pass_mask(
    combined: pd.DataFrame, policy_name: str, mask: pd.DataFrame
) -> pd.DataFrame:
    id_cols = [
        c for c in IDENTITY_COLUMNS if c in combined.columns and c in mask.columns
    ]
    if mask.empty or not id_cols:
        combined[policy_name] = False
        return combined
    selected = mask[id_cols + ["pass"]].rename(columns={"pass": policy_name})
    out = combined.merge(selected, on=id_cols, how="left")
    out[policy_name] = out[policy_name].fillna(False).astype(bool)
    return out


def _summarize_filter_parts(filter_name: str, parts_dir: Path) -> dict[str, Any]:
    rows = 0
    passed = 0
    pass_col_counts: dict[str, dict[str, int]] = defaultdict(
        lambda: {"passed": 0, "total": 0}
    )
    numeric: dict[str, dict[str, float]] = {}

    for part_df in iter_csv_parts(parts_dir):
        rows += len(part_df)
        if "pass" in part_df.columns:
            values = part_df["pass"].fillna(False).astype(bool)
            passed += int(values.sum())
        for col in [c for c in part_df.columns if c.startswith("pass_")]:
            values = part_df[col].fillna(False).astype(bool)
            pass_col_counts[col]["passed"] += int(values.sum())
            pass_col_counts[col]["total"] += len(values)
        for col in part_df.columns:
            if col in IDENTITY_COLUMNS or col.startswith("pass"):
                continue
            values = pd.to_numeric(part_df[col], errors="coerce").dropna()
            if values.empty:
                continue
            cur = numeric.setdefault(
                col,
                {
                    "count": 0.0,
                    "sum": 0.0,
                    "min": float(values.min()),
                    "max": float(values.max()),
                },
            )
            cur["count"] += float(len(values))
            cur["sum"] += float(values.sum())
            cur["min"] = min(cur["min"], float(values.min()))
            cur["max"] = max(cur["max"], float(values.max()))

    summary = {
        "filter": filter_name,
        "total": rows,
        "passed": passed,
        "failed": rows - passed,
        "passed_ratio": (passed / rows) if rows else 0.0,
    }
    for col, counts in sorted(pass_col_counts.items()):
        total = counts["total"]
        summary[f"{col}_passed"] = counts["passed"]
        summary[f"{col}_failed"] = total - counts["passed"]
    for col, stats in sorted(numeric.items()):
        count = stats["count"]
        summary[f"{col}_mean"] = (stats["sum"] / count) if count else 0.0
        summary[f"{col}_min"] = stats["min"]
        summary[f"{col}_max"] = stats["max"]
    return summary


def run_large(
    config: dict,
    stage_dir: str,
    config_struct_filters: dict[str, Any],
    output_dir: Path,
    reporter=None,
) -> None:
    """Run structural filters in streaming large-dataset mode."""
    output_dir = Path(output_dir)
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    chunk_rows = get_chunk_rows(config)
    single_csv_limit = get_single_csv_limit(config)

    filters = list(_enabled_filters(config_struct_filters, config))
    policy_filters = get_calculated_policy_names(filters, config_struct_filters)
    policies_by_calculation = {
        filter_name: [
            policy_name
            for policy_name in policy_filters
            if policy_calculation_filter(policy_name) == filter_name
        ]
        for filter_name in filters
    }
    aligned_enforced = get_aligned_enforced_filters(config, config_struct_filters)
    if aligned_enforced is not None:
        missing_enforced = {
            name
            for name in aligned_enforced
            if policy_calculation_filter(name) not in filters
        }
        if missing_enforced:
            names = ", ".join(sorted(missing_enforced))
            raise ValueError(
                f"Structural filter flag enabled but calculation disabled: {names}"
            )
    enforced_filters = (
        policy_filters
        if aligned_enforced is None
        else [name for name in policy_filters if name in aligned_enforced]
    )
    filter_outputs = should_filter_large_outputs(config)
    filter_writers = {
        name: ShardedCsvWriter(output_dir / "intermediate" / name, config)
        for name in filters
    }
    filtered_writer = ShardedCsvWriter(
        parts_dir_for_csv(output_dir / "filtered_molecules.csv"), config
    )
    liability_writer = ShardedCsvWriter(
        parts_dir_for_csv(output_dir / "structural_liability_profile.csv"), config
    )

    manifest_rows: list[dict[str, Any]] = []
    stage_rows: list[dict[str, Any]] = []
    failure_rows: list[dict[str, Any]] = []
    processed_total = 0
    n_jobs = resolve_n_jobs(config_struct_filters, config)

    logger.info(
        "StructFilters large mode: chunk_rows=%d, filters=%s",
        chunk_rows,
        ",".join(filters),
    )

    for chunk_index, chunk in enumerate(
        _struct_input_chunks(config, stage_dir, chunk_rows), start=1
    ):
        if chunk.empty:
            continue
        processed_total += len(chunk)
        if reporter is not None:
            reporter.progress(
                processed_total,
                processed_total,
                message=f"StructFilters chunk {chunk_index}",
            )

        prepared_payload = prepare_structfilters_input(chunk, None, n_jobs)
        id_cols = [c for c in IDENTITY_COLUMNS if c in chunk.columns]
        combined = chunk[id_cols].copy()
        liability_profile = initialize_structural_liability_profile(chunk)

        for filter_name in filters:
            try:
                apply_func = filter_function_applier(filter_name)
                filter_results = process_prepared_payload(
                    config, prepared_payload, apply_func
                )
            except Exception as exc:
                logger.warning(
                    "StructFilters large mode: filter %s failed on chunk %d: %s",
                    filter_name,
                    chunk_index,
                    exc,
                )
                failure_rows.append(
                    {
                        "chunk": chunk_index,
                        "filter": filter_name,
                        "error": str(exc),
                    }
                )
                for policy_name in policies_by_calculation[filter_name]:
                    combined[policy_name] = False
                continue
            if filter_results is None:
                for policy_name in policies_by_calculation[filter_name]:
                    combined[policy_name] = False
                continue
            model_names = sorted(
                chunk["model_name"].dropna().astype(str).unique().tolist()
            )
            model_name = model_names[0] if len(model_names) == 1 else model_names
            _, final_extended = get_basic_stats(
                config_struct_filters,
                filter_results,
                model_name,
                filter_name=filter_name,
            )
            final_extended, _enforcement_mask = attach_structural_enforcement_pass(
                config_struct_filters, filter_name, final_extended
            )
            filter_writers[filter_name].write(final_extended)
            policy_masks = build_structural_policy_pass_masks(
                config_struct_filters,
                filter_name,
                final_extended,
                default_mask=_enforcement_mask,
            )
            for policy_name, policy_mask in policy_masks.items():
                combined = _merge_named_pass_mask(combined, policy_name, policy_mask)
            liability_profile = merge_structural_liability_profile(
                liability_profile, filter_name, final_extended
            )
            liability_profile = merge_structural_policy_aliases(
                liability_profile, filter_name
            )

        if config_struct_filters.get("write_structural_liability_profile", True):
            liability_profile = finalize_structural_liability_profile(
                liability_profile, aligned_enforced, policy_filters
            )
            liability_writer.write(liability_profile)

        if enforced_filters and filter_outputs:
            pass_all = combined[enforced_filters].all(axis=1)
            passed_ids = combined.loc[pass_all, id_cols]
            filtered = chunk.merge(passed_ids, on=id_cols, how="inner")
        else:
            filtered = chunk.copy()
        filtered_writer.write(filtered[id_cols])

        stage_rows.append(
            {
                "chunk": chunk_index,
                "input": len(chunk),
                "filtered": len(filtered),
            }
        )
        manifest_rows.append(
            {
                "table": "struct_filters.filtered_molecules",
                "chunk": chunk_index,
                "rows": len(filtered),
                "parts_dir": str(filtered_writer.parts_dir),
            }
        )
        logger.info(
            "StructFilters large chunk %d: %d input, %d filtered",
            chunk_index,
            len(chunk),
            len(filtered),
        )

    summary_dir = output_dir / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)
    summary_rows = []
    for filter_name, writer in filter_writers.items():
        summary = _summarize_filter_parts(filter_name, writer.parts_dir)
        summary_rows.append(summary)
        pd.DataFrame([summary]).to_csv(
            summary_dir / f"{filter_name}_metrics.csv", index=False
        )
    if summary_rows:
        pd.DataFrame(summary_rows).to_csv(
            summary_dir / "filter_counts.tsv", sep="\t", index=False
        )
    if failure_rows:
        pd.DataFrame(failure_rows).to_csv(
            summary_dir / "filter_failures.tsv", sep="\t", index=False
        )

    write_stage_counts(summary_dir / "stage_counts.tsv", stage_rows)
    write_manifest(output_dir / "manifest.tsv", manifest_rows)
    if config_struct_filters.get("write_structural_liability_profile", True):
        materialize_csv_if_small(
            liability_writer.parts_dir,
            output_dir / "structural_liability_profile.csv",
            single_csv_limit,
            columns=["smiles", "model_name", "mol_idx"],
        )
    materialize_csv_if_small(
        filtered_writer.parts_dir,
        output_dir / "filtered_molecules.csv",
        single_csv_limit,
        columns=["smiles", "model_name", "mol_idx"],
    )

    logger.info(
        "StructFilters large mode complete: %d input, %d filtered",
        processed_total,
        count_part_rows(filtered_writer.parts_dir) or 0,
    )
