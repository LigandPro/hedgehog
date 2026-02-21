"""JavaScript payload builders for report generation."""

from __future__ import annotations

import html as html_lib
import statistics
from pathlib import Path
from typing import Any

import pandas as pd

from hedgehog.reporting import plots
from hedgehog.utils.constants import STAGE_DIRS


def build_descriptors_js_data(
    desc_detailed: dict[str, Any], available_models: list[str]
) -> dict[str, Any]:
    """Build JSON data for JavaScript descriptors visualization.

    Args:
        desc_detailed: Detailed descriptor data from _get_descriptors_detailed
        available_models: List of model names

    Returns:
        Dictionary with structure for JavaScript plotting
    """
    raw_data = desc_detailed.get("raw_data", [])
    if not raw_data:
        return {}
    thresholds = {
        "MolWt": {"min": 100, "max": 500},
        "LogP": {"min": -2, "max": 5},
        "cLogP": {"min": -2, "max": 5},
        "TPSA": {"min": 20, "max": 140},
        "NumHDonors": {"min": 0, "max": 5},
        "NumHAcceptors": {"min": 0, "max": 10},
        "NumRotatableBonds": {"min": 0, "max": 10},
        "n_rot_bonds": {"min": 0, "max": 10},
        "QED": {"min": 0.3, "max": 1},
        "Fsp3": {"min": 0, "max": 1},
        "n_atoms": {"min": 10, "max": 70},
        "n_heavy_atoms": {"min": 10, "max": 50},
        "n_rings": {"min": 1, "max": 6},
        "n_aroma_rings": {"min": 0, "max": 4},
        "n_rigid_bonds": {"min": 0, "max": 30},
    }
    all_descriptors = set()
    for record in raw_data:
        for key in record.keys():
            if key != "model_name":
                all_descriptors.add(key)
    descriptor_names = sorted(list(all_descriptors))
    all_data = {}
    for desc_name in descriptor_names:
        values = [r.get(desc_name) for r in raw_data if r.get(desc_name) is not None]
        if values:
            all_data[desc_name] = {
                "values": values,
                "mean": statistics.mean(values),
                "median": statistics.median(values),
                "std": statistics.stdev(values) if len(values) > 1 else 0,
                "min": min(values),
                "max": max(values),
            }
    by_model = {}
    for model in available_models:
        model_records = [r for r in raw_data if r.get("model_name") == model]
        if not model_records:
            continue
        model_data = {}
        for desc_name in descriptor_names:
            values = [
                r.get(desc_name) for r in model_records if r.get(desc_name) is not None
            ]
            if values:
                model_data[desc_name] = {
                    "values": values,
                    "mean": statistics.mean(values),
                    "median": statistics.median(values),
                    "std": statistics.stdev(values) if len(values) > 1 else 0,
                    "min": min(values),
                    "max": max(values),
                }
        if model_data:
            by_model[model] = model_data
    compare_data = {
        "is_comparison": True,
        "models": available_models,
        "model_colors": plots.COMPARE_PALETTE[: len(available_models)],
        "data": by_model,
    }
    return {
        "all": all_data,
        "by_model": by_model,
        "__compare__": compare_data,
        "models": available_models,
        "descriptors": descriptor_names,
        "thresholds": thresholds,
    }


def build_filters_js_data(
    base_path: Path, available_models: list[str]
) -> dict[str, Any]:
    """Build JSON data for JavaScript filters visualization.

    Args:
        available_models: List of model names

    Returns:
        Dictionary with structure for JavaScript plotting.
        Sub-metrics are converted from ratios to absolute molecule counts.
    """
    result = {"filters": [], "models": available_models, "filter_data": {}}
    for stage_key in ["struct_filters_post"]:
        stage_dir = base_path / STAGE_DIRS.get(stage_key, "")
        if not stage_dir.exists():
            continue
        for subdir in sorted(stage_dir.iterdir()):
            if not subdir.is_dir() or subdir.name in ["plots", "logs"]:
                continue
            filter_name = subdir.name
            if filter_name in result["filter_data"]:
                continue
            metrics_path = subdir / "metrics.csv"
            if not metrics_path.exists():
                continue
            try:
                df = pd.read_csv(metrics_path)
            except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
                continue
            if "model_name" not in df.columns:
                continue
            filter_info = {
                "models": available_models,
                "num_mol": {},
                "pass_rate": {},
                "sub_metrics": {},
                "sub_metric_names": [],
            }
            if "num_mol" in df.columns:
                for model in available_models:
                    model_df = df[df["model_name"] == model]
                    if not model_df.empty:
                        filter_info["num_mol"][model] = int(model_df["num_mol"].iloc[0])
            exclude_cols = {"model_name", "num_mol", "smiles", "mol_idx"}
            all_cols = [c for c in df.columns if c not in exclude_cols]
            ratio_cols = [
                c
                for c in all_cols
                if c.endswith("_banned_ratio")
                or c == "banned_ratio"
                or c.startswith("banned_ratio_")
                or c.endswith("_ratio")
                or (c == "all_banned_ratio")
                or (c == "any_banned_ratio")
            ]
            if ratio_cols:
                sub_metric_names = []
                col_to_name = {}
                for col in ratio_cols:
                    if col == "banned_ratio":
                        name = "failed"
                    elif col == "all_banned_ratio":
                        name = "all"
                    elif col == "any_banned_ratio":
                        continue
                    elif col.startswith("banned_ratio_"):
                        name = col.replace("banned_ratio_", "")
                    elif col.endswith("_banned_ratio"):
                        name = col.replace("_banned_ratio", "")
                    elif col.endswith("_ratio"):
                        name = col.replace("_ratio", "")
                    else:
                        name = col
                    sub_metric_names.append(name)
                    col_to_name[col] = name
                filter_info["sub_metric_names"] = sub_metric_names
                main_ratio_col = None
                for candidate in ["all_banned_ratio", "banned_ratio"]:
                    if candidate in df.columns:
                        main_ratio_col = candidate
                        break
                has_only_failed = sub_metric_names == ["failed"] or (
                    len(sub_metric_names) == 1 and "failed" in sub_metric_names
                )
                if has_only_failed:
                    filter_info["sub_metric_names"] = ["passed", "failed"]
                for model in available_models:
                    model_df = df[df["model_name"] == model]
                    if not model_df.empty:
                        num_mol = int(model_df["num_mol"].iloc[0])
                        model_metrics = {}
                        for col, name in col_to_name.items():
                            if col in model_df.columns:
                                ratio = model_df[col].iloc[0]
                                if pd.notna(ratio):
                                    model_metrics[name] = int(round(num_mol * ratio))
                                else:
                                    model_metrics[name] = 0
                        if has_only_failed and "failed" in model_metrics:
                            model_metrics["passed"] = num_mol - model_metrics["failed"]
                        filter_info["sub_metrics"][model] = model_metrics
                        if main_ratio_col and main_ratio_col in model_df.columns:
                            main_ratio = model_df[main_ratio_col].iloc[0]
                            if pd.notna(main_ratio):
                                filter_info["pass_rate"][model] = round(
                                    (1 - main_ratio) * 100, 1
                                )
                            else:
                                filter_info["pass_rate"][model] = 100.0
                        else:
                            filter_info["pass_rate"][model] = 100.0
            else:
                if "banned_ratio" not in all_cols:
                    continue
                filter_info["sub_metric_names"] = ["passed", "failed"]
                for model in available_models:
                    model_df = df[df["model_name"] == model]
                    if not model_df.empty:
                        num_mol = int(model_df["num_mol"].iloc[0])
                        ratio = model_df["banned_ratio"].iloc[0]
                        failed = int(round(num_mol * ratio)) if pd.notna(ratio) else 0
                        passed = num_mol - failed
                        filter_info["sub_metrics"][model] = {
                            "passed": passed,
                            "failed": failed,
                        }
                        if pd.notna(ratio):
                            filter_info["pass_rate"][model] = round(
                                (1 - ratio) * 100, 1
                            )
                        else:
                            filter_info["pass_rate"][model] = 100.0
            if filter_info["sub_metrics"]:
                result["filters"].append(filter_name)
                result["filter_data"][filter_name] = filter_info
    return result


def build_descriptor_comparison_data(
    initial_detailed: dict[str, Any], final_detailed: dict[str, Any]
) -> dict[str, Any]:
    """Build comparison data between initial and final descriptors.

    Args:
        initial_detailed: Detailed data from initial descriptors stage
        final_detailed: Detailed data from final descriptors stage

    Returns:
        Dictionary with per-descriptor initial/final value arrays and stats
    """
    initial_raw = initial_detailed.get("raw_data", [])
    final_raw = final_detailed.get("raw_data", [])
    if not initial_raw or not final_raw:
        return {}
    initial_descs = set()
    for record in initial_raw:
        initial_descs.update(k for k in record if k != "model_name")
    final_descs = set()
    for record in final_raw:
        final_descs.update(k for k in record if k != "model_name")
    common_descs = sorted(initial_descs & final_descs)
    if not common_descs:
        return {}
    comparison = {}
    for desc in common_descs:
        init_vals = [r[desc] for r in initial_raw if r.get(desc) is not None]
        final_vals = [r[desc] for r in final_raw if r.get(desc) is not None]
        if init_vals and final_vals:
            comparison[desc] = {
                "initial_values": init_vals,
                "final_values": final_vals,
                "mean_initial": statistics.mean(init_vals),
                "mean_final": statistics.mean(final_vals),
            }
    return {"descriptors": list(comparison.keys()), "data": comparison}


def build_model_options(models: list[str]) -> str:
    """Build HTML option elements for model dropdown.

    Args:
        models: List of model names

    Returns:
        HTML string with option elements
    """
    return "".join(
        f'<option value="{html_lib.escape(str(m), quote=True)}">{html_lib.escape(str(m))}</option>'
        for m in models
    )
