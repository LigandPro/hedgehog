"""Data collection helpers for report generation."""

from __future__ import annotations

import base64
import json
import logging
import statistics
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd

from hedgehog.reporting import moleval_metrics
from hedgehog.utils.constants import STAGE_DIRS, STAGE_DISPLAY_NAMES

logger = logging.getLogger(__name__)

# Key descriptors to show in report
KEY_DESCRIPTORS = [
    "MolWt",
    "LogP",
    "TPSA",
    "NumHDonors",
    "NumHAcceptors",
    "NumRotatableBonds",
]

# Mapping for case-insensitive descriptor lookup (lowercase -> display name)
DESCRIPTOR_ALIASES = {
    "molwt": "MolWt",
    "logp": "LogP",
    "clogp": "cLogP",
    "tpsa": "TPSA",
    "numhdonors": "NumHDonors",
    "numhacceptors": "NumHAcceptors",
    "numrotatablebonds": "NumRotatableBonds",
    "hbd": "NumHDonors",
    "hba": "NumHAcceptors",
    "qed": "QED",
    "fsp3": "Fsp3",
    "n_atoms": "n_atoms",
    "n_heavy_atoms": "n_heavy_atoms",
    "n_het_atoms": "n_het_atoms",
    "n_n_atoms": "n_N_atoms",
    "fn_atoms": "fN_atoms",
    "charged_mol": "charged_mol",
    "sw": "SW",
    "ring_size": "ring_size",
    "n_rings": "n_rings",
    "n_aroma_rings": "n_aroma_rings",
    "n_fused_aromatic_rings": "n_fused_aromatic_rings",
    "n_rigid_bonds": "n_rigid_bonds",
    "n_rot_bonds": "n_rot_bonds",
}

# Columns to exclude from descriptor analysis (not numeric descriptors)
DESCRIPTOR_EXCLUDE_COLS = {
    "smiles",
    "model_name",
    "mol_idx",
    "chars",
    "name",
    "id",
}


def collect_all_data(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
) -> dict[str, Any]:
    """Collect all metrics from stage outputs.

    Returns:
        Dictionary with all collected data
    """
    funnel_data = get_all_funnel_data(
        base_path, config, stages, initial_count, final_count, output_dir
    )
    return {
        "metadata": get_metadata(
            base_path, config, stages, initial_count, final_count, output_dir
        ),
        "summary": get_summary(
            base_path, config, stages, initial_count, final_count, output_dir
        ),
        "funnel": funnel_data["all"],
        "funnel_by_model": funnel_data["by_model"],
        "available_models": funnel_data["models"],
        "stages": get_stage_stats(
            base_path, config, stages, initial_count, final_count, output_dir
        ),
        "models": get_model_stats(
            base_path, config, stages, initial_count, final_count, output_dir
        ),
        "descriptors": get_descriptor_stats(
            base_path, config, stages, initial_count, final_count, output_dir
        ),
        "descriptors_detailed": get_descriptors_detailed(
            base_path, config, stages, initial_count, final_count, output_dir
        ),
        "filters": get_filter_stats(
            base_path, config, stages, initial_count, final_count, output_dir
        ),
        "filters_detailed": get_filters_detailed(
            base_path, config, stages, initial_count, final_count, output_dir
        ),
        "synthesis": get_synthesis_stats(
            base_path, config, stages, initial_count, final_count, output_dir
        ),
        "synthesis_detailed": get_synthesis_detailed(
            base_path, config, stages, initial_count, final_count, output_dir
        ),
        "retrosynthesis": get_retrosynthesis_detailed(
            base_path, config, stages, initial_count, final_count, output_dir
        ),
        "docking": get_docking_stats(
            base_path, config, stages, initial_count, final_count, output_dir
        ),
        "docking_detailed": get_docking_detailed(
            base_path, config, stages, initial_count, final_count, output_dir
        ),
        "docking_filters_detailed": get_docking_filters_detailed(
            base_path, config, stages, initial_count, final_count, output_dir
        ),
        "descriptors_final": get_descriptor_stats(
            base_path,
            config,
            stages,
            initial_count,
            final_count,
            output_dir,
            "descriptors_final",
        ),
        "descriptors_final_detailed": get_descriptors_detailed(
            base_path,
            config,
            stages,
            initial_count,
            final_count,
            output_dir,
            "descriptors_final",
        ),
        "existing_plots": get_existing_plots(
            base_path, config, stages, initial_count, final_count, output_dir
        ),
        "moleval": get_moleval_metrics(
            base_path, config, stages, initial_count, final_count, output_dir
        ),
        "config": get_config_summary(
            base_path, config, stages, initial_count, final_count, output_dir
        ),
    }


def get_metadata(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
) -> dict[str, Any]:
    """Get report metadata."""
    return {
        "generated_at": datetime.now().isoformat(),
        "hedgehog_version": "1.0.0",
        "run_path": str(base_path),
    }


def get_summary(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
) -> dict[str, Any]:
    """Get executive summary statistics."""
    retention_rate = final_count / initial_count if initial_count > 0 else 0
    if stages:
        enabled_stages = [s for s in stages if s.enabled]
        completed_stages = [s for s in stages if s.completed]
        stage_statuses = [
            {
                "name": s.name,
                "enabled": s.enabled,
                "completed": s.completed,
                "status": "completed"
                if s.completed
                else "failed"
                if s.enabled
                else "disabled",
            }
            for s in stages
        ]
    else:
        stage_order = [
            ("mol_prep", "Mol Prep"),
            ("descriptors_initial", "Initial Descriptors"),
            ("struct_filters_post", "Structural Filters"),
            ("synthesis", "Synthesis Analysis"),
            ("docking", "Molecular Docking"),
            ("docking_filters", "Docking Filters"),
            ("descriptors_final", "Final Descriptors"),
        ]
        enabled_stages = []
        completed_stages = []
        stage_statuses = []
        for stage_key, display_name in stage_order:
            stage_dir = base_path / STAGE_DIRS.get(stage_key, "")
            if stage_dir.exists():
                enabled_stages.append(stage_key)
                has_output = any(stage_dir.glob("*.csv")) or any(
                    stage_dir.glob("*/*.csv")
                )
                if has_output:
                    completed_stages.append(stage_key)
                stage_statuses.append(
                    {
                        "name": display_name,
                        "enabled": True,
                        "completed": has_output,
                        "status": "completed" if has_output else "failed",
                    }
                )
    return {
        "initial_molecules": initial_count,
        "final_molecules": final_count,
        "retention_rate": retention_rate,
        "retention_percent": f"{retention_rate * 100:.2f}%",
        "stages_enabled": len(enabled_stages),
        "stages_completed": len(completed_stages),
        "stage_statuses": stage_statuses,
    }


def get_funnel_data(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
) -> list[dict[str, Any]]:
    """Get molecule funnel data through pipeline filtering stages.

    Only includes stages that actually filter molecules:
    Descriptors, Structural Filters, Synthesis, Docking Filters.
    Starts from the raw initial count.
    """
    funnel = [{"stage": "Initial", "count": initial_count}]
    stage_order = [
        ("mol_prep", "Mol Prep"),
        ("descriptors_initial", "Descriptors"),
        ("struct_filters_post", "Structural Filters"),
        ("synthesis", "Synthesis"),
        ("docking_filters", "Docking Filters"),
    ]
    for stage_key, display_name in stage_order:
        count = get_stage_output_count(
            base_path, config, stages, initial_count, final_count, output_dir, stage_key
        )
        if count is not None:
            funnel.append({"stage": display_name, "count": count})
    return funnel


def get_stage_output_count(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
    stage_key: str,
) -> int | None:
    """Get molecule count from stage output file."""
    stage_dir = STAGE_DIRS.get(stage_key)
    if not stage_dir:
        return None
    paths_to_try = [
        base_path / stage_dir / "filtered_molecules.csv",
        base_path / stage_dir / "filtered" / "filtered_molecules.csv",
        base_path / stage_dir / "ligands.csv",
    ]
    for path in paths_to_try:
        if path.exists():
            try:
                df = pd.read_csv(path)
                return len(df)
            except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as e:
                logger.debug("Could not read %s: %s", path, e)
                continue
    return None


def get_stage_output_count_by_model(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
    stage_key: str,
    model_name: str | None = None,
) -> int | None:
    """Get molecule count from stage output file, optionally filtered by model.

    Args:
        stage_key: Stage key from STAGE_DIRS
        model_name: Optional model name to filter by

    Returns:
        Count of molecules, or None if data unavailable
    """
    stage_dir = STAGE_DIRS.get(stage_key)
    if not stage_dir:
        return None
    paths_to_try = [
        base_path / stage_dir / "filtered_molecules.csv",
        base_path / stage_dir / "filtered" / "filtered_molecules.csv",
        base_path / stage_dir / "ligands.csv",
    ]
    for path in paths_to_try:
        if path.exists():
            try:
                df = pd.read_csv(path)
                if model_name and "model_name" in df.columns:
                    df = df[df["model_name"] == model_name]
                return len(df)
            except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as e:
                logger.debug("Could not read %s: %s", path, e)
                continue
    return None


def get_available_models(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
) -> list[str]:
    """Get list of available model names from input or output files.

    Returns:
        List of unique model names
    """
    models = set()
    input_path = base_path / "input" / "sampled_molecules.csv"
    if input_path.exists():
        try:
            df = pd.read_csv(input_path)
            if "model_name" in df.columns:
                models.update(df["model_name"].dropna().unique())
        except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as e:
            logger.debug("Could not read %s: %s", input_path, e)
    output_path = output_dir / "final_molecules.csv"
    if output_path.exists():
        try:
            df = pd.read_csv(output_path)
            if "model_name" in df.columns:
                models.update(df["model_name"].dropna().unique())
        except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as e:
            logger.debug("Could not read %s: %s", output_path, e)
    return sorted(models)


def get_initial_count_by_model(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
    model_name: str | None = None,
) -> int:
    """Get initial molecule count, optionally filtered by model.

    Args:
        model_name: Optional model name to filter by

    Returns:
        Count of initial molecules
    """
    if not model_name:
        return initial_count
    input_path = base_path / "input" / "sampled_molecules.csv"
    if input_path.exists():
        try:
            df = pd.read_csv(input_path)
            if "model_name" in df.columns:
                return len(df[df["model_name"] == model_name])
        except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as e:
            logger.debug("Could not read %s: %s", input_path, e)
    return 0


def get_funnel_data_by_model(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
    model_name: str | None = None,
) -> list[dict[str, Any]]:
    """Get molecule funnel data for a specific model or all models.

    Args:
        model_name: Model name to filter by, or None for all models

    Returns:
        List of funnel stage data
    """
    initial_count = get_initial_count_by_model(
        base_path, config, stages, initial_count, final_count, output_dir, model_name
    )
    funnel = [{"stage": "Initial", "count": initial_count}]
    stage_order = [
        ("mol_prep", "Mol Prep"),
        ("descriptors_initial", "Descriptors"),
        ("struct_filters_post", "Structural Filters"),
        ("synthesis", "Synthesis"),
        ("docking_filters", "Docking Filters"),
    ]
    for stage_key, display_name in stage_order:
        count = get_stage_output_count_by_model(
            base_path,
            config,
            stages,
            initial_count,
            final_count,
            output_dir,
            stage_key,
            model_name,
        )
        if count is not None:
            funnel.append({"stage": display_name, "count": count})
    return funnel


def get_all_funnel_data(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
) -> dict[str, Any]:
    """Get funnel data for all models and per-model breakdown.

    Returns:
        Dictionary with 'all' funnel and 'by_model' funnel data
    """
    result = {
        "all": get_funnel_data(
            base_path, config, stages, initial_count, final_count, output_dir
        ),
        "by_model": {},
        "models": get_available_models(
            base_path, config, stages, initial_count, final_count, output_dir
        ),
    }
    for model in result["models"]:
        result["by_model"][model] = get_funnel_data_by_model(
            base_path, config, stages, initial_count, final_count, output_dir, model
        )
    return result


def get_stage_stats(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
) -> list[dict[str, Any]]:
    """Get pass/fail statistics for each stage."""
    stats = []
    for stage_key, display_name in STAGE_DISPLAY_NAMES.items():
        stage_dir = base_path / STAGE_DIRS.get(stage_key, "")
        if not stage_dir.exists():
            continue
        passed = 0
        failed = 0
        passed_path = stage_dir / "filtered_molecules.csv"
        if not passed_path.exists():
            passed_path = stage_dir / "filtered" / "filtered_molecules.csv"
        failed_path = stage_dir / "failed_molecules.csv"
        if not failed_path.exists():
            failed_path = stage_dir / "filtered" / "failed_molecules.csv"
        if passed_path.exists():
            try:
                passed = len(pd.read_csv(passed_path))
            except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
                pass
        if failed_path.exists():
            try:
                failed = len(pd.read_csv(failed_path))
            except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
                pass
        if passed > 0 or failed > 0:
            stats.append(
                {
                    "stage": display_name,
                    "passed": passed,
                    "failed": failed,
                    "total": passed + failed,
                }
            )
    return stats


def get_model_stats(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
) -> list[dict[str, Any]]:
    """Get per-model statistics."""
    final_path = output_dir / "final_molecules.csv"
    if not final_path.exists():
        return []
    try:
        df = pd.read_csv(final_path)
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
        return []
    if "model_name" not in df.columns:
        return []
    model_stats = []
    for model in df["model_name"].unique():
        final_count = len(df[df["model_name"] == model])
        initial_count = get_initial_model_count(
            base_path, config, stages, initial_count, final_count, output_dir, model
        )
        model_stats.append(
            {
                "model_name": model,
                "initial": initial_count or final_count,
                "final": final_count,
                "losses": get_model_losses(
                    base_path,
                    config,
                    stages,
                    initial_count,
                    final_count,
                    output_dir,
                    model,
                ),
            }
        )
    return model_stats


def get_initial_model_count(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
    model: str,
) -> int | None:
    """Get initial molecule count for a model."""
    input_path = base_path / "input" / "sampled_molecules.csv"
    if not input_path.exists():
        return None
    try:
        df = pd.read_csv(input_path)
        if "model_name" in df.columns:
            return len(df[df["model_name"] == model])
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
        pass
    return None


def get_model_losses(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
    model: str,
) -> dict[str, int]:
    """Get molecule losses per stage for a model."""
    losses = {}
    stage_pairs = [
        ("mol_prep", "mol_prep"),
        ("descriptors", "descriptors_initial"),
        ("struct_filters", "struct_filters_post"),
        ("synthesis", "synthesis"),
        ("docking", "docking"),
        ("docking_filters", "docking_filters"),
    ]
    for loss_key, stage_key in stage_pairs:
        failed_path = base_path / STAGE_DIRS.get(stage_key, "") / "failed_molecules.csv"
        if not failed_path.exists():
            failed_path = (
                base_path
                / STAGE_DIRS.get(stage_key, "")
                / "filtered"
                / "failed_molecules.csv"
            )
        if failed_path.exists():
            try:
                df = pd.read_csv(failed_path)
                if "model_name" in df.columns:
                    losses[loss_key] = len(df[df["model_name"] == model])
                else:
                    losses[loss_key] = 0
            except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
                losses[loss_key] = 0
        else:
            losses[loss_key] = 0
    return losses


def get_descriptor_stats(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
    stage_key: str = "descriptors_initial",
) -> dict[str, Any]:
    """Get descriptor statistics.

    Args:
        stage_key: Stage key from STAGE_DIRS (default: descriptors_initial)
    """
    desc_dir = base_path / STAGE_DIRS[stage_key]
    paths_to_try = [
        desc_dir / "metrics" / "descriptors_all.csv",
        desc_dir / "filtered" / "descriptors_passed.csv",
    ]
    df = None
    for path in paths_to_try:
        if path.exists():
            try:
                df = pd.read_csv(path)
                break
            except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
                continue
    if df is None:
        return {}
    stats = {"distributions": {}, "summary": {}}
    col_map = build_descriptor_column_map(
        base_path, config, stages, initial_count, final_count, output_dir, df.columns
    )
    for desc in KEY_DESCRIPTORS:
        col = col_map.get(desc)
        if col and col in df.columns:
            values = df[col].dropna().tolist()
            if values:
                stats["distributions"][desc] = values
                stats["summary"][desc] = {
                    "mean": float(df[col].mean()),
                    "std": float(df[col].std()),
                    "min": float(df[col].min()),
                    "max": float(df[col].max()),
                }
    return stats


def build_descriptor_column_map(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
    columns: list[str],
) -> dict[str, str]:
    """Build mapping from standard descriptor names to actual column names.

    Args:
        columns: List of column names in the dataframe

    Returns:
        Dict mapping standard name -> actual column name
    """
    col_map = {}
    col_lower = {c.lower(): c for c in columns}
    for desc in KEY_DESCRIPTORS + ["QED", "Fsp3"]:
        if desc in columns:
            col_map[desc] = desc
        elif desc.lower() in col_lower:
            col_map[desc] = col_lower[desc.lower()]
    for alias, standard in DESCRIPTOR_ALIASES.items():
        if alias in col_lower and standard not in col_map:
            col_map[standard] = col_lower[alias]
    return col_map


def get_filter_stats(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
) -> dict[str, Any]:
    """Get filter statistics."""
    filter_stats = {"by_filter": {}, "totals": {}, "by_model": {}}
    for stage_key in ["struct_filters_post"]:
        stage_dir = base_path / STAGE_DIRS.get(stage_key, "")
        if not stage_dir.exists():
            continue
        for subdir in stage_dir.iterdir():
            if not subdir.is_dir():
                continue
            if subdir.name in ["plots", "logs"]:
                continue
            metrics_path = subdir / "metrics.csv"
            if metrics_path.exists():
                try:
                    df = pd.read_csv(metrics_path)
                    if "failed_count" in df.columns:
                        total_failed = df["failed_count"].sum()
                        filter_stats["totals"][subdir.name] = int(total_failed)
                        if "model_name" in df.columns:
                            model_counts = (
                                df.groupby("model_name")["failed_count"].sum().to_dict()
                            )
                            filter_stats["by_filter"][subdir.name] = {
                                str(k): int(v) for k, v in model_counts.items()
                            }
                except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
                    continue
    return filter_stats


def get_synthesis_stats(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
) -> dict[str, Any]:
    """Get synthesis analysis statistics."""
    synth_dir = base_path / STAGE_DIRS["synthesis"]
    scores_path = synth_dir / "synthesis_scores.csv"
    if not scores_path.exists():
        return {}
    try:
        df = pd.read_csv(scores_path)
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
        return {}
    stats = {"distributions": {}, "scatter_data": {}}
    score_columns = ["sa_score", "syba_score", "ra_score", "sc_score"]
    for col in score_columns:
        if col in df.columns:
            values = df[col].dropna().tolist()
            if values:
                stats["distributions"][col] = values
    if "sa_score" in df.columns and "syba_score" in df.columns:
        stats["scatter_data"] = {
            "sa_scores": df["sa_score"].dropna().tolist(),
            "syba_scores": df["syba_score"].dropna().tolist(),
            "model_names": df["model_name"].tolist()
            if "model_name" in df.columns
            else [],
        }
    return stats


def get_model_name_lookup(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
) -> dict[str, str]:
    """Build a lookup table from molecule name to model name.

    Returns:
        Dictionary mapping molecule name/ID to model name
    """
    lookup = {}
    ligands_csv = base_path / STAGE_DIRS["docking"] / "ligands.csv"
    if ligands_csv.exists():
        try:
            df = pd.read_csv(ligands_csv)
            if "model_name" in df.columns:
                for id_col in ["name", "mol_idx", "molecule_id"]:
                    if id_col in df.columns:
                        for _, row in df.iterrows():
                            if pd.notna(row.get(id_col)) and pd.notna(
                                row.get("model_name")
                            ):
                                lookup[str(row[id_col])] = row["model_name"]
                        break
        except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
            pass
    return lookup


def detect_docking_tools(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
) -> list[str]:
    """Detect available docking tools from directory structure.

    Scans the docking directory for subdirectories and SDF files to
    identify which docking tools were used.

    Returns:
        List of tool names (e.g., ["gnina", "smina"])
    """
    docking_dir = base_path / STAGE_DIRS["docking"]
    if not docking_dir.exists():
        return []
    tools = set()
    for subdir in docking_dir.iterdir():
        if subdir.is_dir():
            for sdf_file in subdir.glob("*.sdf"):
                if "out" in sdf_file.name.lower():
                    tools.add(subdir.name)
                    break
    for sdf_file in docking_dir.glob("*_out.sdf"):
        tool_name = sdf_file.stem.replace("_out", "")
        tools.add(tool_name)
    for config_file in docking_dir.glob("*_config.ini"):
        tool_name = config_file.stem.replace("_config", "")
        tools.add(tool_name)
    return sorted(tools)


def parse_docking_sdf(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
    sdf_path: Path,
    model_lookup: dict[str, str] | None = None,
) -> list[dict[str, Any]]:
    """Parse SDF file to extract docking scores and molecule info.

    Args:
        sdf_path: Path to SDF file with docking results
        model_lookup: Optional dict mapping molecule names to model names

    Returns:
        List of dicts with molecule_id, affinity, model_name, and any extra scores
    """
    try:
        from rdkit import Chem
    except ImportError:
        logger.warning("RDKit not available for SDF parsing")
        return []
    if not sdf_path.exists():
        return []
    results = []
    mol_best_scores: dict[str, dict] = {}
    try:
        supplier = Chem.SDMolSupplier(str(sdf_path))
        for mol in supplier:
            if mol is None:
                continue
            mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else None
            if not mol_name:
                continue
            affinity = None
            for prop_name in ["minimizedAffinity", "affinity", "score"]:
                if mol.HasProp(prop_name):
                    try:
                        affinity = float(mol.GetProp(prop_name))
                        break
                    except (ValueError, TypeError):
                        pass
            if affinity is None:
                continue
            model_name = "Unknown"
            if mol.HasProp("s_sm_model_name"):
                model_name = mol.GetProp("s_sm_model_name")
            elif model_lookup and mol_name in model_lookup:
                model_name = model_lookup[mol_name]
            extra_scores = {}
            for prop_name in ["CNNscore", "CNNaffinity", "CNN_VS"]:
                if mol.HasProp(prop_name):
                    try:
                        extra_scores[prop_name] = float(mol.GetProp(prop_name))
                    except (ValueError, TypeError):
                        pass
            record = {
                "molecule_id": mol_name,
                "affinity": affinity,
                "model_name": model_name,
                **extra_scores,
            }
            if mol_name not in mol_best_scores:
                mol_best_scores[mol_name] = record
            elif affinity < mol_best_scores[mol_name]["affinity"]:
                mol_best_scores[mol_name] = record
        results = list(mol_best_scores.values())
    except (OSError, ValueError, TypeError, RuntimeError) as e:
        logger.warning("Error parsing SDF file %s: %s", sdf_path, e)
    return results


def get_docking_stats(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
) -> dict[str, Any]:
    """Get docking statistics."""
    docking_dir = base_path / STAGE_DIRS["docking"]
    detected = detect_docking_tools(
        base_path, config, stages, initial_count, final_count, output_dir
    )
    tools = sorted(set(detected) | {"gnina", "smina"})
    stats = {tool: {} for tool in tools}
    model_lookup = get_model_name_lookup(
        base_path, config, stages, initial_count, final_count, output_dir
    )
    for tool in tools:
        csv_paths = [
            docking_dir / tool / "scores.csv",
            docking_dir / f"{tool}_scores.csv",
            docking_dir / tool / "docking_results.csv",
        ]
        csv_found = False
        for csv_path in csv_paths:
            if csv_path.exists():
                try:
                    df = pd.read_csv(csv_path)
                    for score_col in ["affinity", "score", "minimizedAffinity"]:
                        if score_col in df.columns:
                            stats[tool] = {
                                "scores": df[score_col].dropna().tolist(),
                                "count": len(df),
                            }
                            csv_found = True
                            break
                except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
                    pass
                if csv_found:
                    break
        if not csv_found:
            sdf_paths = [
                docking_dir / tool / f"{tool}_out.sdf",
                docking_dir / tool / "output.sdf",
                docking_dir / f"{tool}_out.sdf",
            ]
            for sdf_path in sdf_paths:
                if sdf_path.exists():
                    records = parse_docking_sdf(
                        base_path,
                        config,
                        stages,
                        initial_count,
                        final_count,
                        output_dir,
                        sdf_path,
                        model_lookup,
                    )
                    if records:
                        scores = [r["affinity"] for r in records]
                        stats[tool] = {"scores": scores, "count": len(records)}
                        break
    return stats


def get_descriptors_detailed(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
    stage_key: str = "descriptors_initial",
) -> dict[str, Any]:
    """Get detailed descriptor data for enhanced visualization.

    Reads ALL numeric columns from the CSV file, not just predefined ones.

    Args:
        stage_key: Stage key from STAGE_DIRS (default: descriptors_initial)

    Returns:
        Dictionary with raw_data, summary_by_model, and key_descriptors
    """
    desc_dir = base_path / STAGE_DIRS[stage_key]
    paths_to_try = [
        desc_dir / "metrics" / "descriptors_all.csv",
        desc_dir / "filtered" / "descriptors_passed.csv",
        desc_dir / "filtered_molecules.csv",
    ]
    df = None
    for path in paths_to_try:
        if path.exists():
            try:
                df = pd.read_csv(path)
                break
            except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
                continue
    if df is None:
        return {}
    result = {
        "raw_data": [],
        "summary_by_model": {},
        "key_descriptors": ["MolWt", "LogP", "TPSA", "QED"],
    }
    col_normalize = {}
    for col in df.columns:
        col_lower = col.lower()
        if col_lower in DESCRIPTOR_EXCLUDE_COLS:
            continue
        if col_lower in DESCRIPTOR_ALIASES:
            col_normalize[col] = DESCRIPTOR_ALIASES[col_lower]
        else:
            col_normalize[col] = col
    numeric_cols = []
    for col in col_normalize.keys():
        if df[col].dtype in ["int64", "float64", "int32", "float32"]:
            numeric_cols.append(col)
    for _, row in df.iterrows():
        record = {}
        if "model_name" in df.columns and pd.notna(row.get("model_name")):
            record["model_name"] = row.get("model_name")
        for actual_col in numeric_cols:
            display_name = col_normalize[actual_col]
            if pd.notna(row.get(actual_col)):
                val = row.get(actual_col)
                if isinstance(val, (int, float)):
                    record[display_name] = float(val)
        if record:
            result["raw_data"].append(record)
    if "model_name" in df.columns:
        for model in df["model_name"].dropna().unique():
            model_df = df[df["model_name"] == model]
            model_summary = {}
            for actual_col in numeric_cols:
                display_name = col_normalize[actual_col]
                values = model_df[actual_col].dropna()
                if len(values) > 0:
                    model_summary[display_name] = float(values.mean())
            if model_summary:
                result["summary_by_model"][model] = model_summary
    return result


def get_filters_detailed(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
) -> dict[str, Any]:
    """Get detailed filter data for enhanced visualization.

    Returns:
        Dictionary with by_filter, banned_ratios, common_alerts_reasons
    """
    result = {
        "by_filter": {},
        "banned_ratios": {},
        "common_alerts_reasons": {},
        "filter_metrics": {},
    }
    for stage_key in ["struct_filters_post"]:
        stage_dir = base_path / STAGE_DIRS.get(stage_key, "")
        if not stage_dir.exists():
            continue
        for subdir in stage_dir.iterdir():
            if not subdir.is_dir() or subdir.name in ["plots", "logs"]:
                continue
            filter_name = subdir.name
            metrics_path = subdir / "metrics.csv"
            if metrics_path.exists():
                try:
                    df = pd.read_csv(metrics_path)
                    if "model_name" in df.columns:
                        if "failed_count" in df.columns:
                            model_counts = (
                                df.groupby("model_name")["failed_count"].sum().to_dict()
                            )
                            result["by_filter"][filter_name] = {
                                str(k): int(v) for k, v in model_counts.items()
                            }
                        if "banned_ratio" in df.columns:
                            ratios = (
                                df.groupby("model_name")["banned_ratio"]
                                .mean()
                                .to_dict()
                            )
                            result["banned_ratios"][filter_name] = {
                                str(k): float(v) for k, v in ratios.items()
                            }
                        elif (
                            "total_count" in df.columns and "failed_count" in df.columns
                        ):
                            model_ratios = {}
                            for model in df["model_name"].unique():
                                model_df = df[df["model_name"] == model]
                                total = model_df["total_count"].sum()
                                failed = model_df["failed_count"].sum()
                                if total > 0:
                                    model_ratios[str(model)] = failed / total
                            if model_ratios:
                                result["banned_ratios"][filter_name] = model_ratios
                    filter_metrics = {}
                    if filter_name == "lilly" and "demerit_score" in df.columns:
                        filter_metrics["avg_demerit_score"] = float(
                            df["demerit_score"].mean()
                        )
                    if filter_name == "molcomplexity":
                        if "bertz" in df.columns:
                            filter_metrics["avg_bertz"] = float(df["bertz"].mean())
                        if "sas" in df.columns:
                            filter_metrics["avg_sas"] = float(df["sas"].mean())
                    if filter_metrics:
                        result["filter_metrics"][filter_name] = filter_metrics
                except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
                    continue
            if filter_name == "common_alerts":
                alerts_path = subdir / "alerts_summary.csv"
                if not alerts_path.exists():
                    alerts_path = subdir / "failed_molecules.csv"
                if alerts_path.exists():
                    try:
                        alerts_df = pd.read_csv(alerts_path)
                        alert_cols = [
                            "alert_type",
                            "reason",
                            "alert_name",
                            "filter_reason",
                        ]
                        for col in alert_cols:
                            if col in alerts_df.columns:
                                reasons = alerts_df[col].value_counts().to_dict()
                                result["common_alerts_reasons"] = {
                                    str(k): int(v) for k, v in reasons.items()
                                }
                                break
                    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
                        pass
    return result


def get_synthesis_detailed(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
) -> dict[str, Any]:
    """Get detailed synthesis data for enhanced visualization.

    Returns:
        Dictionary with raw_data, solved/unsolved counts, time stats, by_model
    """
    synth_dir = base_path / STAGE_DIRS["synthesis"]
    scores_path = synth_dir / "synthesis_extended.csv"
    if not scores_path.exists():
        scores_path = synth_dir / "synthesis_scores.csv"
    if not scores_path.exists():
        return {}
    try:
        df = pd.read_csv(scores_path)
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
        return {}
    result = {
        "raw_data": [],
        "sa_scores": [],
        "syba_scores": [],
        "ra_scores": [],
        "solved_count": 0,
        "unsolved_count": 0,
        "summary": {},
        "by_model": {},
    }
    if "sa_score" in df.columns:
        result["sa_scores"] = df["sa_score"].dropna().tolist()
        if result["sa_scores"]:
            result["summary"]["avg_sa_score"] = float(df["sa_score"].mean())
    if "syba_score" in df.columns:
        result["syba_scores"] = df["syba_score"].dropna().tolist()
        if result["syba_scores"]:
            result["summary"]["avg_syba_score"] = float(df["syba_score"].mean())
    if "ra_score" in df.columns:
        result["ra_scores"] = df["ra_score"].dropna().tolist()
        if result["ra_scores"]:
            result["summary"]["avg_ra_score"] = float(df["ra_score"].dropna().mean())
    solved_col = None
    solved_cols = ["solved", "route_found", "synthesis_solved"]
    for col in solved_cols:
        if col in df.columns:
            solved_col = col
            solved = (
                df[col].sum() if df[col].dtype == bool else df[col].astype(bool).sum()
            )
            result["solved_count"] = int(solved)
            result["unsolved_count"] = len(df) - int(solved)
            result["summary"]["pct_solved"] = (
                100 * solved / len(df) if len(df) > 0 else 0
            )
            break
    time_col = None
    time_cols = ["search_time", "route_time", "synthesis_time"]
    for col in time_cols:
        if col in df.columns:
            times = df[col].dropna().tolist()
            if times:
                time_col = col
                result["summary"]["avg_search_time"] = float(df[col].mean())
                if "model_name" in df.columns:
                    for _, row in df.iterrows():
                        if pd.notna(row.get(col)):
                            result["raw_data"].append(
                                {
                                    "model_name": row.get("model_name", "Unknown"),
                                    "search_time": row.get(col),
                                }
                            )
            break
    if "model_name" in df.columns:
        for model in df["model_name"].dropna().unique():
            model_df = df[df["model_name"] == model]
            model_data = {
                "sa_scores": model_df["sa_score"].dropna().tolist()
                if "sa_score" in df.columns
                else [],
                "syba_scores": model_df["syba_score"].dropna().tolist()
                if "syba_score" in df.columns
                else [],
                "ra_scores": model_df["ra_score"].dropna().tolist()
                if "ra_score" in df.columns
                else [],
                "solved_count": 0,
                "unsolved_count": 0,
                "summary": {},
            }
            if "sa_score" in model_df.columns and len(model_data["sa_scores"]) > 0:
                model_data["summary"]["avg_sa_score"] = float(
                    model_df["sa_score"].mean()
                )
            if "syba_score" in model_df.columns and len(model_data["syba_scores"]) > 0:
                model_data["summary"]["avg_syba_score"] = float(
                    model_df["syba_score"].mean()
                )
            if "ra_score" in model_df.columns and len(model_data["ra_scores"]) > 0:
                model_data["summary"]["avg_ra_score"] = float(
                    model_df["ra_score"].dropna().mean()
                )
            if solved_col and solved_col in model_df.columns:
                model_solved = (
                    model_df[solved_col].sum()
                    if model_df[solved_col].dtype == bool
                    else model_df[solved_col].astype(bool).sum()
                )
                model_data["solved_count"] = int(model_solved)
                model_data["unsolved_count"] = len(model_df) - int(model_solved)
                model_data["summary"]["pct_solved"] = (
                    100 * model_solved / len(model_df) if len(model_df) > 0 else 0
                )
            if time_col and time_col in model_df.columns:
                model_times = model_df[time_col].dropna()
                if len(model_times) > 0:
                    model_data["summary"]["avg_search_time"] = float(model_times.mean())
            result["by_model"][model] = model_data
    return result


def get_retrosynthesis_detailed(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
) -> dict[str, Any]:
    """Get detailed retrosynthesis data from AiZynthFinder results.

    Returns:
        Dictionary with route_scores, steps, precursors, solve_rate, summary
    """
    synth_dir = base_path / STAGE_DIRS["synthesis"]
    json_path = synth_dir / "retrosynthesis_results.json"
    if not json_path.exists():
        return {}
    try:
        with open(json_path) as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError):
        return {}
    if "data" not in data:
        return {}
    result = {
        "route_scores": [],
        "steps": [],
        "precursors": [],
        "solved_count": 0,
        "total_count": 0,
        "summary": {},
    }
    for item in data["data"]:
        result["total_count"] += 1
        if item.get("is_solved", False):
            result["solved_count"] += 1
            if item.get("top_score") is not None:
                result["route_scores"].append(item["top_score"])
            if item.get("number_of_steps") is not None:
                result["steps"].append(item["number_of_steps"])
            if item.get("number_of_precursors") is not None:
                result["precursors"].append(item["number_of_precursors"])
    if result["total_count"] > 0:
        result["summary"]["solve_rate"] = (
            100.0 * result["solved_count"] / result["total_count"]
        )
    if result["route_scores"]:
        result["summary"]["avg_route_score"] = statistics.mean(result["route_scores"])
    if result["steps"]:
        result["summary"]["avg_steps"] = statistics.mean(result["steps"])
    if result["precursors"]:
        result["summary"]["avg_precursors"] = statistics.mean(result["precursors"])
    return result


def get_existing_plots(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
) -> dict[str, str]:
    """Collect existing plot images from stage directories.

    Returns:
        Dictionary mapping plot name -> base64-encoded image data URI
    """
    plots_found = {}
    plot_locations = {
        "descriptors_initial_distribution": STAGE_DIRS["descriptors_initial"]
        + "/plots/descriptors_distribution.png",
        "descriptors_final_distribution": STAGE_DIRS["descriptors_final"]
        + "/plots/descriptors_distribution.png",
        "filters_post_counts": STAGE_DIRS["struct_filters_post"]
        + "/plots/molecule_counts_comparison.png",
        "filters_post_ratios": STAGE_DIRS["struct_filters_post"]
        + "/plots/restriction_ratios_comparison.png",
    }
    for name, rel_path in plot_locations.items():
        full_path = base_path / rel_path
        if full_path.exists():
            try:
                with open(full_path, "rb") as f:
                    img_data = base64.b64encode(f.read()).decode("utf-8")
                plots_found[name] = f"data:image/png;base64,{img_data}"
                logger.debug("Found plot: %s", name)
            except OSError as e:
                logger.warning("Failed to load plot %s: %s", name, e)
    return plots_found


def get_docking_detailed(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
) -> dict[str, Any]:
    """Get detailed docking data for enhanced visualization.

    Supports both CSV and SDF file formats, with automatic tool detection.

    Returns:
        Dictionary with raw_data, top_molecules, summary stats, by_model
    """
    docking_dir = base_path / STAGE_DIRS["docking"]
    detected = detect_docking_tools(
        base_path, config, stages, initial_count, final_count, output_dir
    )
    tools = sorted(set(detected) | {"gnina", "smina"})
    result = {}
    for tool in tools:
        result[tool] = {
            "raw_data": [],
            "top_molecules": [],
            "summary": {},
            "by_model": {},
        }
    model_lookup = get_model_name_lookup(
        base_path, config, stages, initial_count, final_count, output_dir
    )
    for tool in tools:
        records = []
        csv_paths = [
            docking_dir / tool / "scores.csv",
            docking_dir / tool / "docking_results.csv",
            docking_dir / f"{tool}_scores.csv",
        ]
        df = None
        for path in csv_paths:
            if path.exists():
                try:
                    df = pd.read_csv(path)
                    break
                except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
                    continue
        if df is not None:
            score_col = None
            for col in ["affinity", "score", "minimizedAffinity", "docking_score"]:
                if col in df.columns:
                    score_col = col
                    break
            if score_col:
                id_col = None
                for col in ["molecule_id", "mol_id", "name", "smiles"]:
                    if col in df.columns:
                        id_col = col
                        break
                for _, row in df.iterrows():
                    if pd.notna(row.get(score_col)):
                        record = {
                            "molecule_id": row.get(id_col, "Unknown")
                            if id_col
                            else "Unknown",
                            "affinity": float(row[score_col]),
                            "model_name": row.get("model_name", "Unknown")
                            if "model_name" in df.columns
                            else model_lookup.get(str(row.get(id_col, "")), "Unknown"),
                        }
                        records.append(record)
        else:
            sdf_paths = [
                docking_dir / tool / f"{tool}_out.sdf",
                docking_dir / tool / "output.sdf",
                docking_dir / f"{tool}_out.sdf",
            ]
            for sdf_path in sdf_paths:
                if sdf_path.exists():
                    records = parse_docking_sdf(
                        base_path,
                        config,
                        stages,
                        initial_count,
                        final_count,
                        output_dir,
                        sdf_path,
                        model_lookup,
                    )
                    break
        if not records:
            continue
        result[tool]["raw_data"] = records
        scores = [r["affinity"] for r in records]
        result[tool]["summary"] = {
            "avg_affinity": sum(scores) / len(scores),
            "best_affinity": min(scores),
            "count": len(records),
        }
        sorted_records = sorted(records, key=lambda x: x["affinity"])
        result[tool]["top_molecules"] = sorted_records[:10]
        models = set(r["model_name"] for r in records if r["model_name"] != "Unknown")
        for model in models:
            model_records = [r for r in records if r["model_name"] == model]
            if not model_records:
                continue
            model_scores = [r["affinity"] for r in model_records]
            model_sorted = sorted(model_records, key=lambda x: x["affinity"])
            result[tool]["by_model"][model] = {
                "scores": model_scores,
                "summary": {
                    "avg_affinity": sum(model_scores) / len(model_scores),
                    "best_affinity": min(model_scores),
                    "count": len(model_records),
                },
                "top_molecules": model_sorted[:10],
            }
    return result


def get_docking_filters_detailed(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
) -> dict[str, Any]:
    """Get detailed docking filters data for report visualization.

    Reads metrics.csv and filtered_molecules.csv from the docking filters
    stage directory. Dynamically detects enabled filters from pass_* columns.

    Returns:
        Dictionary with total/passed/failed counts, per-filter stats,
        numeric metric distributions, and per-model breakdown.
    """
    df_dir = base_path / STAGE_DIRS["docking_filters"]
    metrics_path = df_dir / "metrics.csv"
    if not metrics_path.exists():
        return {}
    try:
        df = pd.read_csv(metrics_path)
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
        return {}
    if df.empty:
        return {}
    pass_cols = [c for c in df.columns if c.startswith("pass_") and c != "pass"]
    total_poses = len(df)
    passed_poses = int(df["pass"].sum()) if "pass" in df.columns else 0
    pass_rate = round(100.0 * passed_poses / total_poses, 1) if total_poses > 0 else 0.0
    filtered_path = df_dir / "filtered_molecules.csv"
    unique_molecules_passed = 0
    if filtered_path.exists():
        try:
            fdf = pd.read_csv(filtered_path)
            unique_molecules_passed = len(fdf)
        except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
            pass
    aggregation_mode = "all"
    dock_filt_config = config.get("docking_filters", {})
    if isinstance(dock_filt_config, dict):
        agg = dock_filt_config.get("aggregation", {})
        if isinstance(agg, dict):
            aggregation_mode = agg.get("mode", "all")
    per_filter = {}
    for col in pass_cols:
        filter_name = col.replace("pass_", "")
        if col in df.columns:
            total = int(df[col].notna().sum())
            passed = int(df[col].sum())
            per_filter[filter_name] = {
                "passed": passed,
                "total": total,
                "pass_rate": round(100.0 * passed / total, 1) if total > 0 else 0.0,
            }
    metric_columns = [
        "clashes",
        "strain_energy",
        "min_conformer_rmsd",
        "shape_score",
        "n_hbonds",
        "frac_atoms_outside_box",
    ]
    numeric_metrics = {}
    for col in metric_columns:
        if col in df.columns:
            values = df[col].dropna().tolist()
            if values:
                numeric_metrics[col] = values
    by_model = {}
    if "model_name" in df.columns:
        for model in df["model_name"].dropna().unique():
            model_df = df[df["model_name"] == model]
            m_total = len(model_df)
            m_passed = int(model_df["pass"].sum()) if "pass" in model_df.columns else 0
            by_model[str(model)] = {
                "total": m_total,
                "passed": m_passed,
                "pass_rate": round(100.0 * m_passed / m_total, 1)
                if m_total > 0
                else 0.0,
            }
    thresholds = {}
    if isinstance(dock_filt_config, dict):
        pq = dock_filt_config.get("pose_quality", {})
        if isinstance(pq, dict):
            if pq.get("max_clashes") is not None:
                thresholds["clashes"] = {"max": pq["max_clashes"]}
            if pq.get("max_strain_energy") is not None:
                thresholds["strain_energy"] = {"max": pq["max_strain_energy"]}
        cd = dock_filt_config.get("conformer_deviation", {})
        if isinstance(cd, dict):
            if cd.get("max_rmsd_to_conformer") is not None:
                thresholds["min_conformer_rmsd"] = {"max": cd["max_rmsd_to_conformer"]}
        sb = dock_filt_config.get("search_box", {})
        if isinstance(sb, dict):
            if sb.get("max_outside_fraction") is not None:
                thresholds["frac_atoms_outside_box"] = {
                    "max": sb["max_outside_fraction"]
                }
        ss = dock_filt_config.get("shepherd_score", {})
        if isinstance(ss, dict):
            if ss.get("min_shape_score") is not None:
                thresholds["shape_score"] = {"min": ss["min_shape_score"]}
    return {
        "total_poses": total_poses,
        "passed_poses": passed_poses,
        "pass_rate": pass_rate,
        "unique_molecules_passed": unique_molecules_passed,
        "aggregation_mode": aggregation_mode,
        "per_filter": per_filter,
        "numeric_metrics": numeric_metrics,
        "thresholds": thresholds,
        "by_model": by_model,
    }


def get_moleval_metrics(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
) -> dict[str, Any]:
    """Compute MolEval generative metrics across pipeline stages.

    Loads moleval config from config_moleval path, collects SMILES from
    each stage, and computes intrinsic distribution quality metrics.

    Returns:
        Dictionary with 'by_stage', 'stages', 'metrics' keys, or empty dict.
    """
    moleval_config = load_moleval_config(
        base_path, config, stages, initial_count, final_count, output_dir
    )
    if not moleval_config or not moleval_config.get("run", True):
        return {}
    stage_smiles = collect_stage_smiles(
        base_path, config, stages, initial_count, final_count, output_dir
    )
    if not stage_smiles:
        return {}
    try:
        return moleval_metrics.compute_stage_metrics(
            stage_smiles, moleval_config, seed=moleval_config.get("seed", 42)
        )
    except Exception as e:  # noqa: BLE001 — intentional: MolEval metrics may raise diverse errors
        logger.debug("MolEval metrics computation failed: %s", e)
        return {}


def load_stage_config(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
    config_key: str,
) -> dict[str, Any]:
    """Load a stage YAML config by its key in the pipeline config.

    Args:
        config_key: Key in self.config that points to the YAML config file path.

    Returns:
        Parsed config dict, or empty dict if not configured or not found.
    """
    import yaml

    config_path = config.get(config_key, "")
    if not config_path:
        return {}
    config_path = Path(config_path)
    if not config_path.is_absolute():
        for base in [base_path, Path.cwd()]:
            candidate = base / config_path
            if candidate.exists():
                config_path = candidate
                break
    if not config_path.exists():
        logger.debug("Stage config not found: %s (key=%s)", config_path, config_key)
        return {}
    try:
        with open(config_path) as f:
            return yaml.safe_load(f) or {}
    except (OSError, yaml.YAMLError) as e:
        logger.debug("Failed to load stage config %s: %s", config_key, e)
        return {}


def load_moleval_config(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
) -> dict[str, Any]:
    """Load MolEval configuration from the path specified in pipeline config.

    Returns:
        MolEval config dict, or empty dict if not configured.
    """
    return load_stage_config(
        base_path,
        config,
        stages,
        initial_count,
        final_count,
        output_dir,
        "config_moleval",
    )


def collect_stage_smiles(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
) -> dict[str, list[str]]:
    """Collect SMILES from key pipeline stages for MolEval analysis.

    Reads SMILES from CSV files at each pipeline checkpoint.

    Returns:
        Mapping of stage display name to list of SMILES strings.
    """
    stage_paths = [
        ("Input", "input/sampled_molecules.csv"),
        ("MolPrep", "stages/00_mol_prep/filtered_molecules.csv"),
        (
            "Descriptors",
            "stages/01_descriptors_initial/filtered/filtered_molecules.csv",
        ),
        ("StructFilters", "stages/03_structural_filters_post/filtered_molecules.csv"),
        ("Synthesis", "stages/04_synthesis/filtered_molecules.csv"),
        ("DockingFilters", "stages/06_docking_filters/filtered_molecules.csv"),
    ]
    result: dict[str, list[str]] = {}
    for stage_name, rel_path in stage_paths:
        smiles = read_stage_smiles(
            base_path,
            config,
            stages,
            initial_count,
            final_count,
            output_dir,
            base_path / rel_path,
        )
        if smiles:
            result[stage_name] = smiles
    return result


def read_stage_smiles(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
    path: Path,
) -> list[str]:
    """Read SMILES column from a stage CSV file.

    Args:
        path: Path to CSV file.

    Returns:
        List of SMILES strings, or empty list if not available.
    """
    if not path.exists():
        return []
    try:
        df = pd.read_csv(path)
        smiles_col = None
        for col in df.columns:
            if col.lower() == "smiles":
                smiles_col = col
                break
        if smiles_col is None:
            return []
        return df[smiles_col].dropna().astype(str).tolist()
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError) as e:
        logger.debug("Could not read SMILES from %s: %s", path, e)
        return []


def get_config_summary(
    base_path: Path,
    config: dict[str, Any],
    stages: list[Any],
    initial_count: int,
    final_count: int,
    output_dir: Path,
) -> dict[str, Any]:
    """Get configuration summary."""
    folder = config.get("folder_to_save", "")
    if not folder:
        folder = str(base_path)
    if stages:
        stages_enabled = [s.name for s in stages if s.enabled]
    else:
        stages_enabled = []
        stages_dir = base_path / "stages"
        if stages_dir.exists():
            stage_names = {
                "00_mol_prep": "Mol Prep",
                "01_descriptors_initial": "Descriptors (Initial)",
                "03_structural_filters_post": "Structural Filters (Post)",
                "04_synthesis": "Synthesis Analysis",
                "05_docking": "Molecular Docking",
                "06_docking_filters": "Docking Filters",
                "07_descriptors_final": "Descriptors (Final)",
            }
            for dir_name, display_name in stage_names.items():
                if (stages_dir / dir_name).exists():
                    stages_enabled.append(display_name)
    return {"folder_to_save": folder, "stages_enabled": stages_enabled}
