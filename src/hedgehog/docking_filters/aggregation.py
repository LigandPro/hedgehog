"""Pass/fail aggregation, single-pose collapse, and result saving."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem

from hedgehog.configs.logger import logger


def _collapse_to_single_pose(
    mols: list[Chem.Mol],
    results_df: pd.DataFrame,
) -> tuple[list[Chem.Mol], pd.DataFrame]:
    """Keep one best pose per source molecule id.

    Priority: lowest gnina_minimizedAffinity, then lowest pose index.
    """
    if results_df.empty:
        return [], results_df.copy()

    collapsed = results_df.copy()
    source_ids = collapsed["source_mol_idx"].astype(str).str.strip()
    fallback_ids = collapsed["mol_idx"].astype(str)
    collapsed["source_mol_idx"] = source_ids.where(source_ids != "", fallback_ids)

    collapsed["_affinity_sort"] = pd.to_numeric(
        collapsed["gnina_minimizedAffinity"], errors="coerce"
    ).fillna(float("inf"))

    collapsed = collapsed.sort_values(
        by=["source_mol_idx", "_affinity_sort", "mol_idx"],
        ascending=[True, True, True],
        kind="mergesort",
    ).drop_duplicates(subset=["source_mol_idx"], keep="first")

    rows: list[dict[str, Any]] = []
    selected_mols: list[Chem.Mol] = []

    for _, row in collapsed.iterrows():
        pose_idx = int(row["mol_idx"])
        if 0 <= pose_idx < len(mols) and mols[pose_idx] is not None:
            row_dict = row.to_dict()
            row_dict.pop("_affinity_sort", None)
            rows.append(row_dict)
            selected_mols.append(mols[pose_idx])

    collapsed_df = pd.DataFrame(rows)
    if collapsed_df.empty:
        return [], collapsed_df

    collapsed_df = collapsed_df.reset_index(drop=True)
    collapsed_df["mol_idx"] = range(len(collapsed_df))
    return selected_mols, collapsed_df


def aggregate_pass_columns(
    results_df: pd.DataFrame,
    agg_mode: str,
) -> pd.DataFrame:
    """Aggregate pass columns into a single 'pass' column.

    Args:
        results_df: DataFrame with pass_* columns.
        agg_mode: "all" or "any".

    Returns:
        DataFrame with added 'pass' column.
    """
    pass_cols = [c for c in results_df.columns if c.startswith("pass_")]
    results_df[pass_cols] = results_df[pass_cols].fillna(False)

    if agg_mode == "all":
        results_df["pass"] = results_df[pass_cols].all(axis=1)
    else:  # "any"
        results_df["pass"] = results_df[pass_cols].any(axis=1)

    return results_df


def save_results(
    results_df: pd.DataFrame,
    mols: list[Chem.Mol],
    output_dir: Path,
    docking_dir: Path,
    filter_config: dict[str, Any],
    filters_applied: list[str],
) -> pd.DataFrame:
    """Save filtered results to disk.

    Args:
        results_df: DataFrame with all filter results and 'pass' column.
        mols: Full molecule list.
        output_dir: Output directory for this stage.
        docking_dir: Docking stage directory (for ligands.csv lookup).
        filter_config: Filter configuration dict.
        filters_applied: List of applied filter names.

    Returns:
        The results_df (unmodified).
    """
    n_passed = results_df["pass"].sum()
    n_total = len(results_df)
    logger.info("Docking filters complete: %d/%d molecules passed", n_passed, n_total)
    logger.info("Filters applied: %s", ", ".join(filters_applied))

    # Save metrics
    if filter_config.get("aggregation", {}).get("save_metrics", True):
        metrics_path = output_dir / "metrics.csv"
        results_df.to_csv(metrics_path, index=False)
        logger.info("Saved metrics to %s", metrics_path)

    # Save filtered molecules (always create the file to make pipeline state explicit)
    filtered_df = results_df[results_df["pass"]].copy()
    filtered_path = output_dir / "filtered_molecules.csv"

    if filtered_df.empty:
        pd.DataFrame(columns=["smiles", "model_name", "mol_idx"]).to_csv(
            filtered_path, index=False
        )
        logger.info("Saved 0 filtered molecules to %s", filtered_path)
    else:
        pose_indices = filtered_df["mol_idx"].tolist()

        # Use original SMILES from ligands.csv (preserves 2D stereochemistry)
        # instead of generating from 3D poses which can resolve stereo differently.
        ligands_path = docking_dir / "ligands.csv"
        smiles_lookup: dict[str, str] = {}
        if ligands_path.exists():
            lig_df = pd.read_csv(ligands_path)
            smiles_lookup = dict(zip(lig_df["mol_idx"].astype(str), lig_df["smiles"]))

        fallback_smiles = pd.Series(
            [
                Chem.MolToSmiles(mols[i]) if 0 <= i < len(mols) and mols[i] else ""
                for i in pose_indices
            ],
            index=filtered_df.index,
        )
        filtered_df["smiles"] = (
            filtered_df["source_mol_idx"]
            .astype(str)
            .map(smiles_lookup)
            .fillna(fallback_smiles)
        )

        # For downstream pipeline stages, mol_idx should refer to the original molecule id
        # (not the pose index inside the SDF).
        filtered_df["mol_idx"] = filtered_df["source_mol_idx"]
        filtered_df = filtered_df.drop(columns=["source_mol_idx"])

        # Save all passing poses to CSV (pose-level detail)
        all_poses_path = output_dir / "filtered_poses.csv"
        filtered_df.to_csv(all_poses_path, index=False)
        logger.info("Saved %d filtered poses to %s", len(filtered_df), all_poses_path)

        # Deduplicate to unique molecules for downstream stages.
        # Keep the best pose per molecule (lowest minimizedAffinity).
        aff_col = "gnina_minimizedAffinity"
        if aff_col in filtered_df.columns:
            filtered_df = filtered_df.sort_values(aff_col, ascending=True)
        dedup_df = filtered_df.drop_duplicates(subset=["mol_idx"], keep="first")
        dedup_df[["smiles", "model_name", "mol_idx"]].to_csv(filtered_path, index=False)
        logger.info(
            "Saved %d unique molecules to %s (from %d poses)",
            len(dedup_df),
            filtered_path,
            len(filtered_df),
        )

        # Save filtered SDF
        filtered_sdf_path = output_dir / "filtered_poses.sdf"
        writer = Chem.SDWriter(str(filtered_sdf_path))
        for i in pose_indices:
            if 0 <= i < len(mols) and mols[i]:
                writer.write(mols[i])
        writer.close()
        logger.info("Saved filtered poses to %s", filtered_sdf_path)

    # Save failed molecules if configured
    if filter_config.get("aggregation", {}).get("save_failed", False):
        failed_df = results_df[~results_df["pass"]]
        if not failed_df.empty:
            failed_path = output_dir / "failed_molecules.csv"
            failed_df.to_csv(failed_path, index=False)
            logger.info("Saved %d failed molecules to %s", len(failed_df), failed_path)

    return results_df
