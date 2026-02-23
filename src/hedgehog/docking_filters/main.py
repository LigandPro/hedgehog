"""Main entry point for docking filters stage."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem

from hedgehog.configs.logger import load_config, logger
from hedgehog.utils.parallel import resolve_n_jobs
from hedgehog.utils.paths import resolve_existing_path

from .utils import (
    _project_root,
    apply_conformer_deviation_filter,
    apply_interaction_filter,
    apply_pose_quality_filter,
    apply_posebusters_fast_filter,
    apply_search_box_filter,
    apply_shepherd_score_filter,
    apply_symmetry_rmsd_filter,
    load_molecules_from_sdf,
)

apply_posecheck_fast_filter = apply_posebusters_fast_filter


def _get_first_prop_value(mol: Chem.Mol, canonical_names: set[str]) -> str | None:
    """Return the first SDF property value whose (normalized) key matches any canonical name."""
    # Fast path: exact keys
    for name in canonical_names:
        if mol.HasProp(name):
            return mol.GetProp(name)

    # Some toolchains escape underscores in SDF property names (e.g. "s_sm_model\\_name").
    for prop in mol.GetPropNames():
        normalized = prop.replace("\\", "")
        if normalized in canonical_names:
            return mol.GetProp(prop)
    return None


def _get_prop_as_float(mol: Chem.Mol, prop_name: str) -> float | None:
    """Parse an SDF property as float (best-effort)."""
    if not mol.HasProp(prop_name):
        return None
    try:
        return float(mol.GetProp(prop_name))
    except (ValueError, TypeError):
        return None


def _collapse_to_single_pose(
    mols: list[Chem.Mol],
    results_df: pd.DataFrame,
) -> tuple[list[Chem.Mol], pd.DataFrame]:
    """Keep one best pose per source molecule id.

    Priority: lowest gnina_minimizedAffinity, then lowest pose index.
    """
    if results_df.empty:
        return [], results_df.copy()

    required_cols = {"source_mol_idx", "mol_idx", "gnina_minimizedAffinity"}
    missing = required_cols - set(results_df.columns)
    if missing:
        logger.warning(
            "Cannot collapse to single pose: missing columns %s", sorted(missing)
        )
        return list(mols), results_df.copy()

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


def _load_filter_inputs(
    config: dict[str, Any],
) -> dict[str, Any] | None:
    """Load and validate all inputs for docking filters.

    Returns a dict with keys: filter_config, output_dir, docking_dir,
    docking_config, protein_pdb, mols, n_jobs, base_folder.
    Returns None on error.
    """
    base_folder = Path(config["folder_to_save"]).resolve()

    filter_cfg_path = config.get("config_docking_filters")
    if not filter_cfg_path:
        logger.error("Docking filters config path is missing (config_docking_filters)")
        return None

    filter_config = load_config(filter_cfg_path)
    if not filter_config.get("run", False):
        logger.info("Docking filters disabled in config")
        return None

    n_jobs = resolve_n_jobs(filter_config, config)
    for sub_key in (
        "search_box",
        "pose_quality",
        "interactions",
        "shepherd_score",
        "conformer_deviation",
    ):
        sub = filter_config.get(sub_key)
        if isinstance(sub, dict) and "n_jobs" not in sub:
            sub["n_jobs"] = n_jobs

    output_dir = base_folder / "stages" / "06_docking_filters"
    output_dir.mkdir(parents=True, exist_ok=True)

    docking_dir = base_folder / "stages" / "05_docking"

    # Find input SDF
    input_sdf = filter_config.get("input_sdf")
    if input_sdf is None:
        candidates = [
            docking_dir / "smina" / "smina_out.sdf",
            docking_dir / "gnina" / "gnina_out.sdf",
        ]
        input_sdf = next((p for p in candidates if p.exists()), None)
        if input_sdf is None:
            logger.error("No docking output found. Run docking stage first.")
            return None
    else:
        input_sdf = resolve_existing_path(base_folder, input_sdf, _project_root())

    logger.info("Input SDF: %s", input_sdf)

    # Load docking config (optional; used as fallback for receptor path)
    docking_config: dict[str, Any] = {}
    docking_cfg_path = config.get("config_docking")
    if docking_cfg_path:
        try:
            docking_config = load_config(docking_cfg_path)
        except (OSError, ValueError) as e:
            logger.warning(
                "Failed to load docking config (%s): %s", docking_cfg_path, e
            )

    # Find protein PDB only if at least one enabled filter needs it
    pq_config = filter_config.get("pose_quality", {})
    int_config = filter_config.get("interactions", {})
    needs_protein = bool(
        pq_config.get("enabled", True) or int_config.get("enabled", True)
    )

    protein_pdb: Path | None = None
    if needs_protein:
        protein_pdb_raw = filter_config.get("receptor_pdb") or docking_config.get(
            "receptor_pdb"
        )
        if protein_pdb_raw:
            protein_pdb = resolve_existing_path(
                base_folder, protein_pdb_raw, _project_root()
            )
        else:
            prepared_pdb = docking_dir / "_workdir" / "protein_prepared.pdb"
            if not prepared_pdb.exists():
                # Legacy fallback
                prepared_pdb = docking_dir / "protein_prepared.pdb"
            if prepared_pdb.exists():
                protein_pdb = prepared_pdb

        if protein_pdb is None or not protein_pdb.exists():
            logger.error("No protein PDB found (required by enabled filters)")
            return None
        logger.info("Protein PDB: %s", protein_pdb)

    # Load molecules
    try:
        mols = load_molecules_from_sdf(input_sdf)
    except (OSError, ValueError, RuntimeError) as e:
        logger.error("Failed to load molecules: %s", e)
        return None

    if not mols:
        logger.warning("No molecules loaded from SDF")
        return None

    return {
        "filter_config": filter_config,
        "output_dir": output_dir,
        "docking_dir": docking_dir,
        "docking_config": docking_config,
        "protein_pdb": protein_pdb,
        "mols": mols,
        "n_jobs": n_jobs,
        "base_folder": base_folder,
    }


def _extract_pose_identifiers(mols: list[Chem.Mol]) -> pd.DataFrame:
    """Extract identifiers and affinity data from SDF properties.

    Returns DataFrame with columns: mol_idx, model_name, source_mol_idx,
    gnina_minimizedAffinity, gnina_CNNscore, gnina_CNNaffinity.
    """
    model_names: list[str] = []
    mol_idxs: list[str] = []
    gnina_min_aff: list[float | None] = []
    gnina_cnn_score: list[float | None] = []
    gnina_cnn_aff: list[float | None] = []

    for mol in mols:
        model_name = _get_first_prop_value(
            mol, {"model_name", "sm_model_name", "s_sm_model_name"}
        )
        mol_idx = _get_first_prop_value(
            mol, {"mol_idx", "sm_mol_idx", "s_sm_mol_idx", "name"}
        )
        model_names.append(model_name or "")
        mol_idxs.append(mol_idx or "")

        gnina_min_aff.append(_get_prop_as_float(mol, "minimizedAffinity"))
        gnina_cnn_score.append(_get_prop_as_float(mol, "CNNscore"))
        gnina_cnn_aff.append(_get_prop_as_float(mol, "CNNaffinity"))

    return pd.DataFrame(
        {
            "mol_idx": range(len(mols)),
            "model_name": model_names,
            "source_mol_idx": mol_idxs,
            "gnina_minimizedAffinity": gnina_min_aff,
            "gnina_CNNscore": gnina_cnn_score,
            "gnina_CNNaffinity": gnina_cnn_aff,
        }
    )


def _run_single_filter(
    filter_name: str,
    filter_fn,
    mols: list[Chem.Mol],
    active_pose_indices: list[int],
    results_df: pd.DataFrame,
) -> tuple[pd.DataFrame, bool]:
    """Run one filter on active poses with standard error handling.

    Args:
        filter_name: Name for the pass_<name> column.
        filter_fn: Callable(mols_active) -> DataFrame with mol_idx column.
        mols: All molecules.
        active_pose_indices: Indices of molecules to filter.
        results_df: Existing results to merge into.

    Returns:
        (updated_results_df, filter_was_applied)
    """
    if not active_pose_indices:
        results_df[f"pass_{filter_name}"] = False
        return results_df, False

    try:
        mols_active = [mols[i] for i in active_pose_indices]
        filter_df = filter_fn(mols_active)
        filter_df["mol_idx"] = [active_pose_indices[i] for i in filter_df["mol_idx"]]
        results_df = results_df.merge(filter_df, on="mol_idx", how="left")
        return results_df, True
    except Exception as e:  # noqa: BLE001 — intentional: filter failure should not crash pipeline
        logger.warning(
            "Filter '%s' failed, defaulting all poses to pass: %s", filter_name, e
        )
        results_df[f"pass_{filter_name}"] = True
        return results_df, False


def _merge_filter_results(
    results_df: pd.DataFrame,
    agg_mode: str,
    filters_applied: list[str],
) -> pd.DataFrame:
    """Aggregate per-filter pass columns into final pass/fail and log summary."""
    pass_cols = [c for c in results_df.columns if c.startswith("pass_")]
    results_df[pass_cols] = results_df[pass_cols].fillna(False)

    if agg_mode == "all":
        results_df["pass"] = results_df[pass_cols].all(axis=1)
    else:  # "any"
        results_df["pass"] = results_df[pass_cols].any(axis=1)

    n_passed = results_df["pass"].sum()
    n_total = len(results_df)
    logger.info("Docking filters complete: %d/%d molecules passed", n_passed, n_total)
    logger.info("Filters applied: %s", ", ".join(filters_applied))

    return results_df


def _save_filter_outputs(
    results_df: pd.DataFrame,
    mols: list[Chem.Mol],
    output_dir: Path,
    docking_dir: Path,
    filter_config: dict[str, Any],
) -> None:
    """Save metrics, filtered molecules, SDF, and optionally failed molecules."""
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
        # mol_idx is always present (set by pipeline); smiles is added below
        if "mol_idx" not in filtered_df.columns:
            logger.warning("Cannot save filter outputs: missing 'mol_idx' column")
            return
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


def docking_filters_main(config: dict[str, Any], reporter=None) -> pd.DataFrame | None:
    """
    Main entry point for docking filters stage.

    Applies configured filters to docking poses and saves filtered results.

    Args:
        config: Full pipeline configuration dict containing:
            - config_docking_filters: Filter settings
            - folder_to_save: Output directory path
            - config_docking: Docking configuration (for receptor path)

    Returns:
        DataFrame with filtered molecules and metrics, or None if no molecules pass
    """
    # Load and validate inputs
    inputs = _load_filter_inputs(config)
    if inputs is None:
        return None

    filter_config = inputs["filter_config"]
    output_dir = inputs["output_dir"]
    docking_dir = inputs["docking_dir"]
    docking_config = inputs["docking_config"]
    protein_pdb = inputs["protein_pdb"]
    mols = inputs["mols"]
    base_folder = inputs["base_folder"]

    # Extract identifiers and collapse to single pose per molecule
    results_df = _extract_pose_identifiers(mols)

    initial_pose_count = len(mols)
    mols, results_df = _collapse_to_single_pose(mols, results_df)
    if not mols:
        logger.warning("No valid poses remained after single-pose collapse")
        return None
    logger.info(
        "Collapsed docking poses to one pose per molecule: %d -> %d",
        initial_pose_count,
        len(mols),
    )

    # Resolve filter configs
    filters_applied: list[str] = []
    pq_config = filter_config.get("pose_quality", {})
    int_config = filter_config.get("interactions", {})
    sb_config = filter_config.get("search_box", {})
    ss_config = filter_config.get("shepherd_score", {})
    cd_config = filter_config.get("conformer_deviation", {})
    agg_mode = filter_config.get("aggregation", {}).get("mode", "all")

    ref_path_raw = ss_config.get("reference_ligand")
    ref_path = None
    if ref_path_raw:
        ref_path = resolve_existing_path(base_folder, ref_path_raw, _project_root())
    shepherd_enabled = bool(ss_config.get("enabled", False)) and bool(
        ref_path and ref_path.exists()
    )

    # Stage progress setup
    step_names: list[str] = ["search_box"]
    if pq_config.get("enabled", True):
        step_names.append("pose_quality")
    if int_config.get("enabled", True):
        step_names.append("interactions")
    if shepherd_enabled:
        step_names.append("shepherd_score")
    if cd_config.get("enabled", True):
        step_names.append("conformer_deviation")

    stage_total = max(1, len(step_names) * 100)

    def _step_progress(step_index: int, label: str):
        if reporter is None:
            return None

        base = step_index * 100

        def _progress(done: int, total: int) -> None:
            if total <= 0:
                pct = 0
            else:
                pct = int(round((done / total) * 100))
            pct = max(0, min(100, pct))
            reporter.progress(
                base + pct, stage_total, message=f"DockingFilters: {label}"
            )

        return _progress

    def _report(value: int, msg: str) -> None:
        if reporter is not None:
            reporter.progress(value, stage_total, message=msg)

    # Filter 0: Search-box containment (fast) — runs on all mols, no active_pose_indices
    try:
        _report(0, "DockingFilters: search_box")
        sb_df = apply_search_box_filter(mols, base_folder, docking_config, sb_config)
        results_df = results_df.merge(sb_df, on="mol_idx", how="left")
        filters_applied.append("search_box")
    except Exception as e:  # noqa: BLE001 — intentional: filter failure should not crash pipeline
        logger.error("Search-box filter failed: %s", e)
        results_df["pass_search_box"] = True
    finally:
        _report(100, "DockingFilters: search_box")

    # Short-circuit: under aggregation mode "all", skip heavier checks for failed poses
    sb_short_circuit = bool(sb_config.get("short_circuit", True)) and agg_mode == "all"
    active_pose_indices = list(range(len(mols)))
    if sb_short_circuit and "pass_search_box" in results_df.columns:
        active_pose_indices = results_df.loc[
            results_df["pass_search_box"] == True, "mol_idx"  # noqa: E712
        ].tolist()

    # Filter 1: Pose Quality
    if pq_config.get("enabled", True):
        pq_step_idx = step_names.index("pose_quality")
        _report(pq_step_idx * 100, "DockingFilters: pose_quality")

        pq_backend = pq_config.get("backend", "posebusters_fast")

        def _pose_quality_fn(mols_active):
            if pq_backend == "posebusters_fast":
                return apply_posebusters_fast_filter(
                    mols_active,
                    protein_pdb,
                    pq_config,
                    progress_cb=_step_progress(pq_step_idx, "pose_quality"),
                )
            else:  # "posecheck" (legacy)
                return apply_pose_quality_filter(mols_active, protein_pdb, pq_config)

        results_df, applied = _run_single_filter(
            "pose_quality", _pose_quality_fn, mols, active_pose_indices, results_df
        )
        if applied:
            filters_applied.append("pose_quality")
        _report((pq_step_idx + 1) * 100, "DockingFilters: pose_quality")

    # Short-circuit after pose quality
    pq_short_circuit = bool(pq_config.get("short_circuit", True)) and agg_mode == "all"
    if pq_short_circuit and "pass_pose_quality" in results_df.columns:
        pq_mask = results_df["pass_pose_quality"].fillna(False) == True  # noqa: E712
        if "pass_search_box" in results_df.columns:
            sb_mask = results_df["pass_search_box"].fillna(False) == True  # noqa: E712
            pq_mask = pq_mask & sb_mask
        active_pose_indices = results_df.loc[pq_mask, "mol_idx"].tolist()

    # Filter 2: Interactions
    if int_config.get("enabled", True):
        int_step_idx = step_names.index("interactions")
        _report(int_step_idx * 100, "DockingFilters: interactions")

        def _interactions_fn(mols_active):
            return apply_interaction_filter(mols_active, protein_pdb, int_config)

        results_df, applied = _run_single_filter(
            "interactions", _interactions_fn, mols, active_pose_indices, results_df
        )
        if applied:
            filters_applied.append("interactions")
        _report((int_step_idx + 1) * 100, "DockingFilters: interactions")

    # Filter 3: Shepherd-Score
    if ss_config.get("enabled", False):
        if ref_path and ref_path.exists():
            ss_step_idx = (
                step_names.index("shepherd_score")
                if "shepherd_score" in step_names
                else None
            )
            if reporter is not None and ss_step_idx is not None:
                reporter.progress(
                    ss_step_idx * 100,
                    stage_total,
                    message="DockingFilters: shepherd_score",
                )
            try:
                ref_mol = Chem.MolFromMolFile(str(ref_path))
                if ref_mol:
                    if not active_pose_indices:
                        results_df["pass_shepherd_score"] = False
                    else:
                        mols_active = [mols[i] for i in active_pose_indices]
                        ss_df = apply_shepherd_score_filter(
                            mols_active,
                            ref_mol,
                            ss_config,
                            progress_cb=_step_progress(
                                ss_step_idx or 0, "shepherd_score"
                            )
                            if ss_step_idx is not None
                            else None,
                        )
                        ss_df["mol_idx"] = [
                            active_pose_indices[i] for i in ss_df["mol_idx"]
                        ]
                        results_df = results_df.merge(ss_df, on="mol_idx", how="left")
                        filters_applied.append("shepherd_score")
                else:
                    logger.warning(
                        "Failed to load reference molecule for Shepherd-Score"
                    )
                    results_df["pass_shepherd_score"] = True
            except Exception as e:  # noqa: BLE001 — intentional: filter failure should not crash pipeline
                logger.error("Shepherd-Score filter failed: %s", e)
                results_df["pass_shepherd_score"] = True
            finally:
                if reporter is not None and ss_step_idx is not None:
                    reporter.progress(
                        (ss_step_idx + 1) * 100,
                        stage_total,
                        message="DockingFilters: shepherd_score",
                    )
        else:
            logger.info("Shepherd-Score disabled (no reference ligand)")
            results_df["pass_shepherd_score"] = True

    # Filter 4: Conformer Deviation
    if cd_config.get("enabled", True):
        cd_step_idx = step_names.index("conformer_deviation")
        _report(cd_step_idx * 100, "DockingFilters: conformer_deviation")

        cd_backend = cd_config.get("backend", "symmetry_rmsd")

        def _conformer_deviation_fn(mols_active):
            if cd_backend == "symmetry_rmsd":
                return apply_symmetry_rmsd_filter(
                    mols_active,
                    cd_config,
                    progress_cb=_step_progress(cd_step_idx, "conformer_deviation"),
                )
            else:  # "naive" (legacy)
                return apply_conformer_deviation_filter(
                    mols_active,
                    cd_config,
                    progress_cb=_step_progress(cd_step_idx, "conformer_deviation"),
                )

        results_df, applied = _run_single_filter(
            "conformer_deviation",
            _conformer_deviation_fn,
            mols,
            active_pose_indices,
            results_df,
        )
        if applied:
            filters_applied.append("conformer_deviation")
        _report((cd_step_idx + 1) * 100, "DockingFilters: conformer_deviation")

    _report(stage_total, "DockingFilters complete")

    # Aggregate and save
    results_df = _merge_filter_results(results_df, agg_mode, filters_applied)
    _save_filter_outputs(results_df, mols, output_dir, docking_dir, filter_config)

    return results_df
