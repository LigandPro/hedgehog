"""Main entry point for docking filters stage."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem

from hedgehog.configs.logger import load_config, logger

from .utils import (
    apply_conformer_deviation_filter,
    apply_interaction_filter,
    apply_pose_quality_filter,
    apply_search_box_filter,
    apply_shepherd_score_filter,
    load_molecules_from_sdf,
)


def _project_root() -> Path:
    # src/hedgehog/stages/dockingFilters/main.py -> project root
    return Path(__file__).resolve().parent.parent.parent.parent.parent


def _resolve_existing_path(base: Path, path: str | Path) -> Path:
    """Resolve a possibly-relative path to an existing absolute path.

    Important: repository configs typically use paths relative to project root,
    while runtime artifacts are relative to the results folder.
    """
    p = Path(path)
    if p.is_absolute():
        return p

    candidates = [
        (base / p).resolve(),
        (_project_root() / p).resolve(),
        (Path.cwd() / p).resolve(),
    ]
    for c in candidates:
        if c.exists():
            return c

    # Fall back to base-relative absolute path for better error messages
    return (base / p).resolve()


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
    except Exception:
        return None


def docking_filters_main(config: dict[str, Any]) -> pd.DataFrame | None:
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
    base_folder = Path(config["folder_to_save"]).resolve()

    # Get filter config (pipeline config contains file paths, not dicts)
    filter_cfg_path = config.get("config_docking_filters")
    if not filter_cfg_path:
        logger.error("Docking filters config path is missing (config_docking_filters)")
        return None

    filter_config = load_config(filter_cfg_path)
    if not filter_config.get("run", False):
        logger.info("Docking filters disabled in config")
        return None

    # Determine paths
    output_dir = base_folder / "stages" / "06_docking_filters"
    output_dir.mkdir(parents=True, exist_ok=True)

    docking_dir = base_folder / "stages" / "05_docking"

    # Find input SDF
    input_sdf = filter_config.get("input_sdf")
    if input_sdf is None:
        # Try to find docking output from known locations
        candidates = [
            docking_dir / "smina" / "smina_out.sdf",
            docking_dir / "gnina" / "gnina_out.sdf",
        ]
        input_sdf = next((p for p in candidates if p.exists()), None)
        if input_sdf is None:
            logger.error("No docking output found. Run docking stage first.")
            return None
    else:
        input_sdf = _resolve_existing_path(base_folder, input_sdf)

    logger.info(f"Input SDF: {input_sdf}")

    # Load docking config (optional; used as fallback for receptor path)
    docking_config: dict[str, Any] = {}
    docking_cfg_path = config.get("config_docking")
    if docking_cfg_path:
        try:
            docking_config = load_config(docking_cfg_path)
        except Exception as e:
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
            protein_pdb = _resolve_existing_path(base_folder, protein_pdb_raw)
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
        logger.info(f"Protein PDB: {protein_pdb}")

    # Load molecules
    try:
        mols = load_molecules_from_sdf(input_sdf)
    except Exception as e:
        logger.error(f"Failed to load molecules: {e}")
        return None

    if not mols:
        logger.warning("No molecules loaded from SDF")
        return None

    # Extract identifiers (best-effort) from SDF properties
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

    # Initialize results DataFrame
    results_df = pd.DataFrame(
        {
            "mol_idx": range(len(mols)),
            "model_name": model_names,
            "source_mol_idx": mol_idxs,
            "gnina_minimizedAffinity": gnina_min_aff,
            "gnina_CNNscore": gnina_cnn_score,
            "gnina_CNNaffinity": gnina_cnn_aff,
        }
    )

    # Apply filters
    filters_applied = []

    agg_mode = filter_config.get("aggregation", {}).get("mode", "all")

    # Filter 0: Search-box containment (fast)
    sb_config = filter_config.get("search_box", {})
    try:
        sb_df = apply_search_box_filter(mols, base_folder, docking_config, sb_config)
        results_df = results_df.merge(sb_df, on="mol_idx", how="left")
        filters_applied.append("search_box")
    except Exception as e:
        logger.error(f"Search-box filter failed: {e}")
        results_df["pass_search_box"] = True

    # Optional optimization: under aggregation mode "all", if a pose fails search-box
    # containment it cannot pass the overall filter, so we can skip heavier checks.
    sb_short_circuit = bool(sb_config.get("short_circuit", True)) and agg_mode == "all"
    active_pose_indices = list(range(len(mols)))
    if sb_short_circuit and "pass_search_box" in results_df.columns:
        active_pose_indices = results_df.loc[
            results_df["pass_search_box"] == True, "mol_idx"  # noqa: E712
        ].tolist()

    # Filter 1: Pose Quality
    if pq_config.get("enabled", True):
        if active_pose_indices:
            try:
                mols_active = [mols[i] for i in active_pose_indices]
                pq_df = apply_pose_quality_filter(mols_active, protein_pdb, pq_config)
                pq_df["mol_idx"] = [active_pose_indices[i] for i in pq_df["mol_idx"]]
                results_df = results_df.merge(pq_df, on="mol_idx", how="left")
                filters_applied.append("pose_quality")
            except Exception as e:
                logger.error(f"Pose quality filter failed: {e}")
                results_df["pass_pose_quality"] = True
        else:
            results_df["pass_pose_quality"] = False

    # Filter 2: Interactions
    if int_config.get("enabled", True):
        if active_pose_indices:
            try:
                mols_active = [mols[i] for i in active_pose_indices]
                int_df = apply_interaction_filter(mols_active, protein_pdb, int_config)
                int_df["mol_idx"] = [active_pose_indices[i] for i in int_df["mol_idx"]]
                results_df = results_df.merge(int_df, on="mol_idx", how="left")
                filters_applied.append("interactions")
            except Exception as e:
                logger.error(f"Interaction filter failed: {e}")
                results_df["pass_interactions"] = True
        else:
            results_df["pass_interactions"] = False

    # Filter 3: Shepherd-Score
    ss_config = filter_config.get("shepherd_score", {})
    if ss_config.get("enabled", False):
        ref_path = ss_config.get("reference_ligand")
        if ref_path and Path(ref_path).exists():
            try:
                ref_mol = Chem.MolFromMolFile(str(ref_path))
                if ref_mol:
                    if not active_pose_indices:
                        results_df["pass_shepherd_score"] = False
                    else:
                        mols_active = [mols[i] for i in active_pose_indices]
                        ss_df = apply_shepherd_score_filter(
                            mols_active, ref_mol, ss_config
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
            except Exception as e:
                logger.error(f"Shepherd-Score filter failed: {e}")
                results_df["pass_shepherd_score"] = True
        else:
            logger.info("Shepherd-Score disabled (no reference ligand)")
            results_df["pass_shepherd_score"] = True

    # Filter 4: Conformer Deviation
    cd_config = filter_config.get("conformer_deviation", {})
    if cd_config.get("enabled", True):
        if active_pose_indices:
            try:
                mols_active = [mols[i] for i in active_pose_indices]
                cd_df = apply_conformer_deviation_filter(mols_active, cd_config)
                cd_df["mol_idx"] = [active_pose_indices[i] for i in cd_df["mol_idx"]]
                results_df = results_df.merge(cd_df, on="mol_idx", how="left")
                filters_applied.append("conformer_deviation")
            except Exception as e:
                logger.error(f"Conformer deviation filter failed: {e}")
                results_df["pass_conformer_deviation"] = True
        else:
            results_df["pass_conformer_deviation"] = False

    # Aggregate pass columns
    pass_cols = [c for c in results_df.columns if c.startswith("pass_")]
    results_df[pass_cols] = results_df[pass_cols].fillna(False)

    if agg_mode == "all":
        results_df["pass"] = results_df[pass_cols].all(axis=1)
    else:  # "any"
        results_df["pass"] = results_df[pass_cols].any(axis=1)

    # Summary
    n_passed = results_df["pass"].sum()
    n_total = len(results_df)
    logger.info(f"Docking filters complete: {n_passed}/{n_total} molecules passed")
    logger.info(f"Filters applied: {', '.join(filters_applied)}")

    # Save results
    if filter_config.get("aggregation", {}).get("save_metrics", True):
        metrics_path = output_dir / "metrics.csv"
        results_df.to_csv(metrics_path, index=False)
        logger.info(f"Saved metrics to {metrics_path}")

    # Save filtered molecules (always create the file to make pipeline state explicit)
    filtered_df = results_df[results_df["pass"]].copy()
    filtered_path = output_dir / "filtered_molecules.csv"

    if filtered_df.empty:
        pd.DataFrame(columns=["smiles", "model_name", "mol_idx"]).to_csv(
            filtered_path, index=False
        )
        logger.info(f"Saved 0 filtered molecules to {filtered_path}")
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
        logger.info(f"Saved {len(filtered_df)} filtered poses to {all_poses_path}")

        # Deduplicate to unique molecules for downstream stages.
        # Keep the best pose per molecule (lowest minimizedAffinity).
        aff_col = "gnina_minimizedAffinity"
        if aff_col in filtered_df.columns:
            filtered_df = filtered_df.sort_values(aff_col, ascending=True)
        dedup_df = filtered_df.drop_duplicates(subset=["mol_idx"], keep="first")
        dedup_df[["smiles", "model_name", "mol_idx"]].to_csv(filtered_path, index=False)
        logger.info(
            f"Saved {len(dedup_df)} unique molecules to {filtered_path} "
            f"(from {len(filtered_df)} poses)"
        )

        # Save filtered SDF
        filtered_sdf_path = output_dir / "filtered_poses.sdf"
        writer = Chem.SDWriter(str(filtered_sdf_path))
        for i in pose_indices:
            if 0 <= i < len(mols) and mols[i]:
                writer.write(mols[i])
        writer.close()
        logger.info(f"Saved filtered poses to {filtered_sdf_path}")

    # Save failed molecules if configured
    if filter_config.get("aggregation", {}).get("save_failed", False):
        failed_df = results_df[~results_df["pass"]]
        if not failed_df.empty:
            failed_path = output_dir / "failed_molecules.csv"
            failed_df.to_csv(failed_path, index=False)
            logger.info(f"Saved {len(failed_df)} failed molecules to {failed_path}")

    return results_df
