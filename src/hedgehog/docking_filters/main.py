"""Main entry point for docking filters stage."""

from __future__ import annotations

import ast
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem

from hedgehog._constants import CFG_DOCKING, KEY_FOLDER_TO_SAVE
from hedgehog.configs.logger import load_config, logger
from hedgehog.utils.parallel import resolve_n_jobs

from .utils import (
    _resolve_conformer_backend,
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

_INTERACTION_REPORT_FILES = (
    "interaction_events.csv",
    "interaction_residue_summary.csv",
    "interaction_type_summary.csv",
    "interaction_matrix.csv",
    "interaction_report_meta.json",
)


def _project_root() -> Path:
    # src/hedgehog/docking_filters/main.py -> project root
    return Path(__file__).resolve().parent.parent.parent.parent


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


def _run_single_filter(
    *,
    filter_name: str,
    reporter,
    step_progress_fn,
    mols: list[Chem.Mol],
    active_pose_indices: list[int],
    results_df: pd.DataFrame,
    filters_applied: list[str],
    pass_col: str,
    run_fn,
) -> pd.DataFrame:
    """Run a single filter with standard progress/error/short-circuit handling.

    This eliminates the repeated boilerplate per filter in ``docking_filters_main``.

    Args:
        filter_name: Human label for log messages (e.g. ``"pose_quality"``).
        reporter: Optional progress reporter.
        step_progress_fn: Callable ``(label, step_total) -> progress_cb | None``.
        mols: Full molecule list.
        active_pose_indices: Indices of molecules still active after short-circuit.
        results_df: Cumulative results DataFrame (modified in place).
        filters_applied: Accumulator list of successfully applied filter names.
        pass_col: Name of the boolean pass column (e.g. ``"pass_pose_quality"``).
        run_fn: Callable ``(mols_active, progress_cb) -> pd.DataFrame`` that runs the
            actual filter logic.  Must return a DataFrame with ``mol_idx`` and
            ``pass_col``.

    Returns:
        Updated ``results_df``.
    """
    step_total = len(active_pose_indices) if active_pose_indices else len(mols)
    if reporter is not None:
        reporter.progress(
            0,
            step_total,
            message=f"DockingFilters: {filter_name}",
        )
    if active_pose_indices:
        try:
            mols_active = [mols[i] for i in active_pose_indices]
            progress_cb = step_progress_fn(filter_name, step_total)
            df = run_fn(mols_active, progress_cb)
            df["mol_idx"] = [active_pose_indices[i] for i in df["mol_idx"]]
            results_df = results_df.merge(df, on="mol_idx", how="left")
            filters_applied.append(filter_name)
        except Exception as e:
            logger.error("%s filter failed: %s", filter_name, e)
            results_df[pass_col] = True
    else:
        results_df[pass_col] = False
    if reporter is not None:
        reporter.progress(
            step_total,
            step_total,
            message=f"DockingFilters: {filter_name}",
        )
    return results_df


def _cleanup_interaction_reporting_artifacts(output_dir: Path) -> None:
    for filename in _INTERACTION_REPORT_FILES:
        path = output_dir / filename
        if path.exists():
            try:
                path.unlink()
            except OSError as e:
                logger.warning("Could not remove %s: %s", path, e)


def _coerce_bool(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "y"}
    try:
        if pd.isna(value):
            return False
    except Exception as exc:
        logger.debug(
            "pd.isna(%r) raised %r; falling back to generic bool coercion",
            value,
            exc,
        )
    if isinstance(value, (int, float)):
        return value != 0
    return bool(value)


def _normalize_interaction_label(item: Any) -> dict[str, str] | None:
    if not isinstance(item, dict):
        return None

    label = str(item.get("label", "")).strip()
    if not label:
        return None
    residue = str(item.get("residue", "") or "unknown").strip() or "unknown"
    interaction_type = (
        str(item.get("interaction_type", "") or "unknown").strip() or "unknown"
    )
    return {
        "label": label,
        "residue": residue,
        "interaction_type": interaction_type,
    }


def _interaction_labels_from_row(row: pd.Series) -> list[dict[str, str]]:
    labels: list[dict[str, str]] = []
    payload = row.get("interaction_labels_json")
    if isinstance(payload, str) and payload.strip():
        try:
            parsed = json.loads(payload)
            if isinstance(parsed, list):
                for item in parsed:
                    normalized = _normalize_interaction_label(item)
                    if normalized is not None:
                        labels.append(normalized)
        except (TypeError, ValueError, json.JSONDecodeError):
            labels = []

    # Backward-compatible fallback for historical metrics.csv where only the
    # legacy "interactions" CSV string is available.
    if not labels:
        interactions_raw = row.get("interactions")
        if isinstance(interactions_raw, str) and interactions_raw.strip():
            labels.extend(_parse_legacy_interactions(interactions_raw))

    dedup: dict[tuple[str, str, str], dict[str, str]] = {}
    for label in labels:
        key = (label["label"], label["residue"], label["interaction_type"])
        dedup[key] = label
    return sorted(
        dedup.values(),
        key=lambda item: (item["residue"], item["interaction_type"], item["label"]),
    )


def _extract_interaction_pairs_from_literal(parsed: Any) -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    if not isinstance(parsed, (list, tuple)):
        return pairs

    # Single pair form: ("ASP:123", "HBDonor")
    if len(parsed) >= 2 and not any(
        isinstance(item, (list, tuple, dict)) for item in parsed[:2]
    ):
        residue = str(parsed[0]).strip()
        interaction_type = str(parsed[1]).strip()
        if residue or interaction_type:
            pairs.append((residue, interaction_type))
        return pairs

    # Sequence form: [("ASP:123", "HBDonor"), ("SER:45", "Hydrophobic")]
    for item in parsed:
        if not isinstance(item, (list, tuple)) or len(item) < 2:
            continue
        residue = str(item[0]).strip()
        interaction_type = str(item[1]).strip()
        if residue or interaction_type:
            pairs.append((residue, interaction_type))
    return pairs


def _parse_legacy_interactions(interactions_raw: str) -> list[dict[str, str]]:
    text = interactions_raw.strip()
    if not text:
        return []

    pair_tokens: list[tuple[str, str]] = []
    for candidate in (text, f"[{text}]"):
        try:
            parsed = ast.literal_eval(candidate)
        except (ValueError, SyntaxError):
            continue
        pair_tokens = _extract_interaction_pairs_from_literal(parsed)
        if pair_tokens:
            break

    if not pair_tokens:
        tuple_matches = re.findall(
            r"\(\s*['\"]?([^,'\")]+)['\"]?\s*,\s*['\"]?([^,'\")]+)['\"]?\s*\)",
            text,
        )
        pair_tokens = [
            (residue, interaction_type) for residue, interaction_type in tuple_matches
        ]

    rows: list[dict[str, str]] = []
    if pair_tokens:
        for residue, interaction_type in pair_tokens:
            residue_clean = residue.strip().strip("'\"") or "unknown"
            type_clean = interaction_type.strip().strip("'\"") or "unknown"
            rows.append(
                {
                    "label": f"{residue_clean}|{type_clean}",
                    "residue": residue_clean,
                    "interaction_type": type_clean,
                }
            )
        return rows

    for token in text.split(","):
        label = token.strip().strip("'\"")
        if not label:
            continue
        if "|" in label:
            residue, interaction_type = label.split("|", 1)
            residue_clean = residue.strip().strip("'\"") or "unknown"
            type_clean = interaction_type.strip().strip("'\"") or "unknown"
            rows.append(
                {
                    "label": label,
                    "residue": residue_clean,
                    "interaction_type": type_clean,
                }
            )
        else:
            rows.append(
                {
                    "label": label,
                    "residue": "unknown",
                    "interaction_type": "unknown",
                }
            )
    return rows


def _is_interaction_reporting_enabled(int_config: dict[str, Any]) -> bool:
    interaction_reporting_cfg = int_config.get("reporting", {})
    if isinstance(interaction_reporting_cfg, dict):
        return _coerce_bool(interaction_reporting_cfg.get("enabled", True))
    return _coerce_bool(interaction_reporting_cfg)


def _write_interaction_reporting_artifacts(
    *,
    results_df: pd.DataFrame,
    output_dir: Path,
    interaction_filter_enabled: bool,
) -> None:
    events_rows: list[dict[str, Any]] = []

    sorted_df = results_df.sort_values("mol_idx", kind="mergesort")
    for _, row in sorted_df.iterrows():
        labels = _interaction_labels_from_row(row)
        if not labels:
            continue

        raw_mol_idx = row.get("mol_idx")
        mol_idx = int(raw_mol_idx) if pd.notna(raw_mol_idx) else -1
        source_mol_idx = row.get("source_mol_idx")
        model_name = row.get("model_name")
        pass_all = _coerce_bool(row.get("pass", False))
        pass_interactions = _coerce_bool(row.get("pass_interactions", False))

        for label in labels:
            events_rows.append(
                {
                    "mol_idx": mol_idx,
                    "source_mol_idx": ""
                    if pd.isna(source_mol_idx)
                    else str(source_mol_idx),
                    "model_name": "" if pd.isna(model_name) else str(model_name),
                    "pass": pass_all,
                    "pass_interactions": pass_interactions,
                    "residue": label["residue"],
                    "interaction_type": label["interaction_type"],
                    "label": label["label"],
                }
            )

    events_df = pd.DataFrame(
        events_rows,
        columns=[
            "mol_idx",
            "source_mol_idx",
            "model_name",
            "pass",
            "pass_interactions",
            "residue",
            "interaction_type",
            "label",
        ],
    )
    events_df.to_csv(output_dir / "interaction_events.csv", index=False)

    residue_summary_rows: list[dict[str, Any]] = []
    type_summary_rows: list[dict[str, Any]] = []
    matrix_rows: list[dict[str, Any]] = []

    if not events_df.empty:
        for residue, group in events_df.groupby("residue", sort=False):
            residue_summary_rows.append(
                {
                    "residue": residue,
                    "total_events": int(len(group)),
                    "unique_poses": int(group["mol_idx"].nunique()),
                    "unique_passed_poses": int(
                        group.loc[group["pass"].fillna(False), "mol_idx"].nunique()
                    ),
                    "interaction_types": int(group["interaction_type"].nunique()),
                }
            )

        for interaction_type, group in events_df.groupby(
            "interaction_type", sort=False
        ):
            type_summary_rows.append(
                {
                    "interaction_type": interaction_type,
                    "total_events": int(len(group)),
                    "unique_poses": int(group["mol_idx"].nunique()),
                    "unique_residues": int(group["residue"].nunique()),
                }
            )

        for (residue, interaction_type), group in events_df.groupby(
            ["residue", "interaction_type"], sort=False
        ):
            matrix_rows.append(
                {
                    "residue": residue,
                    "interaction_type": interaction_type,
                    "event_count": int(len(group)),
                    "unique_poses": int(group["mol_idx"].nunique()),
                }
            )

    residue_summary_df = pd.DataFrame(
        residue_summary_rows,
        columns=[
            "residue",
            "total_events",
            "unique_poses",
            "unique_passed_poses",
            "interaction_types",
        ],
    )
    residue_summary_df = residue_summary_df.sort_values(
        ["total_events", "residue"], ascending=[False, True], kind="mergesort"
    )
    residue_summary_df.to_csv(
        output_dir / "interaction_residue_summary.csv", index=False
    )

    type_summary_df = pd.DataFrame(
        type_summary_rows,
        columns=[
            "interaction_type",
            "total_events",
            "unique_poses",
            "unique_residues",
        ],
    )
    type_summary_df = type_summary_df.sort_values(
        ["total_events", "interaction_type"], ascending=[False, True], kind="mergesort"
    )
    type_summary_df.to_csv(output_dir / "interaction_type_summary.csv", index=False)

    matrix_df = pd.DataFrame(
        matrix_rows,
        columns=["residue", "interaction_type", "event_count", "unique_poses"],
    )
    matrix_df = matrix_df.sort_values(
        ["residue", "interaction_type"], ascending=[True, True], kind="mergesort"
    )
    matrix_df.to_csv(output_dir / "interaction_matrix.csv", index=False)

    meta = {
        "schema_version": 1,
        "reporting_enabled": True,
        "interaction_filter_enabled": interaction_filter_enabled,
        "total_poses": int(len(results_df)),
        "poses_with_interactions": int(events_df["mol_idx"].nunique())
        if not events_df.empty
        else 0,
        "total_events": int(len(events_df)),
        "unique_residues": int(residue_summary_df["residue"].nunique())
        if not residue_summary_df.empty
        else 0,
        "unique_interaction_types": int(type_summary_df["interaction_type"].nunique())
        if not type_summary_df.empty
        else 0,
    }
    with open(output_dir / "interaction_report_meta.json", "w") as f:
        json.dump(meta, f, indent=2, sort_keys=True)


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
    base_folder = Path(config[KEY_FOLDER_TO_SAVE]).resolve()

    # Get filter config (pipeline config contains file paths, not dicts)
    filter_cfg_path = config.get("config_docking_filters")
    if not filter_cfg_path:
        logger.error("Docking filters config path is missing (config_docking_filters)")
        return None

    filter_config = load_config(filter_cfg_path)
    if not filter_config.get("run", False):
        logger.info("Docking filters disabled in config")
        return None

    # Resolve n_jobs once and propagate to sub-configs that don't override it
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

    logger.info("Input SDF: %s", input_sdf)

    # Load docking config (optional; used as fallback for receptor path)
    docking_config: dict[str, Any] = {}
    docking_cfg_path = config.get(CFG_DOCKING)
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
        logger.info("Protein PDB: %s", protein_pdb)

    # Load molecules
    try:
        mols = load_molecules_from_sdf(input_sdf)
    except Exception as e:
        logger.error("Failed to load molecules: %s", e)
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
            mol, {"mol_idx", "sm_mol_idx", "s_sm_mol_idx", "name", "_Name"}
        )
        model_names.append(model_name or "")
        mol_idxs.append(mol_idx or "")

        gnina_min_aff.append(_get_prop_as_float(mol, "minimizedAffinity"))
        gnina_cnn_score.append(_get_prop_as_float(mol, "CNNscore"))
        gnina_cnn_aff.append(_get_prop_as_float(mol, "CNNaffinity"))

    # Fallback: try to fill missing model_name and mol_idx from ligands.csv
    ligands_csv = docking_dir / "ligands.csv"
    if ligands_csv.exists():
        try:
            lig_df = pd.read_csv(ligands_csv)
            if "mol_idx" in lig_df.columns and "model_name" in lig_df.columns:
                # Create lookup dictionaries by mol_idx and by name
                mol_idx_to_model_name = dict(
                    zip(lig_df["mol_idx"].astype(str), lig_df["model_name"].astype(str))
                )
                name_to_model_name = {}
                name_to_mol_idx = {}
                if "name" in lig_df.columns:
                    name_to_model_name = dict(
                        zip(
                            lig_df["name"].astype(str), lig_df["model_name"].astype(str)
                        )
                    )
                    name_to_mol_idx = dict(
                        zip(lig_df["name"].astype(str), lig_df["mol_idx"].astype(str))
                    )

                # Fill in missing values using mol_idx or name as lookup key
                for i, mol in enumerate(mols):
                    model_name = model_names[i]
                    mol_idx = mol_idxs[i]

                    # Try to get name from molecule if mol_idx is missing
                    mol_name = None
                    if not mol_idx:
                        mol_name = _get_first_prop_value(mol, {"name", "_Name"})

                    # Fill model_name if missing
                    if not model_name:
                        if mol_idx and mol_idx in mol_idx_to_model_name:
                            model_names[i] = mol_idx_to_model_name[mol_idx]
                        elif mol_name and mol_name in name_to_model_name:
                            model_names[i] = name_to_model_name[mol_name]

                    # Fill mol_idx if missing
                    if not mol_idx:
                        if mol_name and mol_name in name_to_mol_idx:
                            mol_idxs[i] = name_to_mol_idx[mol_name]
        except Exception as e:
            logger.debug("Could not load model_name/mol_idx from ligands.csv: %s", e)

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

    # Apply filters
    filters_applied: list[str] = []

    # Stage progress: map each enabled filter step to a 0-100 slot.
    sb_config = filter_config.get("search_box", {})
    ss_config = filter_config.get("shepherd_score", {})
    cd_config = filter_config.get("conformer_deviation", {})

    ref_path = ss_config.get("reference_ligand")

    def _step_progress(label: str, step_total: int):
        if reporter is None:
            return None

        def _progress(done: int, total: int) -> None:
            progress_total = total if total > 0 else step_total
            progress_done = max(0, min(done, progress_total))
            reporter.progress(
                progress_done, progress_total, message=f"DockingFilters: {label}"
            )

        return _progress

    agg_mode = filter_config.get("aggregation", {}).get("mode", "all")

    # Filter 0: Search-box containment (fast)
    search_total = len(mols)
    try:
        if reporter is not None:
            reporter.progress(0, search_total, message="DockingFilters: search_box")
        sb_df = apply_search_box_filter(mols, base_folder, docking_config, sb_config)
        results_df = results_df.merge(sb_df, on="mol_idx", how="left")
        filters_applied.append("search_box")
    except Exception as e:
        logger.error("Search-box filter failed: %s", e)
        results_df["pass_search_box"] = True
    finally:
        if reporter is not None:
            reporter.progress(
                search_total, search_total, message="DockingFilters: search_box"
            )

    # Optional optimization: under aggregation mode "all", if a pose fails search-box
    # containment it cannot pass the overall filter, so we can skip heavier checks.
    sb_short_circuit = bool(sb_config.get("short_circuit", True)) and agg_mode == "all"
    active_pose_indices = list(range(len(mols)))
    if sb_short_circuit and "pass_search_box" in results_df.columns:
        active_pose_indices = results_df.loc[
            results_df["pass_search_box"] == True, "mol_idx"  # noqa: E712
        ].tolist()

    # Filter 1: Pose Quality -- dispatch by backend
    if pq_config.get("enabled", True):
        pq_backend = pq_config.get("backend", "posebusters_fast")

        def _run_pose_quality(mols_active, progress_cb):
            if pq_backend == "posebusters_fast":
                return apply_posebusters_fast_filter(
                    mols_active, protein_pdb, pq_config, progress_cb=progress_cb
                )
            return apply_pose_quality_filter(mols_active, protein_pdb, pq_config)

        results_df = _run_single_filter(
            filter_name="pose_quality",
            reporter=reporter,
            step_progress_fn=_step_progress,
            mols=mols,
            active_pose_indices=active_pose_indices,
            results_df=results_df,
            filters_applied=filters_applied,
            pass_col="pass_pose_quality",
            run_fn=_run_pose_quality,
        )

    # Optional optimization: under aggregation mode "all", if a pose fails pose quality
    # it cannot pass the overall filter, so we can skip heavier checks (interactions,
    # conformer deviation, etc.) for those poses.
    pq_short_circuit = bool(pq_config.get("short_circuit", True)) and agg_mode == "all"
    if pq_short_circuit and "pass_pose_quality" in results_df.columns:
        pq_mask = results_df["pass_pose_quality"].fillna(False) == True  # noqa: E712
        if "pass_search_box" in results_df.columns:
            sb_mask = results_df["pass_search_box"].fillna(False) == True  # noqa: E712
            pq_mask = pq_mask & sb_mask
        active_pose_indices = results_df.loc[pq_mask, "mol_idx"].tolist()

    # Filter 2: Interactions
    if int_config.get("enabled", True):

        def _run_interactions(mols_active, _progress_cb):
            return apply_interaction_filter(mols_active, protein_pdb, int_config)

        results_df = _run_single_filter(
            filter_name="interactions",
            reporter=reporter,
            step_progress_fn=_step_progress,
            mols=mols,
            active_pose_indices=active_pose_indices,
            results_df=results_df,
            filters_applied=filters_applied,
            pass_col="pass_interactions",
            run_fn=_run_interactions,
        )

    # Filter 3: Shepherd-Score
    if ss_config.get("enabled", False):
        if ref_path and Path(ref_path).exists():
            ref_mol = None
            try:
                ref_mol = Chem.MolFromMolFile(str(ref_path))
            except Exception as e:
                logger.error("Shepherd-Score filter failed: %s", e)

            if ref_mol:

                def _run_shepherd(mols_active, progress_cb):
                    return apply_shepherd_score_filter(
                        mols_active, ref_mol, ss_config, progress_cb=progress_cb
                    )

                results_df = _run_single_filter(
                    filter_name="shepherd_score",
                    reporter=reporter,
                    step_progress_fn=_step_progress,
                    mols=mols,
                    active_pose_indices=active_pose_indices,
                    results_df=results_df,
                    filters_applied=filters_applied,
                    pass_col="pass_shepherd_score",
                    run_fn=_run_shepherd,
                )
            else:
                logger.warning("Failed to load reference molecule for Shepherd-Score")
                results_df["pass_shepherd_score"] = True
        else:
            logger.info("Shepherd-Score disabled (no reference ligand)")
            results_df["pass_shepherd_score"] = True

    # Filter 4: Conformer Deviation -- dispatch by backend
    if cd_config.get("enabled", True):
        cd_backend = _resolve_conformer_backend(cd_config)

        def _run_conformer(mols_active, progress_cb):
            if cd_backend == "symmetry_rmsd":
                return apply_symmetry_rmsd_filter(
                    mols_active, cd_config, progress_cb=progress_cb
                )
            return apply_conformer_deviation_filter(
                mols_active, cd_config, progress_cb=progress_cb
            )

        results_df = _run_single_filter(
            filter_name="conformer_deviation",
            reporter=reporter,
            step_progress_fn=_step_progress,
            mols=mols,
            active_pose_indices=active_pose_indices,
            results_df=results_df,
            filters_applied=filters_applied,
            pass_col="pass_conformer_deviation",
            run_fn=_run_conformer,
        )

    if reporter is not None:
        total_molecules = len(mols)
        reporter.progress(
            total_molecules, total_molecules, message="DockingFilters complete"
        )

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
    logger.info("Docking filters complete: %d/%d molecules passed", n_passed, n_total)
    logger.info("Filters applied: %s", ", ".join(filters_applied))

    interaction_reporting_enabled = _is_interaction_reporting_enabled(int_config)

    if interaction_reporting_enabled:
        _write_interaction_reporting_artifacts(
            results_df=results_df,
            output_dir=output_dir,
            interaction_filter_enabled=bool(int_config.get("enabled", True)),
        )
        logger.info("Saved interaction reporting artifacts to %s", output_dir)
    else:
        _cleanup_interaction_reporting_artifacts(output_dir)
        logger.info("Interaction reporting artifacts disabled by config")

    # Save results
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
