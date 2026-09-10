import re
from pathlib import Path

import pandas as pd

from hedgehog.configs.logger import logger
from hedgehog.docking.identity import (
    build_canonical_mol_idx_map,
    resolve_canonical_mol_idx,
)

_PER_MOLECULE_RESULT_RE = re.compile(r"^\d{6}_(.+)_out$")


def _load_rdkit():
    try:
        from rdkit import Chem
    except ImportError as err:
        raise RuntimeError("RDKit not available for result aggregation") from err
    return Chem


def _read_pose_float_props(mol, property_names) -> float | None:
    for prop_name in dict.fromkeys(name for name in property_names if name):
        if not mol.HasProp(prop_name):
            continue
        try:
            return float(mol.GetProp(prop_name))
        except Exception:
            continue
    return None


def _extract_pose_affinity_value(mol) -> float | None:
    """Pick-best ranking uses docking affinity props only, never BALMUS metadata."""
    return _read_pose_float_props(mol, ("minimizedAffinity", "affinity", "score"))


def _extract_pose_score_value(mol, preferred_property: str = "") -> float | None:
    """Read a lower-is-better docking score without conflating score scales."""
    if preferred_property:
        property_names = [preferred_property]
    else:
        property_names = []
        if mol.HasProp("source_score_property"):
            property_names.append(mol.GetProp("source_score_property").strip())
        property_names.extend(("minimizedAffinity", "affinity", "score"))
    return _read_pose_float_props(mol, property_names)


def _build_matcha_metadata_map(ligands_csv: Path | None) -> dict[str, tuple[str, str]]:
    if ligands_csv is None or not ligands_csv.exists():
        return {}

    try:
        df = pd.read_csv(ligands_csv)
    except Exception:
        return {}

    if "mol_idx" not in df.columns:
        return {}

    metadata: dict[str, tuple[str, str]] = {}
    for _, row in df.iterrows():
        mol_idx = str(row.get("mol_idx", "")).strip()
        if not mol_idx:
            continue
        model_name = str(row.get("model_name", "")).strip()
        # Matcha path sanitizer converts "/" to "_", so index both forms.
        for key in (mol_idx, mol_idx.replace("/", "_")):
            if key and key not in metadata:
                metadata[key] = (mol_idx, model_name)
    return metadata


def _read_docking_sdf_molecules(result_file: Path, Chem):
    """Read docking poses and repair a known charged-sulfoxide bond encoding."""
    from rdkit import rdBase

    with rdBase.BlockLogs():
        supplier = Chem.SDMolSupplier(
            str(result_file),
            removeHs=False,
            sanitize=False,
        )
        for mol in supplier:
            if mol is None:
                continue

            for bond in mol.GetBonds():
                if bond.GetBondType() != Chem.BondType.DOUBLE:
                    continue
                begin = bond.GetBeginAtom()
                end = bond.GetEndAtom()
                charged_sulfoxide = (
                    begin.GetSymbol() == "S"
                    and begin.GetFormalCharge() > 0
                    and end.GetSymbol() == "O"
                    and end.GetFormalCharge() < 0
                ) or (
                    end.GetSymbol() == "S"
                    and end.GetFormalCharge() > 0
                    and begin.GetSymbol() == "O"
                    and begin.GetFormalCharge() < 0
                )
                if charged_sulfoxide:
                    bond.SetBondType(Chem.BondType.SINGLE)

            try:
                Chem.SanitizeMol(mol)
            except Exception as exc:
                logger.warning(
                    "Failed to sanitize docking result %s: %s",
                    result_file,
                    exc,
                )
                continue
            yield mol


def _aggregate_sdf_files(
    result_files: list[Path],
    output_sdf: Path,
    *,
    pick_best_pose: bool,
    model_name: str | None = None,
    canonical_map: dict[str, str] | None = None,
) -> int:
    """Aggregate SDF files into a single output SDF."""
    Chem = _load_rdkit()

    output_sdf.parent.mkdir(parents=True, exist_ok=True)

    if not result_files:
        logger.warning("No result files found for aggregation into %s", output_sdf)
        return 0

    def _pick_best_pose(result_file: Path):
        best_mol = None
        best_affinity = float("inf")
        for mol in _read_docking_sdf_molecules(result_file, Chem):
            affinity = _extract_pose_affinity_value(mol)
            affinity_sort = affinity if affinity is not None else float("inf")
            if best_mol is None or affinity_sort < best_affinity:
                best_mol = mol
                best_affinity = affinity_sort
        return best_mol

    def _restore_per_molecule_source_id(
        mol,
        result_file: Path,
        canonical_map: dict[str, str] | None = None,
    ) -> None:
        match = _PER_MOLECULE_RESULT_RE.match(result_file.stem)
        if not match:
            return

        dock_id = match.group(1).strip()
        if not dock_id:
            return

        source_mol_idx = resolve_canonical_mol_idx(dock_id, canonical_map)

        mol.SetProp("_Name", source_mol_idx)
        mol.SetProp("mol_idx", source_mol_idx)
        mol.SetProp("source_mol_idx", source_mol_idx)

    writer = Chem.SDWriter(str(output_sdf))
    count = 0

    for result_file in result_files:
        try:
            if pick_best_pose:
                output_mol = _pick_best_pose(result_file)
            else:
                output_mol = next(
                    _read_docking_sdf_molecules(result_file, Chem),
                    None,
                )
            if output_mol is not None:
                _restore_per_molecule_source_id(
                    output_mol, result_file, canonical_map=canonical_map
                )
                if model_name:
                    output_mol.SetProp("model_name", model_name)
                writer.write(output_mol)
                count += 1
        except Exception as e:
            logger.warning("Failed to read result file %s: %s", result_file, e)
            continue

    writer.close()
    return count


def _aggregate_docking_results(
    results_dir: Path,
    output_sdf: Path,
    model_name: str | None = None,
    ligands_csv: Path | None = None,
) -> int:
    """Aggregate per-molecule docking results into a single SDF file."""
    canonical_map = build_canonical_mol_idx_map(ligands_csv)
    result_files = sorted(results_dir.glob("*_out.sdf"))
    count = _aggregate_sdf_files(
        result_files,
        output_sdf,
        pick_best_pose=True,
        model_name=model_name,
        canonical_map=canonical_map,
    )
    logger.info(
        "Aggregated %d best poses (single pose per molecule) from %d result files into %s",
        count,
        len(result_files),
        output_sdf,
    )
    return count


def _aggregate_matcha_results(
    best_poses_dir: Path, output_sdf: Path, ligands_csv: Path | None = None
) -> int:
    """Aggregate Matcha best-pose SDF files into a single output SDF."""
    Chem = _load_rdkit()
    result_files = sorted(best_poses_dir.glob("*.sdf"))
    if not result_files:
        logger.warning("No Matcha best pose files found in %s", best_poses_dir)
        return 0

    metadata_map = _build_matcha_metadata_map(ligands_csv)
    output_sdf.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(output_sdf))
    count = 0

    for result_file in result_files:
        try:
            output_mol = next(
                _read_docking_sdf_molecules(result_file, Chem),
                None,
            )
            if output_mol is None:
                continue

            stem = result_file.stem
            output_mol.SetProp("_Name", stem)

            mapped = metadata_map.get(stem)
            if mapped:
                mapped_mol_idx, mapped_model_name = mapped
                output_mol.SetProp("mol_idx", mapped_mol_idx)
                if mapped_model_name:
                    output_mol.SetProp("model_name", mapped_model_name)
            elif not output_mol.HasProp("mol_idx"):
                output_mol.SetProp("mol_idx", stem)

            if output_mol.HasProp("minimizedAffinity"):
                output_mol.SetProp("source_score_property", "minimizedAffinity")

            writer.write(output_mol)
            count += 1
        except Exception as e:
            logger.warning("Failed to read Matcha result file %s: %s", result_file, e)
            continue

    writer.close()
    logger.info(
        "Aggregated %d Matcha best poses from %d files into %s",
        count,
        len(result_files),
        output_sdf,
    )
    return count


def _collect_docking_stage_results(
    ligands_dir: Path,
    selected_tools: list[str],
    tool_outputs: dict[str, Path],
    score_thresholds: dict | None = None,
) -> tuple[int, int]:
    """Collect tool poses and apply optional stage-level docking score cutoffs.

    Without score thresholds, a ligand passes when any selected tool produced a
    valid pose (legacy behavior). With aligned thresholds, a ligand must have a
    score and pass the configured cutoff for every calibrated tool.
    """
    Chem = _load_rdkit()
    ligands_csv = ligands_dir / "ligands.csv"
    input_csv = ligands_dir / "input_molecules.csv"
    ligands_df = pd.read_csv(input_csv if input_csv.exists() else ligands_csv)
    identity_cols = ["smiles", "model_name", "mol_idx"]
    missing_cols = [col for col in identity_cols if col not in ligands_df.columns]
    if missing_cols:
        raise ValueError(
            "Docking ligands.csv is missing required columns: "
            + ", ".join(missing_cols)
        )

    ligands_df = ligands_df.drop_duplicates(subset=["mol_idx"], keep="first").copy()
    ligands_df["mol_idx"] = ligands_df["mol_idx"].astype(str)
    canonical_map = build_canonical_mol_idx_map(ligands_csv)
    ligand_ids = set(ligands_df["mol_idx"])
    model_lookup = dict(
        zip(ligands_df["mol_idx"], ligands_df["model_name"].astype(str))
    )
    smiles_lookup = dict(zip(ligands_df["mol_idx"], ligands_df["smiles"].astype(str)))

    configured_thresholds: dict[str, dict] = {}
    raw_thresholds = score_thresholds if isinstance(score_thresholds, dict) else {}
    for tool in selected_tools:
        spec = raw_thresholds.get(tool)
        if isinstance(spec, dict) and ("min" in spec or "max" in spec):
            configured_thresholds[tool] = spec
            continue
        # Read old generated GNINA affinity configs without changing their files.
        legacy = raw_thresholds.get(f"{tool}_minimizedAffinity")
        if isinstance(legacy, dict) and ("min" in legacy or "max" in legacy):
            configured_thresholds[tool] = legacy

    successful_by_tool: dict[str, set[str]] = {tool: set() for tool in selected_tools}
    scores_by_tool: dict[str, dict[str, float]] = {tool: {} for tool in selected_tools}
    combined_sdf = ligands_dir / "docking_out.sdf"
    writer = Chem.SDWriter(str(combined_sdf))
    pose_count = 0

    try:
        for tool in selected_tools:
            output_sdf = tool_outputs.get(tool)
            if output_sdf is None or not output_sdf.exists():
                logger.warning(
                    "No %s pose output available for stage aggregation", tool
                )
                continue

            for mol in Chem.SDMolSupplier(str(output_sdf), removeHs=False):
                if mol is None:
                    continue
                raw_id = ""
                for prop_name in ("source_mol_idx", "mol_idx", "_Name"):
                    if mol.HasProp(prop_name):
                        raw_id = mol.GetProp(prop_name).strip()
                        if raw_id:
                            break
                canonical_id = resolve_canonical_mol_idx(raw_id, canonical_map)
                if canonical_id not in ligand_ids:
                    logger.warning(
                        "Ignoring %s pose with unresolved molecule id %r",
                        tool,
                        raw_id,
                    )
                    continue

                successful_by_tool[tool].add(canonical_id)
                score_property = str(
                    configured_thresholds.get(tool, {}).get("score_property", "")
                ).strip()
                score = _extract_pose_score_value(mol, score_property)
                if score is not None:
                    previous = scores_by_tool[tool].get(canonical_id)
                    if previous is None or score < previous:
                        scores_by_tool[tool][canonical_id] = score
                mol.SetProp("_Name", canonical_id)
                mol.SetProp("mol_idx", canonical_id)
                mol.SetProp("source_mol_idx", canonical_id)
                mol.SetProp("model_name", model_lookup[canonical_id])
                mol.SetProp("input_smiles", smiles_lookup[canonical_id])
                mol.SetProp("docking_tool", tool)
                writer.write(mol)
                pose_count += 1
    finally:
        writer.close()

    metrics_df = ligands_df[identity_cols].copy()
    pass_cols: list[str] = []
    threshold_pass_cols: list[str] = []
    for tool in selected_tools:
        has_pose_col = f"has_pose_{tool}"
        score_col = f"score_{tool}"
        pass_col = f"pass_{tool}"
        pass_cols.append(pass_col)
        metrics_df[has_pose_col] = metrics_df["mol_idx"].isin(successful_by_tool[tool])
        metrics_df[score_col] = metrics_df["mol_idx"].map(scores_by_tool[tool])

        spec = configured_thresholds.get(tool)
        if spec is None:
            metrics_df[pass_col] = metrics_df[has_pose_col]
            continue
        passed = metrics_df[score_col].notna()
        if spec.get("min") is not None:
            passed &= metrics_df[score_col] >= float(spec["min"])
        if spec.get("max") is not None:
            passed &= metrics_df[score_col] <= float(spec["max"])
        metrics_df[pass_col] = passed
        threshold_pass_cols.append(pass_col)

    if threshold_pass_cols:
        metrics_df["pass"] = metrics_df[threshold_pass_cols].all(axis=1)
    else:
        metrics_df["pass"] = metrics_df[pass_cols].any(axis=1) if pass_cols else False
    metrics_df["docking_tools"] = metrics_df["mol_idx"].map(
        lambda mol_idx: ",".join(
            tool for tool in selected_tools if mol_idx in successful_by_tool[tool]
        )
    )
    metrics_df.to_csv(ligands_dir / "docking_results.csv", index=False)

    output_cols = [*identity_cols, "docking_tools"]
    filtered_df = metrics_df.loc[metrics_df["pass"], output_cols].copy()
    failed_df = metrics_df.loc[~metrics_df["pass"], output_cols].copy()
    filtered_df.to_csv(ligands_dir / "filtered_molecules.csv", index=False)
    failed_df.to_csv(ligands_dir / "failed_molecules.csv", index=False)

    if threshold_pass_cols:
        logger.info(
            "Applied docking score thresholds for tools: %s",
            ", ".join(configured_thresholds),
        )
    logger.info(
        "Docking aggregation complete: %d/%d molecules passed using %d tool poses",
        len(filtered_df),
        len(metrics_df),
        pose_count,
    )
    return len(filtered_df), len(failed_df)
