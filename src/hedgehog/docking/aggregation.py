import re
from pathlib import Path

from hedgehog.configs.logger import logger

_PER_MOLECULE_RESULT_RE = re.compile(r"^\d{6}_(.+)_out$")


def _load_rdkit():
    try:
        from rdkit import Chem
    except ImportError as err:
        raise RuntimeError("RDKit not available for result aggregation") from err
    return Chem


def _aggregate_sdf_files(
    result_files: list[Path],
    output_sdf: Path,
    *,
    pick_best_pose: bool,
    name_from_filename: bool = False,
) -> int:
    """Aggregate SDF files into a single output SDF."""
    Chem = _load_rdkit()

    output_sdf.parent.mkdir(parents=True, exist_ok=True)

    if not result_files:
        logger.warning("No result files found for aggregation into %s", output_sdf)
        return 0

    def _extract_pose_affinity(mol) -> float | None:
        for prop_name in ("minimizedAffinity", "affinity", "score"):
            if mol.HasProp(prop_name):
                try:
                    return float(mol.GetProp(prop_name))
                except Exception:
                    continue
        return None

    def _pick_best_pose(result_file: Path):
        best_mol = None
        best_affinity = float("inf")
        for mol in Chem.SDMolSupplier(str(result_file)):
            if mol is None:
                continue
            affinity = _extract_pose_affinity(mol)
            affinity_sort = affinity if affinity is not None else float("inf")
            if best_mol is None or affinity_sort < best_affinity:
                best_mol = mol
                best_affinity = affinity_sort
        return best_mol

    def _restore_per_molecule_source_id(mol, result_file: Path) -> None:
        """Restore stable source molecule identity from per-molecule result filename.

        Per-molecule docking writes files as:
          <6-digit-split-prefix>_<source-mol-id>_out.sdf
        Example:
          000123_t-42_out.sdf -> source id "t-42"

        We persist this id into canonical SDF fields so downstream stages can
        reliably map back to original input rows even if docking tools mutate
        or drop name-like properties.
        """
        match = _PER_MOLECULE_RESULT_RE.match(result_file.stem)
        if not match:
            return

        source_mol_idx = match.group(1).strip()
        if not source_mol_idx:
            return

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
                supplier = Chem.SDMolSupplier(str(result_file), removeHs=False)
                output_mol = next((mol for mol in supplier if mol is not None), None)
            if output_mol is not None:
                _restore_per_molecule_source_id(output_mol, result_file)
                if name_from_filename:
                    output_mol.SetProp("_Name", result_file.stem)
                    if not output_mol.HasProp("mol_idx"):
                        output_mol.SetProp("mol_idx", result_file.stem)
                writer.write(output_mol)
                count += 1
        except Exception as e:
            logger.warning("Failed to read result file %s: %s", result_file, e)
            continue

    writer.close()
    return count


def _aggregate_docking_results(results_dir: Path, output_sdf: Path) -> int:
    """Aggregate per-molecule docking results into a single SDF file."""
    result_files = sorted(results_dir.glob("*_out.sdf"))
    count = _aggregate_sdf_files(result_files, output_sdf, pick_best_pose=True)
    logger.info(
        "Aggregated %d best poses (single pose per molecule) from %d result files into %s",
        count,
        len(result_files),
        output_sdf,
    )
    return count


def _aggregate_matcha_results(best_poses_dir: Path, output_sdf: Path) -> int:
    """Aggregate Matcha best-pose SDF files into a single output SDF."""
    result_files = sorted(best_poses_dir.glob("*.sdf"))
    count = _aggregate_sdf_files(
        result_files,
        output_sdf,
        pick_best_pose=False,
        name_from_filename=True,
    )
    logger.info(
        "Aggregated %d Matcha best poses from %d files into %s",
        count,
        len(result_files),
        output_sdf,
    )
    return count
