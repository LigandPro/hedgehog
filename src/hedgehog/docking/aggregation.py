from pathlib import Path

from hedgehog.configs.logger import logger


def _aggregate_docking_results(results_dir: Path, output_sdf: Path) -> int:
    """Aggregate per-molecule docking results into a single SDF file.

    Args:
        results_dir: Directory containing per-molecule result SDFs (*_out.sdf)
        output_sdf: Path to write aggregated output SDF

    Returns:
        Number of molecules successfully aggregated
    """
    try:
        from rdkit import Chem
    except ImportError as err:
        raise RuntimeError("RDKit not available for result aggregation") from err

    output_sdf.parent.mkdir(parents=True, exist_ok=True)

    result_files = sorted(results_dir.glob("*_out.sdf"))
    if not result_files:
        logger.warning("No result files found in %s", results_dir)
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

    writer = Chem.SDWriter(str(output_sdf))
    count = 0

    for result_file in result_files:
        try:
            best_pose = _pick_best_pose(result_file)
            if best_pose is not None:
                writer.write(best_pose)
                count += 1
        except Exception as e:
            logger.warning("Failed to read result file %s: %s", result_file, e)
            continue

    writer.close()
    logger.info(
        "Aggregated %d best poses (single pose per molecule) from %d result files into %s",
        count,
        len(result_files),
        output_sdf,
    )
    return count
