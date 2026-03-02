from pathlib import Path

from hedgehog._constants import TOOL_GNINA
from hedgehog.configs.logger import logger
from hedgehog.docking.receptor_prep import _resolve_receptor_path


def _resolve_path(path, base_dir):
    """Resolve path to absolute, using base_dir if relative."""
    path_obj = Path(path)
    if path_obj.is_absolute():
        return str(path_obj.resolve())
    return str((base_dir / path).resolve())


def _resolve_autobox_path(autobox_ligand, project_root):
    """Resolve autobox_ligand path to absolute."""
    path = Path(autobox_ligand)
    if path.is_absolute():
        return path if path.exists() else None

    candidate = (project_root / autobox_ligand).resolve()
    if candidate.exists():
        return candidate

    if "data/" in autobox_ligand:
        data_path = autobox_ligand[autobox_ligand.find("data/") :]
        candidate = (project_root / data_path).resolve()
        if candidate.exists():
            return candidate

    return None


def _pdb_atom_coordinates(pdb_path: Path) -> list[tuple[float, float, float]]:
    """Extract atom coordinates from a PDB file (ATOM/HETATM records only)."""
    coords: list[tuple[float, float, float]] = []
    try:
        for line in pdb_path.read_text(encoding="utf-8", errors="ignore").splitlines():
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except Exception:
                continue
            coords.append((x, y, z))
    except OSError:
        return []
    return coords


def _sdf_center(sdf_path: Path) -> tuple[float, float, float] | None:
    """Compute center of the first conformer in an SDF file (mean of atom positions)."""
    try:
        from rdkit import Chem
    except Exception:
        return None

    suppl = Chem.SDMolSupplier(str(sdf_path))
    mol = next((m for m in suppl if m is not None), None)
    if mol is None or mol.GetNumConformers() == 0:
        return None
    conf = mol.GetConformer()
    n = mol.GetNumAtoms()
    sx = sy = sz = 0.0
    for i in range(n):
        p = conf.GetAtomPosition(i)
        sx += float(p.x)
        sy += float(p.y)
        sz += float(p.z)
    return (sx / n, sy / n, sz / n)


def _min_distance_to_point(
    coords: list[tuple[float, float, float]], point: tuple[float, float, float]
) -> float | None:
    """Compute the minimum Euclidean distance from coords to a point."""
    if not coords:
        return None
    px, py, pz = point
    best = None
    for x, y, z in coords:
        dx = x - px
        dy = y - py
        dz = z - pz
        d2 = dx * dx + dy * dy + dz * dz
        if best is None or d2 < best:
            best = d2
    return (best or 0.0) ** 0.5


def _warn_if_autobox_far_from_receptor(cfg: dict, tool_name: str) -> None:
    """Warn if autobox reference ligand seems far away from the receptor coordinates.

    This typically indicates a mismatched coordinate frame between receptor PDB
    and autobox ligand (e.g., different reference, different prepared structure),
    resulting in docking running in the wrong location.
    """
    try:
        receptor_raw = cfg.get("receptor_pdb")
        if not receptor_raw:
            return

        receptor_path = _resolve_receptor_path(receptor_raw)
        if receptor_path is None:
            return

        tool_cfg = cfg.get(f"{tool_name}_config", {}) or {}
        autobox_ligand = tool_cfg.get("autobox_ligand") or cfg.get("autobox_ligand")
        if not autobox_ligand:
            return

        project_root = Path(__file__).parent.parent.parent.parent
        autobox_path = _resolve_autobox_path(str(autobox_ligand), project_root)
        if autobox_path is None:
            return

        center = _sdf_center(Path(autobox_path))
        if center is None:
            return

        protein_coords = _pdb_atom_coordinates(Path(receptor_path))
        min_dist = _min_distance_to_point(protein_coords, center)
        if min_dist is None:
            return

        warn_threshold = float(tool_cfg.get("autobox_receptor_distance_warn", 10.0))
        if min_dist > warn_threshold:
            logger.warning(
                "%s: Autobox reference ligand appears far from receptor (min dist %.2f A). "
                "This may indicate a wrong/mismatched autobox_ligand or receptor coordinate frame.",
                tool_name.upper(),
                min_dist,
            )
    except Exception:
        # Never fail docking due to a warning-only heuristic.
        return


def _count_box_warnings(log_path: Path) -> dict[str, int]:
    """Count GNINA/Vina 'outside box' warnings in a docking log."""
    import re

    patterns = ("outside box", "not within box")
    counts = {"lines": 0, "unique_molecules": 0}
    try:
        text = log_path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return counts

    mol_ids: set[str] = set()
    for line in text.splitlines():
        low = line.lower()
        if not any(p in low for p in patterns):
            continue
        counts["lines"] += 1
        m = re.match(r"^([^|]+?)\s*\|", line)
        if m:
            mol_id = m.group(1)
            if mol_id:
                mol_ids.add(mol_id.strip())
    counts["unique_molecules"] = len(mol_ids)
    return counts


def _gnina_zero_affinity_count(output_sdf: Path) -> tuple[int, int]:
    """Count how many poses have minimizedAffinity == 0.0 in GNINA SDF output."""
    try:
        from rdkit import Chem
    except Exception:
        return (0, 0)

    total = 0
    zero = 0
    suppl = Chem.SDMolSupplier(str(output_sdf))
    for mol in suppl:
        if mol is None:
            continue
        total += 1
        if mol.HasProp("minimizedAffinity"):
            try:
                if float(mol.GetProp("minimizedAffinity")) == 0.0:
                    zero += 1
            except Exception:
                continue
    return (zero, total)


def _emit_post_docking_warnings(
    tool_name: str,
    log_path: Path | None,
    output_sdf: Path | None = None,
) -> None:
    """Emit warnings about common docking failure modes based on tool logs/outputs."""
    try:
        if log_path and log_path.exists():
            counts = _count_box_warnings(log_path)
            if counts["lines"] > 0:
                logger.warning(
                    "%s log contains %d box warning lines across %d molecules. "
                    "This often indicates ligands started outside the search box or a misconfigured box.",
                    tool_name.upper(),
                    counts["lines"],
                    counts["unique_molecules"],
                )

        if tool_name.lower() == TOOL_GNINA and output_sdf and output_sdf.exists():
            zero, total = _gnina_zero_affinity_count(output_sdf)
            if total > 0 and zero > 0:
                logger.warning(
                    "GNINA output has minimizedAffinity == 0.0 for %d/%d poses. "
                    "If many, docking may have effectively failed (e.g., wrong box / no contacts).",
                    zero,
                    total,
                )
    except Exception:
        # Warning-only logic must never fail docking.
        return
