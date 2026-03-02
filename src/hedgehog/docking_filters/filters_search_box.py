"""Search box containment filter for docking poses."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem

from hedgehog.configs.logger import logger
from hedgehog.docking_filters.types import _make_passthrough_df, _resolve_existing_path


def _autobox_from_reference_ligand(
    base_folder: Path, docking_cfg: dict[str, Any]
) -> tuple[tuple[float, float, float], tuple[float, float, float]] | None:
    """Build a Vina/GNINA-style search box from autobox_ligand + autobox_add.

    Returns:
        (center_xyz, size_xyz) or None if insufficient config.
    """
    tool_cfg = docking_cfg.get("gnina_config") or docking_cfg.get("smina_config") or {}
    autobox_ligand = tool_cfg.get("autobox_ligand") or docking_cfg.get("autobox_ligand")
    if not autobox_ligand:
        return None

    ref_path = _resolve_existing_path(base_folder, autobox_ligand)
    if not ref_path.exists():
        return None

    suppl = Chem.SDMolSupplier(str(ref_path))
    ref_mol = next((m for m in suppl if m is not None), None)
    if ref_mol is None or ref_mol.GetNumConformers() == 0:
        return None

    conf = ref_mol.GetConformer()
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    for i in range(ref_mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        xs.append(float(p.x))
        ys.append(float(p.y))
        zs.append(float(p.z))

    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    min_z, max_z = min(zs), max(zs)

    add = float(tool_cfg.get("autobox_add", 0.0) or 0.0)
    min_x -= add
    max_x += add
    min_y -= add
    max_y += add
    min_z -= add
    max_z += add

    center = ((min_x + max_x) / 2, (min_y + max_y) / 2, (min_z + max_z) / 2)
    size = (max_x - min_x, max_y - min_y, max_z - min_z)
    return center, size


def apply_search_box_filter(
    mols: list[Chem.Mol],
    base_folder: Path,
    docking_cfg: dict[str, Any],
    config: dict[str, Any],
) -> pd.DataFrame:
    """Fast filter: ensure docked poses are inside the configured search box.

    Supports explicit center/size, or autobox_ligand-based boxes.

    Returns:
        DataFrame with columns: mol_idx, frac_atoms_outside_box, pass_search_box
    """
    enabled = bool(config.get("enabled", True))
    if not enabled:
        return _make_passthrough_df(
            len(mols), {"frac_atoms_outside_box": 0.0, "pass_search_box": True}
        )

    # Prefer explicit center/size if present, otherwise fall back to autobox reference.
    center = docking_cfg.get("center")
    size = docking_cfg.get("size")
    if not (
        isinstance(center, (list, tuple))
        and len(center) >= 3
        and isinstance(size, (list, tuple))
        and len(size) >= 3
    ):
        box = _autobox_from_reference_ligand(base_folder, docking_cfg)
        if box is None:
            logger.warning(
                "Search-box filter enabled but no box could be resolved; passing all poses."
            )
            return _make_passthrough_df(
                len(mols), {"frac_atoms_outside_box": 0.0, "pass_search_box": True}
            )
        center, size = box

    cx, cy, cz = float(center[0]), float(center[1]), float(center[2])
    sx, sy, sz = float(size[0]), float(size[1]), float(size[2])
    hx, hy, hz = sx / 2.0, sy / 2.0, sz / 2.0

    max_outside = float(config.get("max_outside_fraction", 0.0) or 0.0)

    frac_outside: list[float] = []
    passed: list[bool] = []

    for mol in mols:
        if mol is None or mol.GetNumConformers() == 0:
            frac_outside.append(1.0)
            passed.append(False)
            continue
        conf = mol.GetConformer()
        n_atoms = mol.GetNumAtoms()
        outside = 0
        for i in range(n_atoms):
            p = conf.GetAtomPosition(i)
            x, y, z = float(p.x), float(p.y), float(p.z)
            if (
                x < cx - hx
                or x > cx + hx
                or y < cy - hy
                or y > cy + hy
                or z < cz - hz
                or z > cz + hz
            ):
                outside += 1
        frac = outside / n_atoms if n_atoms else 1.0
        frac_outside.append(frac)
        passed.append(frac <= max_outside)

    df = pd.DataFrame(
        {
            "mol_idx": range(len(mols)),
            "frac_atoms_outside_box": frac_outside,
            "pass_search_box": passed,
        }
    )
    logger.info(
        "Search-box filter: %d/%d passed (max_outside_fraction=%.3f)",
        int(df["pass_search_box"].sum()),
        len(df),
        max_outside,
    )
    return df
