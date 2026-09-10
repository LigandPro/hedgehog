"""Helpers for resolving molecule identity across docking filesystem encodings."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def dock_mol_idx_aliases(mol_idx: str) -> tuple[str, ...]:
    """Return lookup aliases for a molecule id across docking path sanitization.

      Per-molecule docking files sanitize ``_Name`` by replacing characters outside
    ``[A-Za-z0-9_-]`` with underscores. Numeric ``mol_idx`` values like ``9487.0``
      therefore appear as ``9487_0`` in SDF filenames and aggregated pose metadata.
    """
    key = str(mol_idx).strip()
    if not key:
        return ()

    aliases = {key}
    if "." in key:
        aliases.add(key.replace(".", "_"))

    if "_" in key:
        head, tail = key.rsplit("_", 1)
        numeric_head = head.replace(".", "")
        if tail.isdigit() and numeric_head.isdigit():
            aliases.add(f"{head}.{tail}")

    return tuple(aliases)


def build_canonical_mol_idx_map(ligands_csv: Path | None) -> dict[str, str]:
    """Map dock filesystem ids and aliases to canonical ``mol_idx`` values."""
    if ligands_csv is None or not ligands_csv.exists():
        return {}

    try:
        lig_df = pd.read_csv(ligands_csv)
    except Exception:
        return {}

    if "mol_idx" not in lig_df.columns:
        return {}

    canonical_map: dict[str, str] = {}
    name_col = lig_df["name"].astype(str) if "name" in lig_df.columns else None

    for row_idx, mol_idx in enumerate(lig_df["mol_idx"].astype(str)):
        canonical = str(mol_idx).strip()
        if not canonical:
            continue

        keys = set(dock_mol_idx_aliases(canonical))
        if name_col is not None:
            keys.update(dock_mol_idx_aliases(str(name_col.iloc[row_idx]).strip()))

        for key in keys:
            if key:
                canonical_map[key] = canonical

    return canonical_map


def resolve_canonical_mol_idx(
    dock_id: str,
    canonical_map: dict[str, str] | None = None,
) -> str:
    """Resolve a dock-side molecule id to the canonical ``mol_idx`` string."""
    key = str(dock_id).strip()
    if not key:
        return key

    if canonical_map:
        if key in canonical_map:
            return canonical_map[key]
        for alias in dock_mol_idx_aliases(key):
            if alias in canonical_map:
                return canonical_map[alias]

    for alias in dock_mol_idx_aliases(key):
        if alias != key:
            return alias
    return key


def build_smiles_lookup(ligands_csv: Path) -> dict[str, str]:
    """Build a SMILES lookup keyed by canonical and dock-side ``mol_idx`` aliases."""
    canonical_map = build_canonical_mol_idx_map(ligands_csv)
    if not canonical_map:
        return {}

    try:
        lig_df = pd.read_csv(ligands_csv)
    except Exception:
        return {}

    if "mol_idx" not in lig_df.columns or "smiles" not in lig_df.columns:
        return {}

    smiles_by_canonical = dict(
        zip(lig_df["mol_idx"].astype(str).str.strip(), lig_df["smiles"].astype(str))
    )

    lookup: dict[str, str] = {}
    for alias, canonical in canonical_map.items():
        smiles = smiles_by_canonical.get(canonical)
        if smiles is not None:
            lookup[alias] = smiles
    return lookup
