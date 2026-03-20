"""Input loading and property extraction for docking filters."""

from __future__ import annotations

from pathlib import Path

from rdkit import Chem

from hedgehog.configs.logger import logger


def load_molecules_from_sdf(sdf_path: str | Path) -> list[Chem.Mol]:
    """Load molecules from SDF file."""
    sdf_path = Path(sdf_path)
    if not sdf_path.exists():
        raise FileNotFoundError(f"SDF file not found: {sdf_path}")

    suppl = Chem.SDMolSupplier(str(sdf_path))
    mols = [mol for mol in suppl if mol is not None]
    logger.info("Loaded %d molecules from %s", len(mols), sdf_path)
    return mols


def _get_first_prop_value(
    mol: Chem.Mol, canonical_names: set[str] | tuple[str, ...] | list[str]
) -> str | None:
    """Return the first SDF property value whose (normalized) key matches any canonical name."""
    ordered_names = tuple(canonical_names)
    canonical_lookup = set(ordered_names)

    # Fast path: exact keys
    for name in ordered_names:
        if mol.HasProp(name):
            return mol.GetProp(name)

    # Some toolchains escape underscores in SDF property names (e.g. "s_sm_model\\_name").
    for prop in mol.GetPropNames():
        normalized = prop.replace("\\", "")
        if normalized in canonical_lookup:
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
