"""Tests for docking molecule identity alias resolution."""

from pathlib import Path

import pandas as pd

from hedgehog.docking.identity import (
    build_canonical_mol_idx_map,
    build_smiles_lookup,
    dock_mol_idx_aliases,
    resolve_canonical_mol_idx,
)


def test_dock_mol_idx_aliases_for_numeric_float_ids():
    assert set(dock_mol_idx_aliases("9487.0")) == {"9487.0", "9487_0"}
    assert set(dock_mol_idx_aliases("9487_0")) == {"9487_0", "9487.0"}


def test_dock_mol_idx_aliases_preserve_non_numeric_ids():
    assert dock_mol_idx_aliases("LP-0001-00001") == ("LP-0001-00001",)


def test_build_smiles_lookup_supports_dock_side_alias(tmp_path: Path):
    ligands_csv = tmp_path / "ligands.csv"
    pd.DataFrame(
        {
            "smiles": ["c1ccccc1"],
            "name": ["9487.0"],
            "model_name": ["tgmdlm"],
            "mol_idx": ["9487.0"],
        }
    ).to_csv(ligands_csv, index=False)

    canonical_map = build_canonical_mol_idx_map(ligands_csv)
    smiles_lookup = build_smiles_lookup(ligands_csv)

    assert resolve_canonical_mol_idx("9487_0", canonical_map) == "9487.0"
    assert smiles_lookup["9487_0"] == "c1ccccc1"
