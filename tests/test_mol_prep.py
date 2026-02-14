"""Tests for MolPrep stage (Datamol standardization)."""

from pathlib import Path

import datamol as dm
import pandas as pd

from hedgehog.stages.molPrep.utils import run_mol_prep


def _cfg():
    """Return a minimal MolPrep config dict."""
    return {
        "columns": {
            "smiles": "smiles",
            "model_name": "model_name",
            "mol_idx": "mol_idx",
            "smiles_raw": "smiles_raw",
        },
        "steps": {
            "to_mol": {
                "ordered": True,
                "sanitize": False,
                "allow_cxsmiles": True,
                "strict_cxsmiles": True,
                "remove_hs": True,
            },
            "fix_mol": {
                "enabled": True,
                "n_iter": 1,
                "remove_singleton": True,
                "largest_only": False,
            },
            "sanitize_mol": {"enabled": True},
            "remove_salts_solvents": {
                "enabled": True,
                "defn_data": None,
                "defn_format": "smarts",
                "dont_remove_everything": True,
                "sanitize": True,
            },
            "keep_largest_fragment": True,
            "standardize_mol": {
                "enabled": True,
                "disconnect_metals": True,
                "normalize": True,
                "reionize": True,
                "uncharge": True,
                "stereo": True,
            },
            "remove_stereochemistry": True,
            "standardize_smiles": {"enabled": True},
        },
        "filters": {
            "allowed_atoms": ["C", "N", "O", "S", "F", "Cl", "Br", "I", "P", "H"],
            "reject_radicals": True,
            "require_neutral": True,
            "reject_isotopes": True,
            "require_single_fragment": True,
        },
        "output": {"write_duplicates_removed": True},
    }


def test_salts_and_fragments_removed(tmp_path: Path):
    """SMILES with multiple fragments should become single-fragment."""
    df = pd.DataFrame(
        {
            "smiles": ["CCO.Cl"],  # ethanol + chlorine fragment
            "model_name": ["test"],
            "mol_idx": ["LP-0001-00001"],
        }
    )
    out = run_mol_prep(df, _cfg(), tmp_path)
    assert len(out) == 1
    assert "." not in out["smiles"].iloc[0]


def test_charges_removed(tmp_path: Path):
    """Protonated amine should be neutralized."""
    df = pd.DataFrame(
        {
            "smiles": ["C[NH+](C)C"],
            "model_name": ["test"],
            "mol_idx": ["LP-0001-00001"],
        }
    )
    out = run_mol_prep(df, _cfg(), tmp_path)
    assert len(out) == 1
    smi = out["smiles"].iloc[0]
    mol = dm.to_mol(smi, sanitize=True)
    assert mol is not None
    assert not any(a.GetFormalCharge() != 0 for a in mol.GetAtoms())


def test_stereochemistry_removed(tmp_path: Path):
    """Chiral SMILES should lose stereochemical markers."""
    df = pd.DataFrame(
        {
            "smiles": ["C[C@H](O)F"],
            "model_name": ["test"],
            "mol_idx": ["LP-0001-00001"],
        }
    )
    out = run_mol_prep(df, _cfg(), tmp_path)
    assert len(out) == 1
    assert "@" not in out["smiles"].iloc[0]


def test_invalid_smiles_reported(tmp_path: Path):
    """Invalid SMILES should be written to failed_molecules.csv with parse_failed."""
    df = pd.DataFrame(
        {
            "smiles": ["CCO", "invalid_smiles"],
            "model_name": ["test", "test"],
            "mol_idx": ["LP-0001-00001", "LP-0001-00002"],
        }
    )
    out = run_mol_prep(df, _cfg(), tmp_path)
    assert len(out) == 1

    failed = pd.read_csv(tmp_path / "failed_molecules.csv")
    assert "parse_failed" in failed["reason"].astype(str).tolist()
