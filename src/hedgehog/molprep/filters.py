from __future__ import annotations

from rdkit import Chem


def _is_single_fragment(mol: Chem.Mol) -> bool:
    try:
        frags = Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False)
        return len(frags) <= 1
    except Exception:
        return False


def _has_radicals(mol: Chem.Mol) -> bool:
    return any(a.GetNumRadicalElectrons() > 0 for a in mol.GetAtoms())
    

def _has_isotopes(mol: Chem.Mol) -> bool:
    return any(a.GetIsotope() != 0 for a in mol.GetAtoms())


def _allowed_atoms_ok(mol: Chem.Mol, allowed: set[str]) -> bool:
    if not allowed:
        return True
    return all(a.GetSymbol() in allowed for a in mol.GetAtoms())
