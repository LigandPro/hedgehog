"""MCE-18 (Molecular Complexity Estimation) descriptor calculation.

Reference: J. Med. Chem. 2019, DOI: 10.1021/acs.jmedchem.9b00004
"Are We Opening the Door to a New Era of Medicinal Chemistry or
Being Collapsed to a Chemical Singularity?"

Formula: MCE-18 = [AR + NAR + CHIRAL + SPIRO + (SP3 + Cyc - Acyc)/(1 + SP3)] × Q1

Components:
- AR: Presence of aromatic/heteroaromatic rings (0 or 1)
- NAR: Presence of aliphatic/heteroaliphatic rings (0 or 1)
- CHIRAL: Presence of chiral centers (0 or 1)
- SPIRO: Presence of spiro centers (0 or 1)
- SP3: Fraction of sp³ carbons (0-1)
- Cyc: Fraction of cyclic sp³ carbons (0-1)
- Acyc: Fraction of acyclic sp³ carbons (0-1)
- Q1: Normalized quadratic index from adjacency matrix

Typical range: 0-56 for drug-like molecules.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def compute_mce18(mol: Chem.Mol | None) -> float | None:
    """Compute MCE-18 descriptor for a molecule.

    Args:
        mol: RDKit molecule object

    Returns:
        MCE-18 value (typically 0-56 for drug-like molecules),
        or None if calculation fails
    """
    if mol is None:
        return None

    try:
        # Binary structural features
        ar = _has_aromatic_ring(mol)
        nar = _has_aliphatic_ring(mol)
        chiral = _has_chiral_center(mol)
        spiro = _has_spiro_center(mol)

        # sp³ fractions
        sp3 = _calc_sp3_fraction(mol)
        cyc = _calc_cyclic_sp3_fraction(mol)
        acyc = _calc_acyclic_sp3_fraction(mol)

        # Normalized quadratic index
        q1 = _calc_q1_index(mol)

        # MCE-18 formula
        # Part 1: Sum of binary features
        part1 = ar + nar + chiral + spiro

        # Part 2: sp³ complexity ratio
        part2 = sp3 + cyc - acyc

        # Part 3: Normalization factor
        part3 = 1 + sp3

        # Final MCE-18
        mce18 = (part1 + part2 / part3) * q1

        return round(mce18, 4)

    except Exception:
        return None


def _has_aromatic_ring(mol: Chem.Mol) -> int:
    """Check if molecule has aromatic or heteroaromatic rings.

    Returns:
        1 if aromatic ring present, 0 otherwise
    """
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            return 1
    return 0


def _has_aliphatic_ring(mol: Chem.Mol) -> int:
    """Check if molecule has aliphatic or heteroaliphatic rings.

    Returns:
        1 if non-aromatic ring present, 0 otherwise
    """
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if not any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            return 1
    return 0


def _has_chiral_center(mol: Chem.Mol) -> int:
    """Check if molecule has chiral centers.

    Returns:
        1 if chiral center present, 0 otherwise
    """
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    return 1 if chiral_centers else 0


def _has_spiro_center(mol: Chem.Mol) -> int:
    """Check if molecule has spiro centers.

    A spiro center is an atom shared by two rings that doesn't
    share any bonds between the rings.

    Returns:
        1 if spiro center present, 0 otherwise
    """
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    if len(atom_rings) < 2:
        return 0

    # Find atoms that belong to multiple rings
    atom_ring_count: dict[int, int] = {}
    for ring in atom_rings:
        for atom_idx in ring:
            atom_ring_count[atom_idx] = atom_ring_count.get(atom_idx, 0) + 1

    # Check each atom in multiple rings for spiro character
    for atom_idx, count in atom_ring_count.items():
        if count >= 2:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Spiro atoms typically have 4 neighbors (tetrahedral)
            if atom.GetDegree() == 4:
                # Check if this is truly a spiro (rings share only this atom)
                rings_with_atom = [r for r in atom_rings if atom_idx in r]
                if len(rings_with_atom) >= 2:
                    # Check if rings share only this one atom
                    ring1 = set(rings_with_atom[0])
                    ring2 = set(rings_with_atom[1])
                    shared = ring1 & ring2
                    if len(shared) == 1:  # Only the spiro atom is shared
                        return 1
    return 0


def _calc_sp3_fraction(mol: Chem.Mol) -> float:
    """Calculate fraction of sp³ hybridized carbons.

    Returns:
        Fraction of sp³ carbons (0-1)
    """
    return rdMolDescriptors.CalcFractionCSP3(mol)


def _calc_cyclic_sp3_fraction(mol: Chem.Mol) -> float:
    """Calculate fraction of cyclic sp³ carbons.

    Returns:
        Fraction of cyclic sp³ carbons relative to total carbons (0-1)
    """
    ring_info = mol.GetRingInfo()
    ring_atoms = set()
    for ring in ring_info.AtomRings():
        ring_atoms.update(ring)

    total_carbons = 0
    cyclic_sp3_carbons = 0

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            total_carbons += 1
            if atom.GetIdx() in ring_atoms:
                # Check if sp³ (no aromatic, degree suggests tetrahedral)
                if (
                    not atom.GetIsAromatic()
                    and atom.GetHybridization() == Chem.HybridizationType.SP3
                ):
                    cyclic_sp3_carbons += 1

    if total_carbons == 0:
        return 0.0
    return cyclic_sp3_carbons / total_carbons


def _calc_acyclic_sp3_fraction(mol: Chem.Mol) -> float:
    """Calculate fraction of acyclic sp³ carbons.

    Returns:
        Fraction of acyclic sp³ carbons relative to total carbons (0-1)
    """
    ring_info = mol.GetRingInfo()
    ring_atoms = set()
    for ring in ring_info.AtomRings():
        ring_atoms.update(ring)

    total_carbons = 0
    acyclic_sp3_carbons = 0

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            total_carbons += 1
            if atom.GetIdx() not in ring_atoms:
                # Check if sp³
                if atom.GetHybridization() == Chem.HybridizationType.SP3:
                    acyclic_sp3_carbons += 1

    if total_carbons == 0:
        return 0.0
    return acyclic_sp3_carbons / total_carbons


def _calc_q1_index(mol: Chem.Mol) -> float:
    """Calculate normalized quadratic index Q1 (Qindex).

    Q1 is calculated from the first Zagreb index (ZM1) as defined in
    Todeschini & Consonni's Handbook of Molecular Descriptors and
    implemented in DRAGON software.

    Formula: Q1 = 3 - 2*N + ZM1/2
    where ZM1 = sum of squared vertex degrees, N = number of atoms.

    Reference implementations:
    - PyBioMed (topology.py): CalculateQuadratic()
    - kotori-y/MCE-18 (mce18.py): CalculateQ1Index()

    Returns:
        Normalized quadratic index (typically 1-10 for drug-like molecules)
    """
    n_atoms = mol.GetNumAtoms()
    if n_atoms < 2:
        return 1.0

    # First Zagreb index (ZM1): sum of squared vertex degrees
    zm1 = sum(atom.GetDegree() ** 2 for atom in mol.GetAtoms())

    # Normalized quadratic index
    q1 = 3 - 2 * n_atoms + zm1 / 2.0

    return max(q1, 0.0)  # Ensure non-negative
