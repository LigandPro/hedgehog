"""Tests for MCE-18 (Molecular Complexity Estimation) descriptor."""

import pytest
from rdkit import Chem

from hedgehog.utils.mce18 import (
    _calc_acyclic_sp3_fraction,
    _calc_cyclic_sp3_fraction,
    _calc_q1_index,
    _calc_sp3_fraction,
    _has_aliphatic_ring,
    _has_aromatic_ring,
    _has_chiral_center,
    _has_spiro_center,
    compute_mce18,
)


class TestMCE18Components:
    """Tests for individual MCE-18 components."""

    def test_aromatic_ring_benzene(self):
        """Benzene should have aromatic ring."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        assert _has_aromatic_ring(mol) == 1

    def test_aromatic_ring_cyclohexane(self):
        """Cyclohexane should NOT have aromatic ring."""
        mol = Chem.MolFromSmiles("C1CCCCC1")
        assert _has_aromatic_ring(mol) == 0

    def test_aromatic_ring_pyridine(self):
        """Pyridine (heteroaromatic) should have aromatic ring."""
        mol = Chem.MolFromSmiles("c1ccncc1")
        assert _has_aromatic_ring(mol) == 1

    def test_aliphatic_ring_cyclohexane(self):
        """Cyclohexane should have aliphatic ring."""
        mol = Chem.MolFromSmiles("C1CCCCC1")
        assert _has_aliphatic_ring(mol) == 1

    def test_aliphatic_ring_benzene(self):
        """Benzene should NOT have aliphatic ring."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        assert _has_aliphatic_ring(mol) == 0

    def test_aliphatic_ring_tetrahydrofuran(self):
        """Tetrahydrofuran should have aliphatic ring."""
        mol = Chem.MolFromSmiles("C1CCOC1")
        assert _has_aliphatic_ring(mol) == 1

    def test_chiral_center_present(self):
        """Molecule with chiral center should be detected."""
        mol = Chem.MolFromSmiles("C[C@H](O)CC")
        assert _has_chiral_center(mol) == 1

    def test_chiral_center_absent(self):
        """Symmetric molecule should NOT have chiral center."""
        mol = Chem.MolFromSmiles("CC(C)C")
        assert _has_chiral_center(mol) == 0

    def test_spiro_center_present(self):
        """Spiro[5.5]undecane should have spiro center."""
        mol = Chem.MolFromSmiles("C1CCC2(CC1)CCCCC2")
        assert _has_spiro_center(mol) == 1

    def test_spiro_center_absent(self):
        """Fused ring (decalin) should NOT have spiro center."""
        mol = Chem.MolFromSmiles("C1CCC2CCCCC2C1")
        assert _has_spiro_center(mol) == 0

    def test_spiro_center_single_ring(self):
        """Single ring molecule cannot have spiro center."""
        mol = Chem.MolFromSmiles("C1CCCCC1")
        assert _has_spiro_center(mol) == 0

    def test_sp3_fraction_ethane(self):
        """Ethane should have 100% sp3 carbons."""
        mol = Chem.MolFromSmiles("CC")
        assert _calc_sp3_fraction(mol) == 1.0

    def test_sp3_fraction_benzene(self):
        """Benzene should have 0% sp3 carbons."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        assert _calc_sp3_fraction(mol) == 0.0

    def test_cyclic_sp3_fraction_cyclohexane(self):
        """Cyclohexane: all carbons are cyclic sp3."""
        mol = Chem.MolFromSmiles("C1CCCCC1")
        assert _calc_cyclic_sp3_fraction(mol) == 1.0

    def test_cyclic_sp3_fraction_ethane(self):
        """Ethane: no cyclic carbons."""
        mol = Chem.MolFromSmiles("CC")
        assert _calc_cyclic_sp3_fraction(mol) == 0.0

    def test_acyclic_sp3_fraction_ethane(self):
        """Ethane: all carbons are acyclic sp3."""
        mol = Chem.MolFromSmiles("CC")
        assert _calc_acyclic_sp3_fraction(mol) == 1.0

    def test_acyclic_sp3_fraction_cyclohexane(self):
        """Cyclohexane: no acyclic carbons."""
        mol = Chem.MolFromSmiles("C1CCCCC1")
        assert _calc_acyclic_sp3_fraction(mol) == 0.0

    def test_q1_index_positive(self):
        """Q1 index should be positive for any molecule."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        q1 = _calc_q1_index(mol)
        assert q1 > 0

    def test_q1_index_single_atom(self):
        """Q1 index for single atom should be 1.0."""
        mol = Chem.MolFromSmiles("C")
        q1 = _calc_q1_index(mol)
        assert q1 == 1.0


class TestMCE18Values:
    """Tests for complete MCE-18 calculation."""

    def test_benzene_low_complexity(self):
        """Benzene should have low MCE-18 (mostly aromatic, no sp3)."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        mce18 = compute_mce18(mol)
        assert mce18 is not None
        assert mce18 < 10  # Simple aromatic, low complexity

    def test_cyclohexane_has_complexity(self):
        """Cyclohexane has sp3 carbons and aliphatic ring."""
        mol = Chem.MolFromSmiles("C1CCCCC1")
        mce18 = compute_mce18(mol)
        assert mce18 is not None
        assert mce18 > 2  # Has sp3 and NAR components

    def test_spiro_molecule_higher_complexity(self):
        """Spiro molecule should have higher complexity."""
        mol = Chem.MolFromSmiles("C1CCC2(CC1)CCCCC2")
        mce18 = compute_mce18(mol)
        assert mce18 is not None
        assert mce18 > 5  # Spiro center adds complexity

    def test_drug_like_molecule(self):
        """Typical drug-like molecule should have moderate MCE-18."""
        mol = Chem.MolFromSmiles("CC(C)Cc1ccc(C(C)C(=O)O)cc1")  # Ibuprofen
        mce18 = compute_mce18(mol)
        assert mce18 is not None
        assert mce18 > 0

    def test_invalid_molecule_returns_none(self):
        """Invalid molecule should return None."""
        assert compute_mce18(None) is None

    def test_mce18_is_deterministic(self):
        """Same molecule should give same MCE-18."""
        mol = Chem.MolFromSmiles("CCO")
        val1 = compute_mce18(mol)
        val2 = compute_mce18(mol)
        assert val1 == val2

    def test_mce18_returns_float(self):
        """MCE-18 should return a float value."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        mce18 = compute_mce18(mol)
        assert isinstance(mce18, float)

    def test_mce18_rounded_to_4_decimals(self):
        """MCE-18 should be rounded to 4 decimal places."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        mce18 = compute_mce18(mol)
        # Check that result has at most 4 decimal places
        assert mce18 == round(mce18, 4)


class TestMCE18Integration:
    """Integration tests with descriptors module."""

    def test_mce18_in_descriptors(self):
        """MCE-18 should be computed alongside other descriptors."""
        from hedgehog.stages.descriptors.utils import (
            _compute_single_molecule_descriptors,
        )

        mol = Chem.MolFromSmiles("c1ccccc1")
        result = _compute_single_molecule_descriptors(mol, "test", "idx-0")

        assert "mce18" in result
        assert isinstance(result["mce18"], (int, float, type(None)))

    def test_mce18_in_descriptor_key_map(self):
        """MCE-18 should be in the descriptor key map."""
        from hedgehog.stages.descriptors.utils import _DESCRIPTOR_KEY_MAP

        assert "mce18" in _DESCRIPTOR_KEY_MAP
        assert _DESCRIPTOR_KEY_MAP["mce18"] == "mce18"

    def test_mce18_in_moleval_config_map(self):
        """MCE-18 should be in the MolEval config map."""
        from hedgehog.reporting.moleval_metrics import _METRIC_CONFIG_MAP

        assert "MCE18" in _METRIC_CONFIG_MAP
        assert _METRIC_CONFIG_MAP["MCE18"] == "mce18"
