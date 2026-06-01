"""Descriptor key mappings and constants."""

from rdkit import RDLogger, rdBase

# Canonical mapping for descriptor keys (lowercase -> canonical case)
_DESCRIPTOR_KEY_MAP = {
    "logp": "logP",
    "molwt": "molWt",
    "tpsa": "tpsa",
    "hbd": "hbd",
    "hba": "hba",
    "qed": "qed",
    "fsp3": "fsp3",
    "mce18": "mce18",
    "sw": "sw",
    "n_atoms": "n_atoms",
    "n_heavy_atoms": "n_heavy_atoms",
    "n_rot_bonds": "n_rot_bonds",
    "n_rigid_bonds": "n_rigid_bonds",
    "n_rings": "n_rings",
    "fns_atoms": "fNS_atoms",
    "n_o_atoms": "n_O_atoms",
    "n_s_atoms": "n_S_atoms",
    "n_no_atoms": "n_NO_atoms",
    "n_small_rings_3_4": "n_small_rings_3_4",
    "max_acyclic_chain_length": "max_acyclic_chain_length",
    "has_spider_side_chains": "has_spider_side_chains",
    "fraction_ring_system": "fraction_ring_system",
    ".=o": ".=O",
    "c2r": "C2r",
    "c3r": "C3r",
    "car": "Car",
    "cs2": "Cs2",
    "cs3": "Cs3",
    "csp": "Csp",
    "nac": "Nac",
    "nd+": "Nd+",
    "nd0": "Nd0",
    "o_a": "O_a",
    "o_d": "O_d",
    "so2": "SO2",
    "sul": "Sul",
    "hal": "Hal",
}

TYPE_ALIAS_COLUMNS = (
    ".=O",
    "C2r",
    "C3r",
    "Car",
    "Cs2",
    "Cs3",
    "Csp",
    "Nac",
    "Nd+",
    "Nd0",
    "O_a",
    "O_d",
    "SO2",
    "Sul",
    "Hal",
)

STRUCTURAL_ELEMENT_LIMIT_MAP = {
    "n": "n_N_atoms",
    "N": "n_N_atoms",
    "n_n_atoms": "n_N_atoms",
    "n_N_atoms": "n_N_atoms",
    "o": "n_O_atoms",
    "O": "n_O_atoms",
    "n_o_atoms": "n_O_atoms",
    "n_O_atoms": "n_O_atoms",
    "s": "n_S_atoms",
    "S": "n_S_atoms",
    "n_s_atoms": "n_S_atoms",
    "n_S_atoms": "n_S_atoms",
}

# Disable RDKit warnings
RDLogger.DisableLog("rdApp.*")
rdBase.DisableLog("rdApp.*")
