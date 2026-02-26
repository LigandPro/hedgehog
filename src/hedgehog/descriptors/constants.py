"""Descriptor key mappings and constants."""

from rdkit import RDLogger, rdBase

# Canonical mapping for descriptor keys (lowercase -> canonical case)
# Used to prevent confusion between similar names like logP and clogP
_DESCRIPTOR_KEY_MAP = {
    "logp": "logP",
    "clogp": "clogP",
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
}

# Disable RDKit warnings
RDLogger.DisableLog("rdApp.*")
rdBase.DisableLog("rdApp.*")
