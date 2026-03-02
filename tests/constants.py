"""Shared test constants to reduce string duplication (SonarQube S1192)."""

# Test SMILES strings
SMILES_ETHANOL = "CCO"
SMILES_BENZENE = "c1ccccc1"
SMILES_ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"

# Common DataFrame column names
COL_SMILES = "smiles"
COL_MODEL_NAME = "model_name"
COL_MOL_IDX = "mol_idx"

# Common test model name
MODEL_TEST = "test"

# Common file names used in tests
FILE_FILTERED_MOLECULES = "filtered_molecules.csv"

# Config section keys
CFG_STRUCT_FILTERS = "config_structFilters"
CFG_DESCRIPTORS = "config_descriptors"
