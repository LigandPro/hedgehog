"""Shared string constants used across the hedgehog pipeline.

Centralizes frequently repeated literals to reduce duplication (SonarQube S1192).
"""

# Column names used throughout the pipeline
COL_SMILES = "smiles"
COL_MODEL_NAME = "model_name"
COL_MOL_IDX = "mol_idx"
COL_PASS = "pass"

# Config section keys (master config references to sub-config YAML files)
CFG_STRUCT_FILTERS = "config_structFilters"
CFG_DESCRIPTORS = "config_descriptors"
CFG_DOCKING = "config_docking"

# Common config keys
KEY_FOLDER_TO_SAVE = "folder_to_save"

# Docking tool names
TOOL_GNINA = "gnina"
TOOL_SMINA = "smina"
