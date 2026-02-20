"""Shared constants for the HEDGEHOG pipeline."""

# Directory names
DIR_INPUT = "input"
DIR_STAGES = "stages"
DIR_OUTPUT = "output"
DIR_CONFIGS = "configs"

# Stage subdirectories
DIR_MOL_PREP = "stages/00_mol_prep"
DIR_DESCRIPTORS_INITIAL = "stages/01_descriptors_initial"
DIR_STRUCT_FILTERS_POST = "stages/03_structural_filters_post"
DIR_SYNTHESIS = "stages/04_synthesis"
DIR_DOCKING = "stages/05_docking"
DIR_DOCKING_FILTERS = "stages/06_docking_filters"
DIR_DESCRIPTORS_FINAL = "stages/07_descriptors_final"

# Legacy names for backwards compatibility
DIR_DESCRIPTORS = DIR_DESCRIPTORS_INITIAL
DIR_STRUCT_FILTERS = DIR_STRUCT_FILTERS_POST
DIR_FINAL_DESCRIPTORS = DIR_DESCRIPTORS_FINAL
DIR_RUN_CONFIGS = DIR_CONFIGS

# File names
FILE_SAMPLED_MOLECULES = "sampled_molecules.csv"
FILE_FINAL_MOLECULES = "final_molecules.csv"
FILE_FILTERED_MOLECULES = "filtered_molecules.csv"
FILE_PASS_SMILES_TEMPLATE = "filtered_molecules.csv"
FILE_MASTER_CONFIG = "master_config_resolved.yml"
FILE_GNINA_OUTPUT = "gnina_out.sdf"
FILE_SMINA_OUTPUT = "smina_out.sdf"

DOCKING_SCORE_COLUMNS = [
    "gnina_affinity",
    "gnina_cnnscore",
    "gnina_cnnaffinity",
    "gnina_cnn_vs",
    "smina_affinity",
]

# Stage names
STAGE_MOL_PREP = "mol_prep"
STAGE_DESCRIPTORS = "descriptors"
STAGE_STRUCT_FILTERS = "struct_filters"
STAGE_SYNTHESIS = "synthesis"
STAGE_DOCKING = "docking"
STAGE_DOCKING_FILTERS = "docking_filters"
STAGE_FINAL_DESCRIPTORS = "final_descriptors"

# Config keys
CONFIG_MOL_PREP = "config_mol_prep"
CONFIG_DESCRIPTORS = "config_descriptors"
CONFIG_STRUCT_FILTERS = "config_structFilters"
CONFIG_SYNTHESIS = "config_synthesis"
CONFIG_DOCKING = "config_docking"
CONFIG_DOCKING_FILTERS = "config_docking_filters"
CONFIG_RUN_KEY = "run"
CONFIG_TOOLS = "tools"
CONFIG_FOLDER_TO_SAVE = "folder_to_save"

# Command-line override keys
OVERRIDE_SINGLE_STAGE = "_run_single_stage_override"

# Docking tools
DOCKING_TOOL_SMINA = "smina"
DOCKING_TOOL_GNINA = "gnina"
DOCKING_TOOL_BOTH = "both"
DOCKING_RESULTS_DIR_TEMPLATE = {
    DOCKING_TOOL_SMINA: f"{DIR_DOCKING}/smina",
    DOCKING_TOOL_GNINA: f"{DIR_DOCKING}/gnina",
}

FILE_RUN_INCOMPLETE = ".RUN_INCOMPLETE"

# Stage directory mapping (used by report generator)
STAGE_DIRS = {
    "mol_prep": DIR_MOL_PREP,
    "descriptors_initial": DIR_DESCRIPTORS_INITIAL,
    "struct_filters_post": DIR_STRUCT_FILTERS_POST,
    "synthesis": DIR_SYNTHESIS,
    "docking": DIR_DOCKING,
    "docking_filters": DIR_DOCKING_FILTERS,
    "descriptors_final": DIR_DESCRIPTORS_FINAL,
}

# Stage display names (used by report generator)
STAGE_DISPLAY_NAMES = {
    "mol_prep": "Mol Prep",
    "descriptors_initial": "Initial Descriptors",
    "struct_filters_post": "Post-Descriptors Filters",
    "synthesis": "Synthesis Analysis",
    "docking": "Molecular Docking",
    "docking_filters": "Docking Filters",
    "descriptors_final": "Final Descriptors",
}
