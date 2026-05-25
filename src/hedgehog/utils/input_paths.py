"""Centralized input path resolution for pipeline stages."""

from collections.abc import Sequence
from pathlib import Path

# Stage directory basenames (under stages/)
_STAGE_SYNTHESIS = "04_synthesis"
_STAGE_STRUCT_FILTERS_POST = "03_structural_filters_post"
_STAGE_DESCRIPTORS_INITIAL = "02_descriptors_initial"
_STAGE_MOL_PREP = "01_mol_prep"

# Priority order for input sources (hierarchical structure)
INPUT_SOURCE_PRIORITY_NEW = [
    ("stages", _STAGE_SYNTHESIS, "filtered_molecules.csv"),
    ("stages", _STAGE_STRUCT_FILTERS_POST, "filtered_molecules.csv"),
    ("stages", _STAGE_DESCRIPTORS_INITIAL, "filtered", "filtered_molecules.csv"),
    ("stages", _STAGE_MOL_PREP, "filtered_molecules.csv"),
    ("input", "sampled_molecules.csv"),
]

# Flat layout (pre-stages/ hierarchy)
INPUT_SOURCE_PRIORITY_LEGACY = [
    ("Synthesis", "passSynthesisSMILES.csv"),
    ("StructFilters", "passStructFiltersSMILES.csv"),
    ("Descriptors", "passDescriptorsSMILES.csv"),
    ("sampled_molecules.csv",),
]

STAGE_DIRECTORIES = {
    "synthesis": ["stages/04_synthesis", "Synthesis"],
    "struct_filters": [
        "stages/03_structural_filters_post",
        "StructFilters",
    ],
    "descriptors": ["stages/02_descriptors_initial", "Descriptors"],
    "docking": ["stages/05_docking"],
    "mol_prep": ["stages/01_mol_prep"],
}


def get_all_input_candidates(base_folder: Path) -> list[Path]:
    """Get all potential input file paths in priority order."""
    base = Path(base_folder)
    candidates = []

    for parts in INPUT_SOURCE_PRIORITY_NEW:
        candidates.append(base.joinpath(*parts))

    for parts in INPUT_SOURCE_PRIORITY_LEGACY:
        candidates.append(base.joinpath(*parts))

    return candidates


def find_latest_input_source(
    base_folder: Path,
    skip_stages: Sequence[str] | None = None,
) -> Path | None:
    """Find the most recent input source file."""
    base = Path(base_folder)
    skip_stages = set(skip_stages or [])

    skip_dirs = set()
    for stage in skip_stages:
        for dir_path in STAGE_DIRECTORIES.get(stage, []):
            skip_dirs.add(dir_path.lower())

    candidates = get_all_input_candidates(base)

    for candidate in candidates:
        rel_path = str(candidate.relative_to(base)).lower()
        should_skip = any(skip_dir in rel_path for skip_dir in skip_dirs)

        if should_skip:
            continue

        if _file_exists_and_not_empty(candidate):
            return candidate

    return None


def find_sampled_molecules(base_folder: Path) -> Path | None:
    """Find the sampled molecules file."""
    base = Path(base_folder)
    candidates = [
        base / "input" / "sampled_molecules.csv",
        base / "sampled_molecules.csv",
    ]

    for candidate in candidates:
        if _file_exists_and_not_empty(candidate):
            return candidate

    return None


def _file_exists_and_not_empty(file_path: Path) -> bool:
    """Check if a file exists and is not empty."""
    try:
        return file_path.exists() and file_path.stat().st_size > 0
    except OSError:
        return False
