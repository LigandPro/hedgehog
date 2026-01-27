"""Centralized input path resolution for pipeline stages.

This module provides consistent input file discovery across all pipeline stages,
supporting both new hierarchical structure and legacy flat structure.
"""

from collections.abc import Sequence
from pathlib import Path

# Priority order for input sources (new hierarchical structure)
INPUT_SOURCE_PRIORITY_NEW = [
    ("stages", "04_synthesis", "filtered_molecules.csv"),
    ("stages", "03_structural_filters_post", "filtered_molecules.csv"),
    ("stages", "01_descriptors_initial", "filtered", "filtered_molecules.csv"),
    ("stages", "02_structural_filters_pre", "filtered_molecules.csv"),
    ("input", "sampled_molecules.csv"),
]

# Priority order for input sources (legacy flat structure)
INPUT_SOURCE_PRIORITY_LEGACY = [
    ("Synthesis", "passSynthesisSMILES.csv"),
    ("StructFilters", "passStructFiltersSMILES.csv"),
    ("Descriptors", "passDescriptorsSMILES.csv"),
    ("sampled_molecules.csv",),
]

# Stage name to directory mapping for skip logic
STAGE_DIRECTORIES = {
    "synthesis": ["stages/04_synthesis", "Synthesis"],
    "struct_filters": [
        "stages/03_structural_filters_post",
        "stages/02_structural_filters_pre",
        "StructFilters",
    ],
    "descriptors": ["stages/01_descriptors_initial", "Descriptors"],
    "docking": ["stages/05_docking"],
}


def get_all_input_candidates(base_folder: Path) -> list[Path]:
    """Get all potential input file paths in priority order.

    Args:
        base_folder: Base results folder path

    Returns:
        List of Path objects for all candidate input files (new + legacy)
    """
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
    """Find the most recent input source file.

    Searches through predefined priority order of stage outputs to find
    the first existing non-empty file.

    Args:
        base_folder: Base results folder path
        skip_stages: Optional list of stage names to skip (e.g., ['descriptors'])

    Returns:
        Path to the found input file, or None if not found
    """
    base = Path(base_folder)
    skip_stages = set(skip_stages or [])

    # Build list of directories to skip
    skip_dirs = set()
    for stage in skip_stages:
        for dir_path in STAGE_DIRECTORIES.get(stage, []):
            skip_dirs.add(dir_path.lower())

    candidates = get_all_input_candidates(base)

    for candidate in candidates:
        # Check if this path should be skipped
        rel_path = str(candidate.relative_to(base)).lower()
        should_skip = any(skip_dir in rel_path for skip_dir in skip_dirs)

        if should_skip:
            continue

        if _file_exists_and_not_empty(candidate):
            return candidate

    return None


def find_sampled_molecules(base_folder: Path) -> Path | None:
    """Find the sampled molecules file.

    Checks both new and legacy locations.

    Args:
        base_folder: Base results folder path

    Returns:
        Path to sampled_molecules.csv or None if not found
    """
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
    """Check if a file exists and is not empty.

    Args:
        file_path: Path to check

    Returns:
        True if file exists and has size > 0
    """
    try:
        return file_path.exists() and file_path.stat().st_size > 0
    except OSError:
        return False
