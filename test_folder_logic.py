#!/usr/bin/env python3
"""
Test script for hybrid folder selection logic
Tests all scenarios with --reuse, --force-new, --stage, and --mols combinations
"""
from pathlib import Path
import re


def _get_unique_results_folder(base_folder):
    """Generate a unique folder name by appending a number if the folder already exists."""
    base_folder = Path(base_folder)
    if not base_folder.exists():
        return base_folder
    if base_folder.exists() and any(base_folder.iterdir()):
        base_name = base_folder.name
        parent = base_folder.parent
        match = re.match(r'^(.+?)(\d+)$', base_name)
        if match:
            name_without_counter = match.group(1)
            start_counter = int(match.group(2)) + 1
        else:
            name_without_counter = base_name
            start_counter = 1
        counter = start_counter
        while True:
            new_folder = parent / f"{name_without_counter}{counter}"
            if not new_folder.exists() or not any(new_folder.iterdir()):
                return new_folder
            counter += 1
    return base_folder


def simulate_folder_selection(original_folder, reuse_folder=False, force_new_folder=False,
                               stage=None, generated_mols_path=None):
    """
    Simulate the hybrid folder selection logic.

    Returns tuple: (chosen_folder, description)
    """
    original_folder = Path(original_folder)

    # Validate conflicting flags
    if reuse_folder and force_new_folder:
        return None, "ERROR: Cannot use --reuse and --force-new together"

    if reuse_folder:
        # Explicit reuse
        return original_folder, "Reusing folder (--reuse flag)"
    elif force_new_folder:
        # Explicit new
        folder = _get_unique_results_folder(original_folder)
        if folder != original_folder:
            return folder, f"Creating new folder (--force-new flag)"
        return folder, "Using folder (--force-new but folder is empty)"
    else:
        # Automatic logic (Hybrid Variant 4)
        if stage and not generated_mols_path:
            # Stage rerun on existing data → reuse
            return original_folder, "Reusing folder for stage execution (auto)"
        else:
            # Full run OR stage with new molecules → create new
            folder = _get_unique_results_folder(original_folder)
            if folder != original_folder:
                return folder, "Creating new folder (auto)"
            return folder, "Using folder (empty or doesn't exist)"


def test_all_scenarios():
    """Test all usage scenarios."""
    print("=" * 70)
    print("Testing Hybrid Folder Selection Logic (Variant 4 + Explicit Flags)")
    print("=" * 70)
    print()

    scenarios = [
        # (name, original_folder, reuse, force_new, stage, mols)
        ("1. Full run (no existing results)",
         "results/test", False, False, None, None),

        ("2. Full run (existing results → auto-increment)",
         "results/test", False, False, None, None),

        ("3. Stage rerun without --mols (auto-reuse)",
         "results/test", False, False, "docking", None),

        ("4. Stage rerun with --mols (auto-increment)",
         "results/test", False, False, "descriptors", "new_data.csv"),

        ("5. Explicit --reuse flag",
         "results/test", True, False, None, None),

        ("6. Explicit --force-new flag",
         "results/test", False, True, None, None),

        ("7. Stage with --reuse (override auto-reuse)",
         "results/test", True, False, "docking", None),

        ("8. Stage with --force-new (override auto-reuse)",
         "results/test", False, True, "docking", None),

        ("9. Stage with --mols and --reuse",
         "results/test", True, False, "descriptors", "new.csv"),

        ("10. Conflicting flags (--reuse + --force-new)",
         "results/test", True, True, None, None),
    ]

    print("Scenarios (assuming results/test exists and has content):")
    print("-" * 70)

    for name, folder, reuse, force_new, stage, mols in scenarios:
        result_folder, description = simulate_folder_selection(
            folder, reuse, force_new, stage, mols
        )

        # Build command for clarity
        cmd_parts = ["uv run hedge run"]
        if stage:
            cmd_parts.append(f"--stage {stage}")
        if mols:
            cmd_parts.append(f"--mols {mols}")
        if reuse:
            cmd_parts.append("--reuse")
        if force_new:
            cmd_parts.append("--force-new")

        cmd = " ".join(cmd_parts)

        print(f"\n{name}")
        print(f"  Command: {cmd}")
        if result_folder:
            print(f"  Result:  {result_folder}")
        print(f"  Reason:  {description}")

    print("\n" + "=" * 70)
    print("Summary of Logic:")
    print("=" * 70)
    print("""
Priority:
1. --reuse         → Always reuse existing folder
2. --force-new     → Always create new incremented folder
3. Auto (no flags):
   - --stage WITHOUT --mols → Reuse (typical stage rerun)
   - --stage WITH --mols    → New folder (new data)
   - Full run (no --stage)  → New folder (protect existing)

This ensures:
✓ Stage reruns work correctly (find existing artifacts)
✓ New data gets fresh folders (no mixing)
✓ Users have full control via flags when needed
✓ No accidental overwrites
    """)


if __name__ == "__main__":
    test_all_scenarios()
