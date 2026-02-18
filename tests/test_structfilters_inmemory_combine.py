from pathlib import Path

import pandas as pd
import yaml

from hedgehog.stages.structFilters.utils import (
    combine_filter_results_in_memory,
    filter_data,
)


def _write_yaml(path: Path, payload: dict) -> None:
    path.write_text(yaml.safe_dump(payload), encoding="utf-8")


def _make_pass_mask(input_df: pd.DataFrame, passed_indices: set[int]) -> pd.DataFrame:
    mask = input_df[["smiles", "model_name", "mol_idx"]].copy()
    mask["pass"] = mask["mol_idx"].apply(lambda idx: idx in passed_indices)
    return mask


def test_inmemory_combine_matches_legacy_filtered_output(tmp_path):
    all_input = pd.DataFrame(
        {
            "smiles": ["CCO", "CCN", "CCC", "c1ccccc1"],
            "model_name": ["m1", "m1", "m1", "m1"],
            "mol_idx": [0, 1, 2, 3],
        }
    )

    legacy_root = tmp_path / "legacy"
    legacy_stage_dir = "stages/struct_filters_test"
    legacy_stage_path = legacy_root / legacy_stage_dir
    legacy_stage_path.mkdir(parents=True, exist_ok=True)

    sampled = legacy_root / "sampled_molecules.csv"
    all_input.to_csv(sampled, index=False)

    filter1 = all_input[all_input["mol_idx"].isin({0, 1, 2})].copy()
    filter2 = all_input[all_input["mol_idx"].isin({1, 2, 3})].copy()

    (legacy_stage_path / "filter_a").mkdir(parents=True, exist_ok=True)
    (legacy_stage_path / "filter_b").mkdir(parents=True, exist_ok=True)
    filter1.to_csv(legacy_stage_path / "filter_a" / "filtered_molecules.csv", index=False)
    filter2.to_csv(legacy_stage_path / "filter_b" / "filtered_molecules.csv", index=False)

    struct_cfg = legacy_root / "config_structFilters.yml"
    _write_yaml(struct_cfg, {"calculate_filter_a": True, "calculate_filter_b": True})
    legacy_config = {
        "folder_to_save": str(legacy_root),
        "generated_mols_path": str(sampled),
        "config_structFilters": str(struct_cfg),
    }
    filter_data(legacy_config, legacy_stage_dir)
    legacy_filtered = pd.read_csv(legacy_stage_path / "filtered_molecules.csv")

    modern_root = tmp_path / "modern"
    modern_stage_path = modern_root / "stages" / "struct_filters_test"
    modern_stage_path.mkdir(parents=True, exist_ok=True)

    pass_masks = {
        "filter_a": _make_pass_mask(all_input, {0, 1, 2}),
        "filter_b": _make_pass_mask(all_input, {1, 2, 3}),
    }
    modern_filtered, modern_failed = combine_filter_results_in_memory(
        modern_stage_path, all_input, pass_masks
    )
    modern_filtered = pd.read_csv(modern_stage_path / "filtered_molecules.csv")

    expected_legacy = legacy_filtered.sort_values("mol_idx").reset_index(drop=True)
    expected_modern = modern_filtered.sort_values("mol_idx").reset_index(drop=True)
    pd.testing.assert_frame_equal(expected_modern, expected_legacy)

    assert "pass_filter_a" in modern_failed.columns
    assert "pass_filter_b" in modern_failed.columns
