from pathlib import Path

import pandas as pd
import pytest
import yaml

from hedgehog.struct_filters import main as structfilters_main


def _write_yaml(path: Path, payload: dict) -> None:
    path.write_text(yaml.safe_dump(payload), encoding="utf-8")


def test_single_stage_struct_filters_falls_back_to_mol_prep_when_descriptors_missing(
    tmp_path,
):
    descriptors_cfg = tmp_path / "config_descriptors.yml"
    _write_yaml(descriptors_cfg, {"run": True})

    mol_prep_output = tmp_path / "stages" / "00_mol_prep" / "filtered_molecules.csv"
    mol_prep_output.parent.mkdir(parents=True, exist_ok=True)
    mol_prep_output.write_text(
        "smiles,model_name,mol_idx\nCCO,m1,0\n", encoding="utf-8"
    )

    input_csv = tmp_path / "input.csv"
    input_csv.write_text("smiles,model_name\nCCO,m1\n", encoding="utf-8")

    config = {
        "config_descriptors": str(descriptors_cfg),
        "generated_mols_path": str(input_csv),
        "_run_single_stage_override": "struct_filters",
    }

    resolved = structfilters_main._get_input_path(
        config,
        "stages/03_structural_filters_post",
        tmp_path,
    )

    assert resolved == str(mol_prep_output)


def test_struct_filters_raises_when_descriptors_enabled_and_output_missing(tmp_path):
    descriptors_cfg = tmp_path / "config_descriptors.yml"
    _write_yaml(descriptors_cfg, {"run": True})

    input_csv = tmp_path / "input.csv"
    input_csv.write_text("smiles,model_name\nCCO,m1\n", encoding="utf-8")

    config = {
        "config_descriptors": str(descriptors_cfg),
        "generated_mols_path": str(input_csv),
    }

    with pytest.raises(FileNotFoundError):
        structfilters_main._get_input_path(
            config,
            "stages/03_structural_filters_post",
            tmp_path,
        )


def test_load_input_data_assigns_lp_mol_idx_when_missing(tmp_path):
    input_csv = tmp_path / "input.csv"
    pd.DataFrame({"smiles": ["CCO"], "model_name": ["m1"]}).to_csv(
        input_csv, index=False
    )

    input_df, model_name = structfilters_main._load_input_data(
        str(input_csv), run_base=tmp_path
    )

    assert model_name == "m1"
    assert input_df["mol_idx"].astype(str).str.startswith("LP-").all()


def test_load_input_data_keeps_existing_mol_idx_values(tmp_path):
    input_csv = tmp_path / "input.csv"
    pd.DataFrame(
        {
            "smiles": ["CCO", "CCC"],
            "model_name": ["m1", "m1"],
            "mol_idx": ["keep-1", "keep-2"],
        }
    ).to_csv(input_csv, index=False)

    input_df, _ = structfilters_main._load_input_data(str(input_csv), run_base=tmp_path)

    assert input_df["mol_idx"].tolist() == ["keep-1", "keep-2"]
