from pathlib import Path

import pandas as pd
import yaml
from rdkit import Chem

from hedgehog.config_alignment import finalize_global_alignment
from tests.test_config_alignment import _base_master


def test_global_alignment_protects_one_shared_target_cohort(tmp_path):
    master = _base_master(tmp_path)
    master["alignment"]["descriptor_bounds_mode"] = "expand"
    master["_continue_mode"] = True
    master["_continue_completed_stages"] = ["mol_prep"]
    master["_run_stage_selection_override"] = [
        "mol_prep",
        "descriptors",
        "struct_filters",
        "docking",
    ]

    targets = tmp_path / "targets.csv"
    target_run = tmp_path / "target_run"
    target_ids = [f"mol-{index}" for index in range(10)]
    target_smiles = ["CCO"] * 10
    targets.write_text(
        "smiles,mol_idx\n"
        + "".join(
            f"{smiles},{mol_idx}\n"
            for smiles, mol_idx in zip(target_smiles, target_ids, strict=True)
        ),
        encoding="utf-8",
    )
    sampled = target_run / "input" / "sampled_molecules.csv"
    sampled.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "smiles": target_smiles,
            "model_name": ["KRAS"] * 10,
            "mol_idx": target_ids,
        }
    ).to_csv(sampled, index=False)

    descriptor_dir = target_run / "stages" / "02_descriptors_initial" / "metrics"
    descriptor_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "smiles": target_smiles,
            "mol_idx": target_ids,
            "metric": [0, 100, *range(2, 10)],
            "Hal": [0] * 10,
            "n_N_atoms": [0] * 10,
            "n_NO_atoms": [1] * 10,
        }
    ).to_csv(descriptor_dir / "descriptors_all.csv", index=False)

    structural_dir = (
        target_run / "stages" / "03_structural_filters_post" / "common_alerts"
    )
    structural_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "smiles": target_smiles,
            "model_name": ["KRAS"] * 10,
            "mol_idx": target_ids,
            "pass_Rule A": [False, True, *([True] * 8)],
            "pass_Rule B": [True, False, *([True] * 8)],
            "pass_any": [True] * 10,
        }
    ).to_csv(structural_dir / "extended.csv", index=False)

    docking_dir = target_run / "stages" / "05_docking"
    docking_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "smiles": target_smiles,
            "model_name": ["KRAS"] * 10,
            "mol_idx": target_ids,
        }
    ).to_csv(docking_dir / "input_molecules.csv", index=False)
    writer = Chem.SDWriter(str(docking_dir / "docking_out.sdf"))
    try:
        for index, mol_idx in enumerate(target_ids):
            for tool, offset in {"smina": -12, "gnina": -11, "matcha": -10}.items():
                mol = Chem.MolFromSmiles("CCO")
                mol.SetProp("mol_idx", mol_idx)
                mol.SetProp("docking_tool", tool)
                mol.SetDoubleProp("minimizedAffinity", offset + index)
                writer.write(mol)
    finally:
        writer.close()

    aligned, _master_path, audit_path = finalize_global_alignment(
        master,
        target_run,
        tmp_path / "alignment",
        str(targets),
        80,
    )
    for runtime_key in (
        "_continue_mode",
        "_continue_completed_stages",
        "_run_stage_selection_override",
        "_run_single_stage_override",
    ):
        assert runtime_key not in aligned

    audit = yaml.safe_load(audit_path.read_text())
    assert audit["selection_method"] == "global_protected_target_cohort"
    assert audit["descriptor_bounds_mode"] == "expand"
    assert audit["stages"]["descriptors"]["bounds_mode"] == "expand"
    assert audit["global_guarantee"]["status"] == "verified"
    assert audit["global_guarantee"]["required_retained_molecules"] == 8
    for stage in ("descriptors", "docking"):
        assert audit["stages"][stage]["protected_retained_molecules"] == 8
    assert "synthesis" not in audit["stages"]
    assert "synthesis_bounds_mode" not in audit
    assert audit["stages"]["struct_filters"]["status"] == "source_config_preserved"
    assert "protected_retained_molecules" not in audit["stages"]["struct_filters"]
    assert (
        Path(aligned["config_structFilters"]).read_text()
        == Path(master["config_structFilters"]).read_text()
    )
    descriptors = yaml.safe_load(Path(aligned["config_descriptors"]).read_text())
    assert descriptors["borders"]["metric_min"] == -10
    assert descriptors["borders"]["metric_max"] == 200
    synthesis = yaml.safe_load(Path(aligned["config_synthesis"]).read_text())
    docking_filters = yaml.safe_load(
        Path(aligned["config_docking_filters"]).read_text()
    )
    assert synthesis["filter_solved_only"] is True
    assert (
        Path(aligned["config_synthesis"]).read_text()
        == Path(master["config_synthesis"]).read_text()
    )
    assert "synthesis_bounds_mode" not in aligned.get("alignment", {})
    source_docking_filters = yaml.safe_load(
        Path(master["config_docking_filters"]).read_text()
    )
    assert docking_filters == source_docking_filters
    assert (
        Path(aligned["config_docking_filters"]).read_text()
        == Path(master["config_docking_filters"]).read_text()
    )
    assert audit["stages"]["docking_filters"]["status"] == ("source_config_preserved")
    protected = pd.read_csv(audit["global_guarantee"]["protected_cohort_path"])
    assert len(protected) == 8


def test_global_alignment_honors_descriptor_only_stage_selection(tmp_path):
    master = _base_master(tmp_path)
    master["_run_stage_selection_override"] = ["mol_prep", "descriptors"]

    targets = tmp_path / "targets.csv"
    target_run = tmp_path / "target_run"
    target_ids = [f"mol-{index}" for index in range(10)]
    target_smiles = ["CCO"] * 10
    targets.write_text(
        "smiles,mol_idx\n"
        + "".join(
            f"{smiles},{mol_idx}\n"
            for smiles, mol_idx in zip(target_smiles, target_ids, strict=True)
        ),
        encoding="utf-8",
    )
    sampled = target_run / "input" / "sampled_molecules.csv"
    sampled.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "smiles": target_smiles,
            "model_name": ["target"] * 10,
            "mol_idx": target_ids,
        }
    ).to_csv(sampled, index=False)

    descriptor_dir = target_run / "stages" / "02_descriptors_initial" / "metrics"
    descriptor_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "smiles": target_smiles,
            "mol_idx": target_ids,
            "metric": range(10),
            "Hal": [0] * 10,
            "n_N_atoms": [0] * 10,
            "n_NO_atoms": [1] * 10,
        }
    ).to_csv(descriptor_dir / "descriptors_all.csv", index=False)

    _aligned, _master_path, audit_path = finalize_global_alignment(
        master,
        target_run,
        tmp_path / "alignment",
        str(targets),
        80,
    )

    audit = yaml.safe_load(audit_path.read_text())
    assert audit["global_guarantee"]["required_retained_molecules"] == 8
    assert audit["global_guarantee"]["verified_stages"] == [
        "mol_prep",
        "descriptors",
    ]
    assert audit["stages"]["descriptors"]["protected_retained_molecules"] == 8
    assert audit["stages"]["struct_filters"]["status"] == "pending_metrics"
