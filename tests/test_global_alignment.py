from pathlib import Path

import pandas as pd
import yaml
from rdkit import Chem

from hedgehog.config_alignment import finalize_global_alignment
from tests.test_config_alignment import _base_master


def test_global_alignment_protects_one_shared_target_cohort(tmp_path):
    master = _base_master(tmp_path)
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

    synthesis_dir = target_run / "stages" / "04_synthesis"
    synthesis_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "smiles": target_smiles,
            "mol_idx": target_ids,
            "sa_score": [100, 0, *range(2, 10)],
        }
    ).to_csv(synthesis_dir / "synthesis_scores.csv", index=False)

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

    audit = yaml.safe_load(audit_path.read_text())
    assert audit["selection_method"] == "global_protected_target_cohort"
    assert audit["global_guarantee"]["status"] == "verified"
    assert audit["global_guarantee"]["required_retained_molecules"] == 8
    for stage in ("descriptors", "struct_filters", "synthesis", "docking"):
        assert audit["stages"][stage]["protected_retained_molecules"] == 8
    synthesis = yaml.safe_load(Path(aligned["config_synthesis"]).read_text())
    docking_filters = yaml.safe_load(
        Path(aligned["config_docking_filters"]).read_text()
    )
    assert synthesis["filter_solved_only"] is False
    assert docking_filters["run"] is False
    protected = pd.read_csv(audit["global_guarantee"]["protected_cohort_path"])
    assert len(protected) == 8
