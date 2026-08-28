"""Tests for the LigandPro/docking screening adapter."""

import gzip
import json

from rdkit import Chem
from rdkit.Chem import AllChem

from hedgehog.docking.docking_repository import (
    _prepare_screening_dataset,
    _write_matcha_outputs,
)


def _molecule(name: str) -> Chem.Mol:
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(mol, randomSeed=7)
    mol.SetProp("_Name", name)
    return mol


def test_prepare_screening_dataset_creates_multiligand_layout(tmp_path):
    """Input SDF records should become one target's multiligand dataset."""
    receptor = tmp_path / "receptor.pdb"
    receptor.write_text("END\n", encoding="utf-8")
    ligands = tmp_path / "ligands.sdf"
    writer = Chem.SDWriter(str(ligands))
    writer.write(_molecule("ligand one"))
    writer.write(_molecule("ligand/one"))
    writer.close()
    dataset = tmp_path / "dataset"

    ligand_ids = _prepare_screening_dataset(
        receptor,
        ligands,
        dataset,
        "target",
        (1.0, 2.0, 3.0),
    )

    assert ligand_ids == ["ligand_one", "ligand_one_1"]
    assert (dataset / "proteins" / "target.pdb").is_file()
    assert (dataset / "ligands" / "target" / "ligand_one.sdf").is_file()
    assert json.loads((dataset / "pocket_centers.json").read_text()) == {
        "target": [1.0, 2.0, 3.0]
    }


def test_write_matcha_outputs_prefers_openmm_and_converts_affinity(tmp_path):
    """The highest-ranked final-stage pose should match Hedgehog's score contract."""
    predictions = tmp_path / "predictions.sdf.gz"
    with gzip.open(predictions, "wt") as output:
        writer = Chem.SDWriter(output)
        for stage, final_score, affinity in (
            ("Model 2", 9.0, 4.0),
            ("OpenMM", 1.0, 6.5),
            ("OpenMM", 2.0, 7.25),
        ):
            mol = _molecule("source")
            mol.SetProp("ligand_name", "ligand_1")
            mol.SetProp("stage", stage)
            mol.SetProp("final_score", str(final_score))
            mol.SetProp("cnn_affinity", str(affinity))
            writer.write(mol)
        writer.close()

    count = _write_matcha_outputs(predictions, tmp_path / "run")

    assert count == 1
    best_path = tmp_path / "run" / "best_poses" / "ligand_1.sdf"
    best = next(
        mol
        for mol in Chem.SDMolSupplier(str(best_path), removeHs=False, sanitize=False)
        if mol is not None
    )
    assert best.GetProp("stage") == "OpenMM"
    assert float(best.GetProp("cnn_affinity")) == 7.25
    assert float(best.GetProp("minimizedAffinity")) == -7.25


def test_write_matcha_outputs_prefers_lowest_balmus_score(tmp_path):
    """BALMUS should outrank prior stages and expose its loss-like score."""
    predictions = tmp_path / "predictions.sdf.gz"
    with gzip.open(predictions, "wt") as output:
        writer = Chem.SDWriter(output)
        for stage, final_score, balmus_score in (
            ("Model 1", 100.0, None),
            ("BALMUS", -4.0, 4.0),
            ("BALMUS", -2.5, 2.5),
        ):
            mol = _molecule("source")
            mol.SetProp("ligand_name", "ligand_1")
            mol.SetProp("stage", stage)
            mol.SetProp("final_score", str(final_score))
            if balmus_score is not None:
                mol.SetProp("balmus_score", str(balmus_score))
            writer.write(mol)
        writer.close()

    assert _write_matcha_outputs(predictions, tmp_path / "run") == 1

    best_path = tmp_path / "run" / "best_poses" / "ligand_1.sdf"
    best = next(
        mol
        for mol in Chem.SDMolSupplier(str(best_path), removeHs=False, sanitize=False)
        if mol is not None
    )
    assert best.GetProp("stage") == "BALMUS"
    assert float(best.GetProp("balmus_score")) == 2.5
    assert float(best.GetProp("minimizedAffinity")) == 2.5
    assert best.GetProp("source_score_property") == "balmus_score"
