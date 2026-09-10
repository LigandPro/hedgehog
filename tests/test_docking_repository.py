"""Tests for the LigandPro/docking screening adapter."""

import gzip
import json
from argparse import Namespace

from rdkit import Chem
from rdkit.Chem import AllChem

from hedgehog.docking.docking_repository import (
    _prepare_screening_dataset,
    _screening_command,
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


def test_screening_skips_internal_rankers_before_empirical_rescoring(
    tmp_path, monkeypatch
):
    """The repository should generate poses without BALMUS or CNN score ranking."""
    checkpoint_root = tmp_path / "checkpoints"
    training_config = checkpoint_root / "run" / "config.yaml"
    training_config.parent.mkdir(parents=True)
    training_config.write_text("seed: 42\n", encoding="utf-8")
    monkeypatch.setattr(
        "hedgehog.docking.docking_repository.shutil.which",
        lambda executable: f"/usr/bin/{executable}",
    )
    args = Namespace(
        repo=tmp_path / "docking",
        checkpoint_root=checkpoint_root,
        checkpoint_run="run",
        target_name="target",
    )

    command = _screening_command(args, tmp_path / "dataset", tmp_path / "results")

    assert "+inference.stages.gnina.skip=true" in command
    assert "+inference.stages.balmus.skip=true" in command
    assert (
        f"inference.stages.docking-1.training_config={training_config.resolve()}"
        in command
    )
    assert not any(
        option.startswith(
            (
                "screening.data_workers=",
                "pipe.samples_per_complex=",
                "pipe.sample_timeout_seconds=",
                "inference.concurrency=",
                "inference.batching.batch_size=",
                "inference.stages.docking-1.checkpoint=",
            )
        )
        for option in command
    )


def _fake_gnina_minimization(_gnina_bin, _receptor, input_sdf, output_sdf, _cpu):
    writer = Chem.SDWriter(str(output_sdf))
    for mol in Chem.SDMolSupplier(str(input_sdf), removeHs=False, sanitize=False):
        if mol is None:
            continue
        mol.SetProp("minimizedAffinity", mol.GetProp("expected_minimized_affinity"))
        writer.write(mol)
    writer.close()


def test_write_matcha_outputs_selects_lowest_gnina_affinity(tmp_path, monkeypatch):
    """Every generated pose should be minimized and ranked by GNINA affinity."""
    monkeypatch.setattr(
        "hedgehog.docking.docking_repository._run_gnina_minimization",
        _fake_gnina_minimization,
    )
    predictions = tmp_path / "predictions.sdf.gz"
    with gzip.open(predictions, "wt") as output:
        writer = Chem.SDWriter(output)
        for stage, final_score, affinity in (
            ("Model 2", 9.0, -6.0),
            ("OpenMM", 1.0, -8.0),
            ("OpenMM", 2.0, -7.0),
        ):
            mol = _molecule("source")
            mol.SetProp("ligand_name", "ligand_1")
            mol.SetProp("stage", stage)
            mol.SetProp("final_score", str(final_score))
            mol.SetProp("expected_minimized_affinity", str(affinity))
            writer.write(mol)
        writer.close()

    count = _write_matcha_outputs(
        predictions,
        tmp_path / "run",
        receptor=tmp_path / "receptor.pdb",
        gnina_bin=tmp_path / "gnina",
    )

    assert count == 1
    best_path = tmp_path / "run" / "best_poses" / "ligand_1.sdf"
    best = next(
        mol
        for mol in Chem.SDMolSupplier(str(best_path), removeHs=False, sanitize=False)
        if mol is not None
    )
    assert best.GetProp("stage") == "OpenMM"
    assert float(best.GetProp("minimizedAffinity")) == -8.0
    assert best.GetProp("source_score_property") == "minimizedAffinity"
    assert (tmp_path / "run" / "generated_poses" / "ligand_1_poses.sdf").is_file()
    assert (tmp_path / "run" / "all_poses" / "ligand_1_poses.sdf").is_file()


def test_write_matcha_outputs_keeps_balmus_only_as_metadata(tmp_path, monkeypatch):
    """BALMUS may remain metadata, but it must not select or filter the pose."""
    monkeypatch.setattr(
        "hedgehog.docking.docking_repository._run_gnina_minimization",
        _fake_gnina_minimization,
    )
    predictions = tmp_path / "predictions.sdf.gz"
    with gzip.open(predictions, "wt") as output:
        writer = Chem.SDWriter(output)
        for stage, final_score, balmus_score, affinity in (
            ("Model 1", 100.0, None, -7.0),
            ("BALMUS", -4.0, 4.0, -12.0),
            ("BALMUS", -2.5, 2.5, -8.0),
        ):
            mol = _molecule("source")
            mol.SetProp("ligand_name", "ligand_1")
            mol.SetProp("stage", stage)
            mol.SetProp("final_score", str(final_score))
            mol.SetProp("expected_minimized_affinity", str(affinity))
            if balmus_score is not None:
                mol.SetProp("balmus_score", str(balmus_score))
            writer.write(mol)
        writer.close()

    assert (
        _write_matcha_outputs(
            predictions,
            tmp_path / "run",
            receptor=tmp_path / "receptor.pdb",
            gnina_bin=tmp_path / "gnina",
        )
        == 1
    )

    best_path = tmp_path / "run" / "best_poses" / "ligand_1.sdf"
    best = next(
        mol
        for mol in Chem.SDMolSupplier(str(best_path), removeHs=False, sanitize=False)
        if mol is not None
    )
    assert best.GetProp("stage") == "BALMUS"
    assert float(best.GetProp("balmus_score")) == 4.0
    assert float(best.GetProp("minimizedAffinity")) == -12.0
    assert best.GetProp("source_score_property") == "minimizedAffinity"


def test_write_matcha_outputs_records_one_failed_ligand(tmp_path, monkeypatch):
    """A single GNINA failure must not discard other Matcha results."""

    def minimize(_gnina_bin, _receptor, input_sdf, output_sdf, _cpu):
        if input_sdf.name.startswith("bad_"):
            raise RuntimeError("invalid aromatic ligand")
        _fake_gnina_minimization(_gnina_bin, _receptor, input_sdf, output_sdf, _cpu)

    monkeypatch.setattr(
        "hedgehog.docking.docking_repository._run_gnina_minimization", minimize
    )
    predictions = tmp_path / "predictions.sdf.gz"
    with gzip.open(predictions, "wt") as output:
        writer = Chem.SDWriter(output)
        for ligand_id in ("good_ligand", "bad_ligand"):
            mol = _molecule("source")
            mol.SetProp("ligand_name", ligand_id)
            mol.SetProp("expected_minimized_affinity", "-7.0")
            writer.write(mol)
        writer.close()

    run_dir = tmp_path / "run"
    assert (
        _write_matcha_outputs(
            predictions,
            run_dir,
            receptor=tmp_path / "receptor.pdb",
            gnina_bin=tmp_path / "gnina",
        )
        == 1
    )
    assert (run_dir / "best_poses/good_ligand.sdf").is_file()
    failures = json.loads((run_dir / "minimization_failures.json").read_text())
    assert failures[0]["ligand_id"] == "bad_ligand"
    assert failures[0]["phase"] == "gnina_minimization"
