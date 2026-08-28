"""Run the LigandPro/docking screening API behind Hedgehog's Matcha contract."""

import argparse
import gzip
import json
import os
import shutil
import subprocess
import time
from collections import defaultdict
from pathlib import Path

from rdkit import Chem


def _safe_ligand_id(name: str, index: int, seen: dict[str, int]) -> str:
    ligand_id = name.strip().replace(" ", "_").replace("/", "_").replace("\\", "_")
    ligand_id = ligand_id or f"mol_{index}"
    duplicate_index = seen.get(ligand_id, 0)
    seen[ligand_id] = duplicate_index + 1
    return ligand_id if duplicate_index == 0 else f"{ligand_id}_{duplicate_index}"


def _read_ligands(path: Path) -> list[tuple[str, Chem.Mol]]:
    sources = sorted(path.glob("*.sdf")) if path.is_dir() else [path]
    ligands: list[tuple[str, Chem.Mol]] = []
    seen: dict[str, int] = {}
    record_index = 0

    for source in sources:
        supplier = Chem.SDMolSupplier(str(source), removeHs=False, sanitize=False)
        for mol in supplier:
            if mol is None:
                record_index += 1
                continue
            raw_name = mol.GetProp("_Name") if mol.HasProp("_Name") else source.stem
            ligand_id = _safe_ligand_id(raw_name, record_index, seen)
            ligands.append((ligand_id, mol))
            record_index += 1
    return ligands


def _autobox_center(path: Path) -> tuple[float, float, float]:
    supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=False)
    mol = next((candidate for candidate in supplier if candidate is not None), None)
    if mol is None or mol.GetNumConformers() == 0:
        raise ValueError(f"Cannot determine a docking center from {path}")
    coordinates = mol.GetConformer().GetPositions()
    center = coordinates.mean(axis=0)
    return float(center[0]), float(center[1]), float(center[2])


def _prepare_screening_dataset(
    receptor: Path,
    ligands_path: Path,
    dataset_root: Path,
    target_name: str,
    center: tuple[float, float, float],
) -> list[str]:
    if dataset_root.exists():
        shutil.rmtree(dataset_root)

    proteins_dir = dataset_root / "proteins"
    ligands_dir = dataset_root / "ligands" / target_name
    proteins_dir.mkdir(parents=True)
    ligands_dir.mkdir(parents=True)
    shutil.copyfile(receptor, proteins_dir / f"{target_name}.pdb")

    ligand_ids = []
    for ligand_id, mol in _read_ligands(ligands_path):
        destination = ligands_dir / f"{ligand_id}.sdf"
        writer = Chem.SDWriter(str(destination))
        writer.SetKekulize(False)
        writer.write(mol)
        writer.close()
        ligand_ids.append(ligand_id)

    if not ligand_ids:
        raise ValueError(f"No readable ligands found in {ligands_path}")

    centers_path = dataset_root / "pocket_centers.json"
    centers_path.write_text(
        json.dumps({target_name: list(center)}, indent=2),
        encoding="utf-8",
    )
    (dataset_root / "ligand_ids.json").write_text(
        json.dumps(ligand_ids, indent=2),
        encoding="utf-8",
    )
    return ligand_ids


def _score_value(mol: Chem.Mol, property_name: str) -> float | None:
    if not mol.HasProp(property_name):
        return None
    try:
        return float(mol.GetProp(property_name))
    except ValueError:
        return None


def _select_best_pose(poses: list[Chem.Mol]) -> Chem.Mol:
    stage_priority = {"Model 1": 1, "Model 2": 2, "OpenMM": 3, "BALMUS": 4}
    best_stage = max(
        (stage_priority.get(mol.GetProp("stage"), 0) for mol in poses),
        default=0,
    )
    candidates = [
        mol
        for mol in poses
        if stage_priority.get(mol.GetProp("stage"), 0) == best_stage
    ]

    def ranking_key(mol: Chem.Mol) -> tuple[float, float]:
        final_score = _score_value(mol, "final_score")
        cnn_affinity = _score_value(mol, "cnn_affinity")
        return (
            float("-inf") if final_score is None else final_score,
            float("-inf") if cnn_affinity is None else cnn_affinity,
        )

    return max(candidates, key=ranking_key)


def _write_matcha_outputs(predictions_path: Path, run_dir: Path) -> int:
    poses_by_ligand: dict[str, list[Chem.Mol]] = defaultdict(list)
    with gzip.open(predictions_path, "rb") as sdf_file:
        supplier = Chem.ForwardSDMolSupplier(sdf_file, removeHs=False, sanitize=False)
        for mol in supplier:
            if mol is None or not mol.HasProp("ligand_name"):
                continue
            poses_by_ligand[mol.GetProp("ligand_name")].append(Chem.Mol(mol))

    best_dir = run_dir / "best_poses"
    all_dir = run_dir / "all_poses"
    best_dir.mkdir(parents=True, exist_ok=True)
    all_dir.mkdir(parents=True, exist_ok=True)

    written = 0
    for ligand_id, poses in sorted(poses_by_ligand.items()):
        all_writer = Chem.SDWriter(str(all_dir / f"{ligand_id}_poses.sdf"))
        for pose in poses:
            all_writer.write(pose)
        all_writer.close()

        best = _select_best_pose(poses)
        cnn_affinity = _score_value(best, "cnn_affinity")
        balmus_score = _score_value(best, "balmus_score")
        if balmus_score is not None:
            minimized_affinity = balmus_score
            source_score_property = "balmus_score"
        elif cnn_affinity is not None:
            minimized_affinity = -cnn_affinity
            source_score_property = "cnn_affinity"
        else:
            raise ValueError(
                f"Best pose for {ligand_id} has neither balmus_score nor cnn_affinity"
            )
        best.SetProp("_Name", ligand_id)
        best.SetProp("mol_idx", ligand_id)
        best.SetProp("affinity", str(minimized_affinity))
        best.SetProp("minimizedAffinity", str(minimized_affinity))
        best.SetProp("source_score_property", source_score_property)

        writer = Chem.SDWriter(str(best_dir / f"{ligand_id}.sdf"))
        writer.write(best)
        writer.close()
        written += 1
    return written


def _screening_command(
    args: argparse.Namespace, dataset_root: Path, results_dir: Path
) -> list[str]:
    return [
        args.uv_bin,
        "run",
        "--project",
        str(args.repo),
        "screening",
        f"user.results={args.checkpoint_root}",
        f"user.cache={dataset_root.parent / 'cache'}",
        f"inference.stages.docking-1.training_config={args.training_config}",
        f"inference.stages.docking-1.exp_name={args.checkpoint_run}",
        f"inference.stages.docking-1.checkpoint={args.checkpoint_name}",
        "inference.stages.docking-1.integrator.name=single_step_coord",
        "inference.stages.docking-1.integrator.num_steps=1",
        "inference.stages.docking-1.collect_extra_scores=false",
        "+inference.stages.docking-2.skip=true",
        "+inference.stages.openmm.skip=true",
        "+inference.stages.gnina.skip=true",
        "+inference.stages.balmus.skip=false",
        f"datasets.dekois2.root={dataset_root}",
        f"datasets.dekois2.target={args.target_name}",
        f"datasets.dekois2.pocket_centers={dataset_root / 'pocket_centers.json'}",
        f"screening.inference_results={results_dir}",
        "screening.skip_existing=false",
        f"screening.data_workers={args.data_workers}",
        f"pipe.samples_per_complex={args.n_samples}",
        f"pipe.sample_timeout_seconds={args.sample_timeout_seconds}",
        f"inference.concurrency={args.concurrency}",
        f"inference.batching.batch_size={args.batch_size}",
    ]


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--uv-bin", required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--receptor", type=Path, required=True)
    parser.add_argument("--ligands", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--run-name", default="matcha_run")
    parser.add_argument("--training-config", type=Path, required=True)
    parser.add_argument("--checkpoint-root", type=Path, required=True)
    parser.add_argument("--checkpoint-run", required=True)
    parser.add_argument("--checkpoint-name", default="checkpoint-latest")
    parser.add_argument("--target-name", default="hedgehog_target")
    parser.add_argument("--autobox-ligand", type=Path)
    parser.add_argument("--center-x", type=float)
    parser.add_argument("--center-y", type=float)
    parser.add_argument("--center-z", type=float)
    parser.add_argument("--n-samples", type=int, default=20)
    parser.add_argument("--sample-timeout-seconds", type=float, default=180.0)
    parser.add_argument("--data-workers", type=int, default=16)
    parser.add_argument("--concurrency", type=int, default=64)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--gpus")
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    run_dir = args.out / args.run_name
    work_dir = run_dir / "work"
    dataset_root = work_dir / "screening_dataset"
    results_dir = work_dir / "screening_results"

    if args.autobox_ligand:
        center = _autobox_center(args.autobox_ligand)
    elif None not in (args.center_x, args.center_y, args.center_z):
        center = float(args.center_x), float(args.center_y), float(args.center_z)
    else:
        raise ValueError(
            "The docking repository backend requires autobox_ligand or a center"
        )

    started = time.perf_counter()
    ligand_ids = _prepare_screening_dataset(
        args.receptor,
        args.ligands,
        dataset_root,
        args.target_name,
        center,
    )
    command = _screening_command(args, dataset_root, results_dir)
    environment = os.environ.copy()
    if args.gpus:
        environment["CUDA_VISIBLE_DEVICES"] = args.gpus
    subprocess.run(command, cwd=args.repo, env=environment, check=True)

    predictions_path = results_dir / "dekois2" / args.target_name / "predictions.sdf.gz"
    if not predictions_path.exists():
        raise FileNotFoundError(
            f"Docking screening output is missing: {predictions_path}"
        )
    written = _write_matcha_outputs(predictions_path, run_dir)
    if written == 0:
        raise RuntimeError("Docking screening produced no usable poses")

    timing = {
        "backend": "LigandPro/docking screening",
        "input_ligands": len(ligand_ids),
        "output_ligands": written,
        "total_sec": time.perf_counter() - started,
        "command": command,
    }
    (run_dir / "run_timing.json").write_text(
        json.dumps(timing, indent=2),
        encoding="utf-8",
    )
    print(f"Docking screening complete: {written}/{len(ligand_ids)} ligands")


if __name__ == "__main__":
    main()
