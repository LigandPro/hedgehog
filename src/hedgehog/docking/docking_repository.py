"""Run the LigandPro/docking screening API behind Hedgehog's Matcha contract."""

import argparse
import concurrent.futures
import gzip
import json
import shutil
import subprocess
import time
from collections import defaultdict
from pathlib import Path

from rdkit import Chem, rdBase


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


def _count_sdf_records(path: Path) -> int:
    with path.open("rb") as sdf_file:
        return sum(line.strip() == b"$$$$" for line in sdf_file)


def _run_gnina_minimization(
    gnina_bin: Path,
    receptor: Path,
    input_sdf: Path,
    output_sdf: Path,
    cpu: int,
) -> None:
    command = [
        str(gnina_bin),
        "--receptor",
        str(receptor),
        "--ligand",
        str(input_sdf),
        "--cnn_scoring",
        "none",
        "--minimize",
        "--cpu",
        str(cpu),
        "--out",
        str(output_sdf),
    ]
    completed = subprocess.run(
        command,
        check=False,
        text=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
    )
    if completed.returncode != 0 or not output_sdf.is_file():
        error = completed.stderr.strip() or "GNINA did not create an output SDF"
        raise RuntimeError(f"GNINA minimization failed for {input_sdf.name}: {error}")

    expected_records = _count_sdf_records(input_sdf)
    output_records = _count_sdf_records(output_sdf)
    if output_records != expected_records:
        raise RuntimeError(
            f"GNINA minimized {output_records}/{expected_records} poses "
            f"for {input_sdf.name}"
        )


def _best_minimized_pose(path: Path) -> Chem.Mol:
    with rdBase.BlockLogs():
        poses = [
            mol
            for mol in Chem.SDMolSupplier(str(path), removeHs=False, sanitize=False)
            if mol is not None and _score_value(mol, "minimizedAffinity") is not None
        ]
    if not poses:
        raise ValueError(f"GNINA output has no minimizedAffinity values: {path}")
    return min(poses, key=lambda mol: _score_value(mol, "minimizedAffinity"))


def _write_matcha_outputs(
    predictions_path: Path,
    run_dir: Path,
    *,
    receptor: Path,
    gnina_bin: Path,
    gnina_cpu: int = 1,
    gnina_jobs: int = 1,
) -> int:
    poses_by_ligand: dict[str, list[Chem.Mol]] = defaultdict(list)
    with gzip.open(predictions_path, "rb") as sdf_file:
        supplier = Chem.ForwardSDMolSupplier(sdf_file, removeHs=False, sanitize=False)
        for mol in supplier:
            if mol is None or not mol.HasProp("ligand_name"):
                continue
            poses_by_ligand[mol.GetProp("ligand_name")].append(Chem.Mol(mol))

    generated_dir = run_dir / "generated_poses"
    all_dir = run_dir / "all_poses"
    best_dir = run_dir / "best_poses"
    failures_path = run_dir / "minimization_failures.json"
    failures_path.unlink(missing_ok=True)
    for directory in (generated_dir, all_dir, best_dir):
        if directory.exists():
            shutil.rmtree(directory)
        directory.mkdir(parents=True)

    minimization_jobs: list[tuple[str, Path, Path]] = []
    for ligand_id, poses in sorted(poses_by_ligand.items()):
        generated_path = generated_dir / f"{ligand_id}_poses.sdf"
        minimized_path = all_dir / f"{ligand_id}_poses.sdf"
        all_writer = Chem.SDWriter(str(generated_path))
        for pose in poses:
            all_writer.write(pose)
        all_writer.close()
        minimization_jobs.append((ligand_id, generated_path, minimized_path))

    def minimize(job: tuple[str, Path, Path]) -> tuple[str, Path]:
        ligand_id, generated_path, minimized_path = job
        _run_gnina_minimization(
            gnina_bin,
            receptor,
            generated_path,
            minimized_path,
            max(1, gnina_cpu),
        )
        return ligand_id, minimized_path

    minimized_outputs: list[tuple[str, Path]] = []
    failures: list[dict[str, str]] = []
    total_jobs = len(minimization_jobs)
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, gnina_jobs)) as pool:
        futures = {pool.submit(minimize, job): job[0] for job in minimization_jobs}
        for completed, future in enumerate(
            concurrent.futures.as_completed(futures),
            start=1,
        ):
            ligand_id = futures[future]
            try:
                minimized_outputs.append(future.result())
            except Exception as exc:
                failures.append(
                    {
                        "ligand_id": ligand_id,
                        "phase": "gnina_minimization",
                        "error": str(exc),
                    }
                )
                print(
                    f"WARNING: Matcha post-minimization failed for {ligand_id}: {exc}"
                )
            if completed % 25 == 0 or completed == total_jobs:
                print(f"GNINA minimized Matcha poses: {completed}/{total_jobs}")

    written = 0
    for ligand_id, minimized_path in sorted(minimized_outputs):
        try:
            best = _best_minimized_pose(minimized_path)
        except Exception as exc:
            failures.append(
                {"ligand_id": ligand_id, "phase": "best_pose", "error": str(exc)}
            )
            print(f"WARNING: Matcha best-pose selection failed for {ligand_id}: {exc}")
            continue
        best.SetProp("_Name", ligand_id)
        best.SetProp("mol_idx", ligand_id)
        best.SetProp("source_score_property", "minimizedAffinity")

        writer = Chem.SDWriter(str(best_dir / f"{ligand_id}.sdf"))
        writer.write(best)
        writer.close()
        written += 1

    if failures:
        failures_path.write_text(
            json.dumps(failures, indent=2) + "\n", encoding="utf-8"
        )
    if not written:
        raise RuntimeError("Matcha post-processing produced no usable best poses")
    return written


def _screening_command(
    args: argparse.Namespace, dataset_root: Path, results_dir: Path
) -> list[str]:
    uv_executable = shutil.which("uv")
    if uv_executable is None:
        raise FileNotFoundError("The docking backend requires uv on PATH.")
    training_config = (
        args.checkpoint_root / args.checkpoint_run / "config.yaml"
    ).resolve()
    return [
        uv_executable,
        "run",
        "--project",
        str(args.repo),
        "screening",
        f"user.results={args.checkpoint_root}",
        f"user.cache={dataset_root.parent / 'cache'}",
        f"inference.stages.docking-1.training_config={training_config}",
        f"inference.stages.docking-1.exp_name={args.checkpoint_run}",
        "inference.stages.docking-1.integrator.name=single_step_coord",
        "inference.stages.docking-1.integrator.num_steps=1",
        "inference.stages.docking-1.collect_extra_scores=false",
        "+inference.stages.docking-2.skip=true",
        "+inference.stages.openmm.skip=true",
        "+inference.stages.gnina.skip=true",
        "+inference.stages.balmus.skip=true",
        f"datasets.dekois2.root={dataset_root}",
        f"datasets.dekois2.target={args.target_name}",
        f"datasets.dekois2.pocket_centers={dataset_root / 'pocket_centers.json'}",
        f"screening.inference_results={results_dir}",
        "screening.skip_existing=false",
    ]


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--receptor", type=Path, required=True)
    parser.add_argument("--ligands", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--run-name", default="matcha_run")
    parser.add_argument("--checkpoint-root", type=Path, required=True)
    parser.add_argument("--checkpoint-run", required=True)
    parser.add_argument("--target-name", default="hedgehog_target")
    parser.add_argument("--autobox-ligand", type=Path)
    parser.add_argument("--center-x", type=float)
    parser.add_argument("--center-y", type=float)
    parser.add_argument("--center-z", type=float)
    parser.add_argument("--gnina-bin", type=Path, required=True)
    parser.add_argument("--gnina-cpu", type=int, default=1)
    parser.add_argument("--gnina-jobs", type=int, default=1)
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
    preparation_started = time.perf_counter()
    ligand_ids = _prepare_screening_dataset(
        args.receptor,
        args.ligands,
        dataset_root,
        args.target_name,
        center,
    )
    dataset_preparation_seconds = time.perf_counter() - preparation_started
    command = _screening_command(args, dataset_root, results_dir)
    screening_started = time.perf_counter()
    subprocess.run(command, cwd=args.repo, check=True)
    screening_seconds = time.perf_counter() - screening_started

    predictions_path = results_dir / "dekois2" / args.target_name / "predictions.sdf.gz"
    if not predictions_path.exists():
        raise FileNotFoundError(
            f"Docking screening output is missing: {predictions_path}"
        )
    conversion_started = time.perf_counter()
    written = _write_matcha_outputs(
        predictions_path,
        run_dir,
        receptor=args.receptor,
        gnina_bin=args.gnina_bin,
        gnina_cpu=args.gnina_cpu,
        gnina_jobs=args.gnina_jobs,
    )
    output_conversion_seconds = time.perf_counter() - conversion_started

    timing = {
        "backend": "LigandPro/docking screening",
        "input_ligands": len(ligand_ids),
        "output_ligands": written,
        "dataset_preparation_sec": dataset_preparation_seconds,
        "screening_sec": screening_seconds,
        "output_conversion_sec": output_conversion_seconds,
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
