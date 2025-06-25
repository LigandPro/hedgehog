#!/usr/bin/env python
from pathlib import Path
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
from rich.console import Console
from rich.logging import RichHandler
from rich.progress import (
    Progress,
    SpinnerColumn,
    TextColumn,
    BarColumn,
    TaskProgressColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)
import typer

console = Console()
logging.basicConfig(
    level="INFO",
    format="%(message)s",
    handlers=[RichHandler(console=console, show_time=False, show_level=False)],
)
logger = logging.getLogger(__name__)
app = typer.Typer()

def clean_dir_name(name):
    for suffix in [
        "_10000_mols_smiles_left_mols",
        "_smiles_left_mols",
        "_left_mols",
    ]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name

def run_smina_task(
    protein_path: Path,
    ligand_path: Path,
    pocket_file: Path,
    result_subdir: Path,
    smina_path: Path,
    exhaustiveness: int,
    num_modes: int,
    cpu: int,
    seed: int,
    autobox_add: int,
):
    ligand_name = ligand_path.stem
    result_file = result_subdir / f"{protein_path.stem}__{ligand_name}.mol2"
    log_file = result_subdir / f"{protein_path.stem}__{ligand_name}.log"
    if result_file.exists():
        return None, result_subdir.name
    command = [
        str(smina_path),
        "--receptor", str(protein_path),
        "--ligand", str(ligand_path),
        "--autobox_ligand", str(pocket_file),
        "--autobox_add", str(autobox_add),
        "-o", str(result_file),
        "--log", str(log_file),
        "--exhaustiveness", str(exhaustiveness),
        "--cpu", str(cpu),
        "--seed", str(seed),
    ]
    try:
        with open(log_file, "w") as f:
            subprocess.run(command, check=True, stdout=f, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        return f"Failed: {protein_path.stem} × {ligand_name}", result_subdir.name
    return None, result_subdir.name

@app.command()
def main(
    protein_files: list[Path] = typer.Option(["PROTEIN1.pdbqt"]),
    ligand_dir: Path = typer.Option(Path("pdbqt_ready")),
    pocket_file: Path = typer.Option(Path("YOUR_LIGAND.sdf")),
    result_dir: Path = typer.Option(Path("docking_results")),
    smina_path: Path = typer.Option(Path("smina/smina.static")),
    exhaustiveness: int = typer.Option(12),
    num_modes: int = typer.Option(20),
    cpu: int = typer.Option(24),
    seed: int = typer.Option(42),
    autobox_add: int = typer.Option(6),
    n_jobs: int = typer.Option(16),
):
    result_dir.mkdir(parents=True, exist_ok=True)
    ligand_subdirs = [p for p in ligand_dir.iterdir() if p.is_dir()] or [ligand_dir]
    all_tasks = []
    cleaned_stats = []

    for protein_path in protein_files:
        protein_name = protein_path.stem
        protein_result_dir = result_dir / protein_name
        protein_result_dir.mkdir(parents=True, exist_ok=True)
        for subdir in ligand_subdirs:
            subdir_name_clean = clean_dir_name(subdir.name)
            result_subdir = protein_result_dir / subdir_name_clean
            result_subdir.mkdir(parents=True, exist_ok=True)
            ligands = sorted(subdir.rglob("*.pdbqt"))
            subdir_tasks = [
                (
                    protein_path,
                    ligand_path,
                    pocket_file,
                    result_subdir,
                    smina_path,
                    exhaustiveness,
                    num_modes,
                    cpu,
                    seed,
                    autobox_add,
                )
                for ligand_path in ligands
                if not (result_subdir / f"{protein_name}__{ligand_path.stem}.mol2").exists()
            ]
            all_tasks.extend(subdir_tasks)
            total = len(ligands)
            done = total - len(subdir_tasks)
            cleaned_stats.append((protein_name, subdir_name_clean, done, total))


    max_prot_len = max((len(x[0]) for x in cleaned_stats), default=0)
    max_dir_len = max((len(x[1]) for x in cleaned_stats), default=0)
    for prot, ligdir, done, total in cleaned_stats:
        console.print(f"{prot:<{max_prot_len}}  {ligdir:<{max_dir_len}}:  {done:>5}/{total:<5} completed")

    global_total = len(all_tasks)
    if global_total == 0:
        console.print("Nothing to run. All results exist.")
        return

    with Progress(
        SpinnerColumn(),
        TextColumn("Docking progress: {task.fields[dir]}"),
        BarColumn(bar_width=None),
        TaskProgressColumn("{task.completed}/{task.total}"),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=console,
    ) as progress:
        task_id = progress.add_task("Total", total=global_total, dir="")
        with ProcessPoolExecutor(max_workers=n_jobs) as executor:
            futures = {executor.submit(run_smina_task, *t): t[3].name for t in all_tasks}
            for future in as_completed(futures):
                result, dirname = future.result()
                progress.update(task_id, dir=dirname)
                if result:
                    logger.info(result)
                progress.advance(task_id)
    console.print("All docking complete.")

if __name__ == "__main__":
    app()
