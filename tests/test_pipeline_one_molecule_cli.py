"""End-to-end CLI smoke for one molecule across all pipeline stages."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest
import yaml

PASSING_SMILES = "O=C(c1ccc(F)cc1)N1CCC(N2C(=O)CCc3ccccc32)CC1"
SUCCESS_MARKER = "Pipeline completed: 7/7 stages successful"


def _load_yaml(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def _dump_yaml(path: Path, payload: dict) -> None:
    path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")


def _build_single_core_config(
    repo_root: Path, tmp_path: Path, input_csv: Path, output_dir: Path
) -> Path:
    default_cfg_dir = repo_root / "src" / "hedgehog" / "configs"
    cfg_dir = tmp_path / "configs"
    cfg_dir.mkdir(parents=True, exist_ok=True)

    mol_prep_cfg = _load_yaml(default_cfg_dir / "config_mol_prep.yml")
    mol_prep_cfg["n_jobs"] = 1
    mol_prep_cfg_path = cfg_dir / "config_mol_prep.yml"
    _dump_yaml(mol_prep_cfg_path, mol_prep_cfg)

    descriptors_cfg = _load_yaml(default_cfg_dir / "config_descriptors.yml")
    descriptors_cfg["n_jobs"] = 1
    descriptors_cfg["batch_size"] = 1
    descriptors_cfg_path = cfg_dir / "config_descriptors.yml"
    _dump_yaml(descriptors_cfg_path, descriptors_cfg)

    struct_filters_cfg = _load_yaml(default_cfg_dir / "config_structFilters.yml")
    struct_filters_cfg["parse_input_n_jobs"] = 1
    struct_filters_cfg["common_alerts_auto_n_jobs"] = False
    struct_filters_cfg["common_alerts_small_input_n_jobs"] = 1
    struct_filters_cfg["common_alerts_large_input_n_jobs"] = 1
    struct_filters_cfg_path = cfg_dir / "config_structFilters.yml"
    _dump_yaml(struct_filters_cfg_path, struct_filters_cfg)

    synthesis_cfg = _load_yaml(default_cfg_dir / "config_synthesis.yml")
    synthesis_cfg["n_jobs"] = 1
    synthesis_cfg_path = cfg_dir / "config_synthesis.yml"
    _dump_yaml(synthesis_cfg_path, synthesis_cfg)

    docking_cfg = _load_yaml(default_cfg_dir / "config_docking.yml")
    docking_cfg["tools"] = "gnina"
    docking_cfg["run_in_background"] = False
    docking_cfg["gnina_parallel_jobs"] = 1
    docking_cfg["receptor_pdb"] = str(
        repo_root / "src" / "hedgehog" / "configs" / "examples" / "7EW9_apo.pdb"
    )
    docking_cfg["gnina_config"]["autobox_ligand"] = str(
        repo_root / "src" / "hedgehog" / "configs" / "examples" / "05C_from_7EW9.sdf"
    )
    docking_cfg["gnina_config"]["cpu"] = 1
    docking_cfg["gnina_config"]["no_gpu"] = True
    docking_cfg_path = cfg_dir / "config_docking.yml"
    _dump_yaml(docking_cfg_path, docking_cfg)

    docking_filters_cfg = _load_yaml(default_cfg_dir / "config_docking_filters.yml")
    docking_filters_cfg["n_jobs"] = 1
    docking_filters_cfg_path = cfg_dir / "config_docking_filters.yml"
    _dump_yaml(docking_filters_cfg_path, docking_filters_cfg)

    moleval_cfg = _load_yaml(default_cfg_dir / "config_moleval.yml")
    moleval_cfg["run"] = True
    moleval_cfg["n_jobs"] = 1
    moleval_cfg_path = cfg_dir / "config_moleval.yml"
    _dump_yaml(moleval_cfg_path, moleval_cfg)

    master_cfg = _load_yaml(default_cfg_dir / "config.yml")
    master_cfg["generated_mols_path"] = str(input_csv)
    master_cfg["folder_to_save"] = str(output_dir)
    master_cfg["n_jobs"] = 1
    master_cfg["sample_size"] = 1
    master_cfg["save_sampled_mols"] = True
    master_cfg["ligand_preparation_tool"] = None
    master_cfg["protein_preparation_tool"] = None
    master_cfg["config_mol_prep"] = str(mol_prep_cfg_path)
    master_cfg["config_descriptors"] = str(descriptors_cfg_path)
    master_cfg["config_structFilters"] = str(struct_filters_cfg_path)
    master_cfg["config_synthesis"] = str(synthesis_cfg_path)
    master_cfg["config_docking"] = str(docking_cfg_path)
    master_cfg["config_docking_filters"] = str(docking_filters_cfg_path)
    master_cfg["config_moleval"] = str(moleval_cfg_path)

    master_cfg_path = cfg_dir / "config.yml"
    _dump_yaml(master_cfg_path, master_cfg)
    return master_cfg_path


@pytest.mark.integration
@pytest.mark.slow
def test_cli_full_pipeline_for_single_molecule_single_core(tmp_path: Path):
    """Run `uv run hedgehog` end-to-end for one known passing molecule."""
    repo_root = Path(__file__).resolve().parents[1]

    if shutil.which("uv") is None:
        pytest.skip("uv is not available in PATH")
    if shutil.which("gnina") is None:
        pytest.skip("gnina is not available in PATH")

    aizynth_cfg = (
        repo_root
        / "modules"
        / "retrosynthesis"
        / "aizynthfinder"
        / "public"
        / "config.yml"
    )
    if not aizynth_cfg.exists():
        pytest.skip(
            "AiZynthFinder is not configured; run `hedgehog setup aizynthfinder`"
        )

    input_csv = tmp_path / "one_molecule.csv"
    input_csv.write_text(
        f"smiles,model_name\n{PASSING_SMILES},smoke_model\n", encoding="utf-8"
    )
    output_dir = tmp_path / "run_one_molecule"
    config_path = _build_single_core_config(repo_root, tmp_path, input_csv, output_dir)

    env = os.environ.copy()
    env.update(
        {
            "MOLSCORE_NJOBS": "1",
            "SLURM_CPUS_PER_TASK": "1",
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
            "VECLIB_MAXIMUM_THREADS": "1",
            "AIZYNTH_NPROC": "1",
            "HEDGEHOG_AUTO_INSTALL": "1",
            "HEDGEHOG_PLAIN_OUTPUT": "1",
            "HEDGEHOG_STRICT_RASCORE": "1",
        }
    )

    cmd = [
        "uv",
        "run",
        "hedgehog",
        "--config",
        str(config_path),
        "--mols",
        str(input_csv),
        "--out",
        str(output_dir),
        "--auto-install",
    ]
    proc = subprocess.run(
        cmd,
        cwd=repo_root,
        env=env,
        capture_output=True,
        text=True,
        timeout=900,
        check=False,
    )
    assert proc.returncode == 0, (
        "CLI pipeline run failed.\n"
        f"stdout:\n{proc.stdout[-4000:]}\n"
        f"stderr:\n{proc.stderr[-4000:]}"
    )

    run_logs = sorted(output_dir.glob("run_*.log"))
    assert run_logs, "Expected run log to be created"

    log_text = run_logs[-1].read_text(encoding="utf-8", errors="ignore")
    assert SUCCESS_MARKER in log_text
    expected_stage_lines = [
        "mol_prep: \u2713 COMPLETED",
        "descriptors: \u2713 COMPLETED",
        "struct_filters: \u2713 COMPLETED",
        "synthesis: \u2713 COMPLETED",
        "docking: \u2713 COMPLETED",
        "docking_filters: \u2713 COMPLETED",
        "final_descriptors: \u2713 COMPLETED",
    ]
    for expected_line in expected_stage_lines:
        assert expected_line in log_text

    assert "Failed to load RAScore model" not in log_text
    assert "falling back to legacy worker backend" not in log_text

    assert "MolPrep workers: 1" in log_text
    assert "Descriptors workers: 1" in log_text
    assert "cpu_per_process=1, parallel_jobs=1" in log_text

    final_csv = output_dir / "output" / "final_molecules.csv"
    assert final_csv.exists(), "Expected final_molecules.csv to exist"
    final_df = pd.read_csv(final_csv)
    assert len(final_df) == 1

    assert not (output_dir / ".RUN_INCOMPLETE").exists()
