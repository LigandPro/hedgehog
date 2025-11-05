"""End-to-end tests for HEDGE pipeline execution."""

import subprocess
from pathlib import Path

import pandas as pd
import pytest


@pytest.mark.e2e
@pytest.mark.slow
def test_hedge_run_with_descriptors_only(test_molecules_csv, temp_output_dir, fixtures_dir):
    """Test that 'hedge run' executes successfully with descriptors stage."""
    # Create a custom config for this test that points to temp output
    test_config = temp_output_dir / "config.yml"
    config_content = f"""generated_mols_path: {test_molecules_csv}
folder_to_save: {temp_output_dir}
n_jobs: 2
sample_size: 10
save_sampled_mols: True

config_descriptors: {fixtures_dir}/test_config_descriptors.yml
config_structFilters: {fixtures_dir}/test_config_structFilters.yml
config_synthesis: {fixtures_dir}/test_config_synthesis.yml
config_docking: {fixtures_dir}/test_config_docking.yml
"""
    test_config.write_text(config_content)

    # Temporarily update the main config path to point to our test config
    main_config_path = Path("src/hedge/configs/config.yml")
    backup_path = Path("src/hedge/configs/config.yml.backup")

    # Backup original config
    if main_config_path.exists():
        main_config_path.rename(backup_path)

    try:
        # Copy test config to main config location
        main_config_path.write_text(config_content)

        # Run hedge with test molecules
        result = subprocess.run(
            ["uv", "run", "hedge", "run", "--mols", str(test_molecules_csv)],
            capture_output=True,
            text=True,
            timeout=60,
        )

        # Check that command executed successfully
        assert result.returncode == 0, f"Command failed with stderr: {result.stderr}"

        # Check that there are no ERROR messages in output
        assert "ERROR" not in result.stdout.upper(), f"Found ERROR in output: {result.stdout}"

        # Check that pipeline completed message appears
        assert "Pipeline completed" in result.stdout or "Ligand Pro thanks you" in result.stdout

        # Check that output files were created
        assert (temp_output_dir / "sampledMols.csv").exists(), "sampledMols.csv not created"

        # Check descriptors output (since it's enabled in test config)
        descriptors_dir = temp_output_dir / "Descriptors"
        assert descriptors_dir.exists(), "Descriptors directory not created"

        # Verify molecules were processed
        sampled_df = pd.read_csv(temp_output_dir / "sampledMols.csv")
        assert len(sampled_df) > 0, "No molecules in sampledMols.csv"
        assert "smiles" in sampled_df.columns
        assert "model_name" in sampled_df.columns

    finally:
        # Restore original config
        if backup_path.exists():
            if main_config_path.exists():
                main_config_path.unlink()
            backup_path.rename(main_config_path)


@pytest.mark.e2e
def test_hedge_info_command():
    """Test that 'hedge info' command works."""
    result = subprocess.run(
        ["uv", "run", "hedge", "info"],
        capture_output=True,
        text=True,
        timeout=10,
    )

    assert result.returncode == 0
    assert "descriptors" in result.stdout.lower()
    assert "pipeline" in result.stdout.lower() or "stage" in result.stdout.lower()


@pytest.mark.e2e
def test_hedge_version_command():
    """Test that 'hedge version' command works."""
    result = subprocess.run(
        ["uv", "run", "hedge", "version"],
        capture_output=True,
        text=True,
        timeout=10,
    )

    assert result.returncode == 0
    assert "1.0.0" in result.stdout or "HEDGE" in result.stdout


@pytest.mark.e2e
def test_hedge_help_command():
    """Test that 'hedge --help' command works."""
    result = subprocess.run(
        ["uv", "run", "hedge", "--help"],
        capture_output=True,
        text=True,
        timeout=10,
    )

    assert result.returncode == 0
    assert "hedge" in result.stdout.lower() or "run" in result.stdout.lower()
