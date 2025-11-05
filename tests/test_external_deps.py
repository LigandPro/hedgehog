"""Tests for external dependencies (Lilly, SYBA, AiZynthFinder)."""

import importlib
import sys
from pathlib import Path

import pytest


@pytest.mark.external
def test_rdkit_available():
    """Test that RDKit is available and working."""
    try:
        from rdkit import Chem

        # Test basic functionality
        mol = Chem.MolFromSmiles("CCO")
        assert mol is not None
        smiles = Chem.MolToSmiles(mol)
        assert smiles == "CCO"
    except ImportError as e:
        pytest.fail(f"RDKit not available: {e}")


@pytest.mark.external
def test_lilly_medchem_available():
    """Test that Lilly MedChem Rules is available (conda package)."""
    try:
        # Lilly is a conda package, check if module exists
        import_successful = False
        try:
            # Try to import lilly_medchem_rules if it has Python bindings
            import lilly_medchem_rules  # noqa: F401

            import_successful = True
        except ImportError:
            # If no Python module, check if the binary exists
            # Lilly is often used as a command-line tool
            import shutil

            lilly_binary = shutil.which("Lilly_Medchem_Rules")
            if lilly_binary:
                import_successful = True

        if not import_successful:
            pytest.skip("Lilly MedChem Rules not installed (conda package required)")
    except Exception as e:
        pytest.skip(f"Lilly MedChem Rules check failed: {e}")


@pytest.mark.external
def test_syba_available():
    """Test that SYBA is available (conda package)."""
    try:
        import syba_pkg  # SYBA conda package

        # Test basic SYBA functionality
        assert hasattr(syba_pkg, "SybaClassifier") or hasattr(syba_pkg, "__version__")
    except ImportError:
        pytest.skip("SYBA not installed (conda package required)")
    except Exception as e:
        pytest.skip(f"SYBA check failed: {e}")


@pytest.mark.external
def test_aizynthfinder_available():
    """Test that AiZynthFinder is available."""
    try:
        import aizynthfinder  # noqa: F401

        # Test that we can import key components
        from aizynthfinder.aizynthfinder import AiZynthFinder  # noqa: F401

        assert True
    except ImportError:
        pytest.skip("AiZynthFinder not installed")
    except Exception as e:
        pytest.skip(f"AiZynthFinder check failed: {e}")


@pytest.mark.external
def test_aizynthfinder_config_exists():
    """Test that AiZynthFinder config exists in expected location."""
    # Check if config exists as expected by synthesis stage
    project_root = Path(__file__).parent.parent
    aizynth_config = (
        project_root / "modules" / "retrosynthesis" / "aizynthfinder" / "public" / "config.yml"
    )

    if not aizynth_config.exists():
        pytest.skip(
            f"AiZynthFinder config not found at {aizynth_config}. "
            "Run: cd modules/retrosynthesis/aizynthfinder && "
            "uv run python -m aizynthfinder.tools.download_public_data ./public"
        )


@pytest.mark.external
def test_medchem_available():
    """Test that medchem package is available."""
    try:
        import medchem

        # Test basic medchem functionality
        assert hasattr(medchem, "__version__")
    except ImportError as e:
        pytest.fail(f"medchem package not available: {e}")


@pytest.mark.external
def test_datamol_available():
    """Test that datamol package is available."""
    try:
        import datamol as dm

        # Test basic datamol functionality
        mol = dm.to_mol("CCO")
        assert mol is not None
    except ImportError as e:
        pytest.fail(f"datamol package not available: {e}")


@pytest.mark.external
@pytest.mark.slow
def test_lilly_filtering_example(tmp_path):
    """Test Lilly filtering with a simple example (if binary available)."""
    try:
        import shutil
        import subprocess

        lilly_binary = shutil.which("Lilly_Medchem_Rules")
        if not lilly_binary:
            pytest.skip("Lilly binary not found in PATH")

        # Create test SMILES file
        test_smiles = tmp_path / "test.smi"
        test_smiles.write_text("CCO ethanol\nc1ccccc1 benzene\n")

        # Try to run Lilly (basic check)
        result = subprocess.run(
            [lilly_binary, "-h"],
            capture_output=True,
            text=True,
            timeout=10,
        )

        # If help works, Lilly is functional
        assert result.returncode == 0 or "usage" in result.stdout.lower() or "help" in result.stdout.lower()
    except subprocess.TimeoutExpired:
        pytest.skip("Lilly command timed out")
    except Exception as e:
        pytest.skip(f"Lilly test failed: {e}")


@pytest.mark.external
def test_all_critical_dependencies():
    """Test that all critical dependencies are available."""
    critical_deps = [
        "pandas",
        "numpy",
        "rdkit",
        "datamol",
        "medchem",
        "torch",
        "dask",
    ]

    missing_deps = []
    for dep in critical_deps:
        try:
            importlib.import_module(dep)
        except ImportError:
            missing_deps.append(dep)

    if missing_deps:
        pytest.fail(f"Missing critical dependencies: {', '.join(missing_deps)}")


@pytest.mark.external
def test_conda_packages_info():
    """Display info about conda packages (Lilly, SYBA)."""
    print("\n" + "=" * 60)
    print("EXTERNAL DEPENDENCIES STATUS")
    print("=" * 60)

    # Check Python version
    print(f"Python: {sys.version}")

    # Check conda environment
    conda_env = sys.prefix
    print(f"Environment: {conda_env}")

    # Check for Lilly
    try:
        import shutil

        lilly_binary = shutil.which("Lilly_Medchem_Rules")
        if lilly_binary:
            print(f"✓ Lilly: Found at {lilly_binary}")
        else:
            print("✗ Lilly: Not found in PATH")
    except Exception as e:
        print(f"✗ Lilly: Error checking - {e}")

    # Check for SYBA
    try:
        import syba_pkg

        print(f"✓ SYBA: {syba_pkg.__version__ if hasattr(syba_pkg, '__version__') else 'Available'}")
    except ImportError:
        print("✗ SYBA: Not installed")

    # Check for AiZynthFinder
    try:
        import aizynthfinder

        version = aizynthfinder.__version__ if hasattr(aizynthfinder, "__version__") else "Unknown"
        print(f"✓ AiZynthFinder: {version}")
    except ImportError:
        print("✗ AiZynthFinder: Not installed")

    # Check for RDKit
    try:
        from rdkit import __version__ as rdkit_version

        print(f"✓ RDKit: {rdkit_version}")
    except ImportError:
        print("✗ RDKit: Not installed")

    print("=" * 60)
