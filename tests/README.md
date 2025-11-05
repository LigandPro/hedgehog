# HEDGE Tests

Comprehensive test suite for the HEDGE (Hierarchical Evaluation of Drug GEnerators) pipeline.

## Test Structure

```
tests/
├── conftest.py              # Shared fixtures and configuration
├── fixtures/                # Test data and configurations
│   ├── test_molecules.csv
│   ├── test_config.yml
│   └── test_config_*.yml
├── test_e2e_pipeline.py     # End-to-end integration tests
├── test_data_prep.py        # Data preparation unit tests
├── test_mol_index.py        # Molecular indexing unit tests
├── test_pipeline.py         # Pipeline orchestration tests
└── test_cli.py              # CLI and preprocessing tests
```

## Running Tests

### Prerequisites

HEDGE uses `uv` with system-site-packages enabled venv. This allows tests to access all pre-installed packages without re-downloading.

**IMPORTANT:** Always use `UV_NO_SYNC=1` with `uv run` to prevent dependency re-installation:

```bash
# Run tests with uv (no package sync)
UV_NO_SYNC=1 uv run pytest
```

### Run All Tests

```bash
# Run all tests
UV_NO_SYNC=1 uv run pytest
```

### Run Specific Test Categories

```bash
# Unit tests only (fast, ~6 seconds, 47 tests)
UV_NO_SYNC=1 uv run pytest -m unit

# Integration tests only
UV_NO_SYNC=1 uv run pytest -m integration

# End-to-end tests only
UV_NO_SYNC=1 uv run pytest -m e2e

# External dependencies tests (Lilly, SYBA, AiZynthFinder)
UV_NO_SYNC=1 uv run pytest -m external

# Exclude slow tests
UV_NO_SYNC=1 uv run pytest -m "not slow"

# Exclude external dependencies tests (fastest)
UV_NO_SYNC=1 uv run pytest -m "not external"
```

### Run Specific Test Files

```bash
UV_NO_SYNC=1 uv run pytest tests/test_data_prep.py
UV_NO_SYNC=1 uv run pytest tests/test_e2e_pipeline.py
UV_NO_SYNC=1 uv run pytest tests/test_external_deps.py
```

### Run with Coverage

```bash
UV_NO_SYNC=1 uv run pytest --cov=src/hedge --cov-report=html
```

This will generate a coverage report in `htmlcov/index.html`.

### Quick Test Run (Recommended)

```bash
# Fast unit tests only, skip external dependencies
UV_NO_SYNC=1 uv run pytest -m "unit and not external" -v

# Output: 47 passed in ~6 seconds ✓
```

### Check External Dependencies Status

```bash
# This will show which external tools are installed
UV_NO_SYNC=1 uv run pytest tests/test_external_deps.py::test_conda_packages_info -v -s
```

**Example Output:**
```
============================================================
EXTERNAL DEPENDENCIES STATUS
============================================================
Python: 3.11.14
Environment: /home/user/hedge/.venv
✓ RDKit: 2025.09.1
✓ medchem: Available
✓ datamol: Available
✗ Lilly: Not found in PATH (optional - conda package)
✗ SYBA: Not installed (optional - conda package)
✗ AiZynthFinder: Not installed (optional)
============================================================
```

**Note:** External dependencies (Lilly, SYBA, AiZynthFinder) are optional. Core 47 unit tests will run without them.

## Test Categories

### Unit Tests (`@pytest.mark.unit`)
Fast, isolated tests for individual functions and classes:
- Data preparation utilities
- Molecular index assignment
- Pipeline components
- CLI helpers

### Integration Tests (`@pytest.mark.integration`)
Tests for component interactions:
- Pipeline stage coordination
- Configuration loading
- File I/O operations

### End-to-End Tests (`@pytest.mark.e2e`)
Full pipeline execution tests:
- `hedge run` command execution
- Complete workflow validation
- Output file generation

### Slow Tests (`@pytest.mark.slow`)
Tests that may take longer to execute:
- Full pipeline runs
- Large dataset processing

### External Dependencies Tests (`@pytest.mark.external`)
Tests for external tools and libraries:
- Lilly MedChem Rules (conda package)
- SYBA (conda package)
- AiZynthFinder (retrosynthesis)
- RDKit integration
- These tests require conda environment to be activated

## Test Fixtures

### Common Fixtures (in `conftest.py`)

- `fixtures_dir` - Path to test fixtures directory
- `test_molecules_csv` - Path to test molecules CSV
- `test_config_path` - Path to test configuration
- `temp_output_dir` - Temporary directory for test outputs
- `sample_smiles_data` - Sample SMILES DataFrame
- `sample_multi_model_data` - Multi-model sample data
- `sample_csv_file` - Temporary CSV file with sample data
- `sample_config_dict` - Sample configuration dictionary

## Writing New Tests

### Test Naming Convention

- Test files: `test_*.py`
- Test functions: `test_*`
- Test classes: `Test*`

### Example Test

```python
import pytest

@pytest.mark.unit
def test_my_function(sample_smiles_data, temp_output_dir):
    """Test description."""
    result = my_function(sample_smiles_data, temp_output_dir)
    assert result is not None
    assert len(result) > 0
```

### Marking Tests

```python
@pytest.mark.unit          # Unit test
@pytest.mark.integration   # Integration test
@pytest.mark.e2e           # End-to-end test
@pytest.mark.slow          # Slow-running test
```

## Coverage Goals

- **Overall**: 70%+
- **Critical modules** (data_prep, pipeline, mol_index): 80%+
- **Unit tests**: Fast (<1s per test)
- **E2E tests**: May take 30-60s

## Continuous Integration

Tests are designed to run in CI/CD pipelines:
- All tests should be deterministic
- Use fixtures for file I/O
- Mock external tool dependencies
- Clean up temporary files

## Troubleshooting

### Tests Fail with Import Errors

Ensure hedge is installed:
```bash
uv pip install -e .
```

### Tests Fail with Missing Fixtures

Check that fixture files exist:
```bash
ls tests/fixtures/
```

### Coverage Report Not Generated

Install coverage dependency:
```bash
uv pip install pytest-cov
```
