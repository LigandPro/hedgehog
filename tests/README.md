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

### Install Test Dependencies

```bash
uv pip install -e ".[test]"
```

### Run All Tests

```bash
pytest
```

### Run Specific Test Categories

```bash
# Unit tests only
pytest -m unit

# Integration tests only
pytest -m integration

# End-to-end tests only
pytest -m e2e

# Exclude slow tests
pytest -m "not slow"
```

### Run Specific Test Files

```bash
pytest tests/test_data_prep.py
pytest tests/test_e2e_pipeline.py
```

### Run with Coverage

```bash
pytest --cov=src/hedge --cov-report=html
```

This will generate a coverage report in `htmlcov/index.html`.

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
