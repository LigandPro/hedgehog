"""Tests for docking/utils.py."""

import pandas as pd

from hedgehog.stages.docking.utils import (
    _aggregate_docking_results,
    _find_latest_input_source,
    _prepare_ligands_dataframe,
    _split_sdf_to_molecules,
)


class TestFindLatestInputSource:
    """Tests for _find_latest_input_source function."""

    def test_synthesis_output_new_structure(self, tmp_path):
        """Should find synthesis output in new structure."""
        synthesis_dir = tmp_path / "stages" / "04_synthesis"
        synthesis_dir.mkdir(parents=True)
        (synthesis_dir / "filtered_molecules.csv").write_text("smiles\nCCO")

        result = _find_latest_input_source(tmp_path)
        assert result is not None
        assert "synthesis" in str(result).lower()

    def test_structural_filters_output(self, tmp_path):
        """Should find structural filters output."""
        sf_dir = tmp_path / "stages" / "03_structural_filters_post"
        sf_dir.mkdir(parents=True)
        (sf_dir / "filtered_molecules.csv").write_text("smiles\nCCO")

        result = _find_latest_input_source(tmp_path)
        assert result is not None
        assert "structural_filters" in str(result).lower()

    def test_descriptors_output(self, tmp_path):
        """Should find descriptors output."""
        desc_dir = tmp_path / "stages" / "01_descriptors_initial" / "filtered"
        desc_dir.mkdir(parents=True)
        (desc_dir / "filtered_molecules.csv").write_text("smiles\nCCO")

        result = _find_latest_input_source(tmp_path)
        assert result is not None
        assert "descriptors" in str(result).lower()

    def test_sampled_molecules_input(self, tmp_path):
        """Should find sampled molecules in input directory."""
        input_dir = tmp_path / "input"
        input_dir.mkdir(parents=True)
        (input_dir / "sampled_molecules.csv").write_text("smiles\nCCO")

        result = _find_latest_input_source(tmp_path)
        assert result is not None
        assert "sampled_molecules" in str(result)

    def test_priority_order(self, tmp_path):
        """Should prioritize synthesis over descriptors."""
        # Create both synthesis and descriptors outputs
        synthesis_dir = tmp_path / "stages" / "04_synthesis"
        synthesis_dir.mkdir(parents=True)
        (synthesis_dir / "filtered_molecules.csv").write_text("smiles\nCCO")

        desc_dir = tmp_path / "stages" / "01_descriptors_initial" / "filtered"
        desc_dir.mkdir(parents=True)
        (desc_dir / "filtered_molecules.csv").write_text("smiles\nCC")

        result = _find_latest_input_source(tmp_path)
        assert result is not None
        # Synthesis should have higher priority
        assert "synthesis" in str(result).lower()

    def test_no_input_found(self, tmp_path):
        """Should return None when no input found."""
        result = _find_latest_input_source(tmp_path)
        assert result is None

    def test_legacy_structure(self, tmp_path):
        """Should find legacy flat structure files."""
        (tmp_path / "Synthesis").mkdir()
        (tmp_path / "Synthesis" / "passSynthesisSMILES.csv").write_text("smiles\nCCO")

        result = _find_latest_input_source(tmp_path)
        assert result is not None


class TestPrepareLigandsDataframe:
    """Tests for _prepare_ligands_dataframe function."""

    def test_valid_smiles(self, tmp_path):
        """All valid SMILES should be written."""
        df = pd.DataFrame(
            {
                "smiles": ["c1ccccc1", "CCO", "CC"],
                "model_name": ["test", "test", "test"],
                "mol_idx": ["t-0", "t-1", "t-2"],
            }
        )
        output_csv = tmp_path / "ligands.csv"
        stats = _prepare_ligands_dataframe(df, output_csv)

        assert stats["written"] == 3
        assert stats["skipped"] == 0
        assert stats["total"] == 3
        assert output_csv.exists()

    def test_invalid_smiles_skipped(self, tmp_path):
        """Invalid SMILES should be skipped."""
        df = pd.DataFrame(
            {
                "smiles": ["c1ccccc1", "invalid_smiles", "CCO"],
                "model_name": ["test", "test", "test"],
                "mol_idx": ["t-0", "t-1", "t-2"],
            }
        )
        output_csv = tmp_path / "ligands.csv"
        stats = _prepare_ligands_dataframe(df, output_csv)

        assert stats["written"] == 2
        assert stats["skipped"] == 1
        assert stats["total"] == 3

    def test_all_invalid_smiles(self, tmp_path):
        """All invalid SMILES - should write empty file."""
        df = pd.DataFrame(
            {
                "smiles": ["invalid1", "invalid2"],
                "model_name": ["test", "test"],
                "mol_idx": ["t-0", "t-1"],
            }
        )
        output_csv = tmp_path / "ligands.csv"
        stats = _prepare_ligands_dataframe(df, output_csv)

        assert stats["written"] == 0
        assert stats["skipped"] == 2

    def test_output_csv_columns(self, tmp_path):
        """Output CSV should have correct columns."""
        df = pd.DataFrame(
            {
                "smiles": ["c1ccccc1"],
                "model_name": ["test"],
                "mol_idx": ["t-0"],
            }
        )
        output_csv = tmp_path / "ligands.csv"
        _prepare_ligands_dataframe(df, output_csv)

        result = pd.read_csv(output_csv)
        assert "smiles" in result.columns
        assert "name" in result.columns
        assert "model_name" in result.columns
        assert "mol_idx" in result.columns

    def test_creates_parent_directories(self, tmp_path):
        """Should create parent directories if they don't exist."""
        df = pd.DataFrame(
            {
                "smiles": ["c1ccccc1"],
                "model_name": ["test"],
                "mol_idx": ["t-0"],
            }
        )
        output_csv = tmp_path / "nested" / "deep" / "ligands.csv"
        stats = _prepare_ligands_dataframe(df, output_csv)

        assert stats["written"] == 1
        assert output_csv.exists()

    def test_empty_dataframe(self, tmp_path):
        """Empty dataframe should result in empty output."""
        df = pd.DataFrame(
            {
                "smiles": [],
                "model_name": [],
                "mol_idx": [],
            }
        )
        output_csv = tmp_path / "ligands.csv"
        stats = _prepare_ligands_dataframe(df, output_csv)

        assert stats["written"] == 0
        assert stats["skipped"] == 0
        assert stats["total"] == 0


class TestFindLatestInputSourcePriority:
    """Additional tests for input source priority."""

    def test_struct_filters_over_descriptors(self, tmp_path):
        """Structural filters should be preferred over descriptors."""
        # Create both outputs
        sf_dir = tmp_path / "stages" / "03_structural_filters_post"
        sf_dir.mkdir(parents=True)
        (sf_dir / "filtered_molecules.csv").write_text("smiles\nCCO")

        desc_dir = tmp_path / "stages" / "01_descriptors_initial" / "filtered"
        desc_dir.mkdir(parents=True)
        (desc_dir / "filtered_molecules.csv").write_text("smiles\nCC")

        result = _find_latest_input_source(tmp_path)
        assert result is not None
        assert "structural_filters" in str(result).lower()

    def test_sampled_molecules_fallback(self, tmp_path):
        """Should fall back to sampled_molecules if no stage outputs."""
        input_dir = tmp_path / "input"
        input_dir.mkdir()
        (input_dir / "sampled_molecules.csv").write_text("smiles\nCCO")

        result = _find_latest_input_source(tmp_path)
        assert result is not None
        assert "sampled_molecules" in str(result)

    def test_legacy_descriptors_structure(self, tmp_path):
        """Should find legacy Descriptors directory."""
        legacy_dir = tmp_path / "Descriptors"
        legacy_dir.mkdir()
        (legacy_dir / "passDescriptorsSMILES.csv").write_text("smiles\nCCO")

        result = _find_latest_input_source(tmp_path)
        assert result is not None


class TestLigandNaming:
    """Tests for ligand name generation."""

    def test_name_column_format(self, tmp_path):
        """Name column should combine mol_idx and model_name."""
        df = pd.DataFrame(
            {
                "smiles": ["c1ccccc1", "CCO"],
                "model_name": ["model_a", "model_b"],
                "mol_idx": ["LP-0001-00001", "LP-0002-00001"],
            }
        )
        output_csv = tmp_path / "ligands.csv"
        _prepare_ligands_dataframe(df, output_csv)

        result = pd.read_csv(output_csv)
        assert "name" in result.columns
        # Names should be unique and contain mol_idx info
        assert len(result["name"].unique()) == 2

    def test_preserves_all_columns(self, tmp_path):
        """Should preserve identity columns in output."""
        df = pd.DataFrame(
            {
                "smiles": ["c1ccccc1"],
                "model_name": ["test"],
                "mol_idx": ["t-0"],
            }
        )
        output_csv = tmp_path / "ligands.csv"
        _prepare_ligands_dataframe(df, output_csv)

        result = pd.read_csv(output_csv)
        assert "smiles" in result.columns
        assert "model_name" in result.columns
        assert "mol_idx" in result.columns


class TestPerMoleculeArchitecture:
    """Tests for per-molecule SDF splitting and aggregation."""

    def test_split_sdf_to_molecules(self, tmp_path):
        """Should split multi-molecule SDF into individual files."""
        from rdkit import Chem

        # Create a multi-molecule SDF
        sdf_path = tmp_path / "multi.sdf"
        writer = Chem.SDWriter(str(sdf_path))

        mol1 = Chem.MolFromSmiles("CCO")
        mol1.SetProp("_Name", "ethanol")
        writer.write(mol1)

        mol2 = Chem.MolFromSmiles("c1ccccc1")
        mol2.SetProp("_Name", "benzene")
        writer.write(mol2)

        mol3 = Chem.MolFromSmiles("CC")
        mol3.SetProp("_Name", "ethane")
        writer.write(mol3)

        writer.close()

        # Split the SDF
        molecules_dir = tmp_path / "molecules"
        result = _split_sdf_to_molecules(sdf_path, molecules_dir)

        assert len(result) == 3
        assert molecules_dir.exists()
        assert (molecules_dir / "ethanol.sdf").exists()
        assert (molecules_dir / "benzene.sdf").exists()
        assert (molecules_dir / "ethane.sdf").exists()

    def test_split_sdf_with_unnamed_molecules(self, tmp_path):
        """Should handle molecules without names."""
        from rdkit import Chem

        sdf_path = tmp_path / "unnamed.sdf"
        writer = Chem.SDWriter(str(sdf_path))

        mol1 = Chem.MolFromSmiles("CCO")
        # No name set
        writer.write(mol1)

        mol2 = Chem.MolFromSmiles("CC")
        writer.write(mol2)

        writer.close()

        molecules_dir = tmp_path / "molecules"
        result = _split_sdf_to_molecules(sdf_path, molecules_dir)

        assert len(result) == 2
        # Should use default naming
        assert any("mol_" in str(f) for f in result)

    def test_aggregate_docking_results(self, tmp_path):
        """Should aggregate per-molecule results into single SDF."""
        from rdkit import Chem

        # Create per-molecule result files
        results_dir = tmp_path / "results"
        results_dir.mkdir()

        mol1 = Chem.MolFromSmiles("CCO")
        mol1.SetProp("_Name", "mol1")
        mol1.SetProp("score", "-5.2")
        writer1 = Chem.SDWriter(str(results_dir / "mol1_out.sdf"))
        writer1.write(mol1)
        writer1.close()

        mol2 = Chem.MolFromSmiles("c1ccccc1")
        mol2.SetProp("_Name", "mol2")
        mol2.SetProp("score", "-6.1")
        writer2 = Chem.SDWriter(str(results_dir / "mol2_out.sdf"))
        writer2.write(mol2)
        writer2.close()

        # Aggregate results
        output_sdf = tmp_path / "aggregated.sdf"
        count = _aggregate_docking_results(results_dir, output_sdf)

        assert count == 2
        assert output_sdf.exists()

        # Verify aggregated content
        suppl = Chem.SDMolSupplier(str(output_sdf))
        mols = [m for m in suppl if m is not None]
        assert len(mols) == 2

    def test_aggregate_empty_results_dir(self, tmp_path):
        """Should handle empty results directory gracefully."""
        results_dir = tmp_path / "empty_results"
        results_dir.mkdir()

        output_sdf = tmp_path / "output.sdf"
        count = _aggregate_docking_results(results_dir, output_sdf)

        assert count == 0

    def test_aggregate_with_invalid_files(self, tmp_path):
        """Should skip invalid SDF files during aggregation."""
        from rdkit import Chem

        results_dir = tmp_path / "results"
        results_dir.mkdir()

        # Valid molecule
        mol1 = Chem.MolFromSmiles("CCO")
        mol1.SetProp("_Name", "valid")
        writer = Chem.SDWriter(str(results_dir / "valid_out.sdf"))
        writer.write(mol1)
        writer.close()

        # Invalid file (empty)
        (results_dir / "invalid_out.sdf").write_text("")

        output_sdf = tmp_path / "output.sdf"
        count = _aggregate_docking_results(results_dir, output_sdf)

        # Should have 1 valid molecule
        assert count == 1

    def test_split_sdf_sanitizes_names(self, tmp_path):
        """Should sanitize molecule names with special characters."""
        from rdkit import Chem

        sdf_path = tmp_path / "special.sdf"
        writer = Chem.SDWriter(str(sdf_path))

        mol = Chem.MolFromSmiles("CCO")
        mol.SetProp("_Name", "mol/with:special*chars")
        writer.write(mol)
        writer.close()

        molecules_dir = tmp_path / "molecules"
        result = _split_sdf_to_molecules(sdf_path, molecules_dir)

        assert len(result) == 1
        # Name should be sanitized (no special chars)
        filename = result[0].name
        assert "/" not in filename
        assert ":" not in filename
        assert "*" not in filename
