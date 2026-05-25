"""Tests for large-dataset streaming helpers."""

import json

import pandas as pd

import hedgehog.struct_filters.large as struct_large
import hedgehog.synthesis.large as synthesis_large
from hedgehog.large_dataset import (
    ShardedCsvWriter,
    StreamingMolIdxAssigner,
    count_part_rows,
    iter_csv_parts,
    iter_input_chunks,
    materialize_csv_if_small,
    parts_dir_for_csv,
)
from hedgehog.pipeline import (
    DIR_DESCRIPTORS_INITIAL,
    DIR_STRUCT_FILTERS_POST,
    DIR_SYNTHESIS,
    FILE_FILTERED_MOLECULES,
    STAGE_STRUCT_FILTERS,
    DataChecker,
    MolecularAnalysisPipeline,
    PipelineStageRunner,
)
from hedgehog.struct_filters.large import _struct_input_chunks


def test_sharded_writer_materializes_small_csv(tmp_path):
    output_csv = tmp_path / "filtered_molecules.csv"
    writer = ShardedCsvWriter(
        parts_dir_for_csv(output_csv), {"large_dataset_output_format": "csv.gz"}
    )

    writer.write(
        pd.DataFrame(
            {
                "smiles": ["CCO", "CCN"],
                "model_name": ["model_a", "model_a"],
                "mol_idx": ["LP-0001-00001", "LP-0001-00002"],
            }
        )
    )
    writer.write(
        pd.DataFrame(
            {
                "smiles": ["CCC"],
                "model_name": ["model_b"],
                "mol_idx": ["LP-0002-00001"],
            }
        )
    )

    total = materialize_csv_if_small(writer.parts_dir, output_csv, row_limit=10)

    assert total == 3
    assert count_part_rows(writer.parts_dir) == 3
    assert pd.read_csv(output_csv)["smiles"].tolist() == ["CCO", "CCN", "CCC"]


def test_materialize_empty_parts_writes_header_only_csv(tmp_path):
    output_csv = tmp_path / "filtered_molecules.csv"
    parts_dir = parts_dir_for_csv(output_csv)
    parts_dir.mkdir()

    total = materialize_csv_if_small(
        parts_dir,
        output_csv,
        row_limit=10,
        columns=["smiles", "model_name", "mol_idx"],
    )

    assert total == 0
    df = pd.read_csv(output_csv)
    assert df.empty
    assert df.columns.tolist() == ["smiles", "model_name", "mol_idx"]


def test_iter_csv_parts_skips_empty_files(tmp_path):
    empty_csv = tmp_path / "empty.csv"
    empty_csv.touch()

    assert list(iter_csv_parts(empty_csv)) == []


def test_iter_input_chunks_expands_globbed_csv_inputs(tmp_path):
    first = tmp_path / "a.csv"
    second = tmp_path / "b.csv"
    first.write_text("smiles\nCCO\n", encoding="utf-8")
    second.write_text("smiles\nCCN\n", encoding="utf-8")

    chunks = list(iter_input_chunks(str(tmp_path / "*.csv"), chunk_rows=1))

    assert [chunk["smiles"].tolist() for chunk in chunks] == [["CCO"], ["CCN"]]
    assert [chunk["model_name"].tolist() for chunk in chunks] == [["a"], ["b"]]


def test_iter_input_chunks_reads_headerless_txt_as_smi(tmp_path):
    input_txt = tmp_path / "mols.txt"
    input_txt.write_text("CCO model_a\nCCN model_b\n", encoding="utf-8")

    chunks = list(iter_input_chunks(input_txt, chunk_rows=10))

    assert len(chunks) == 1
    assert chunks[0][["smiles", "model_name"]].to_dict("records") == [
        {"smiles": "CCO", "model_name": "model_a"},
        {"smiles": "CCN", "model_name": "model_b"},
    ]


def test_streaming_mol_idx_assigner_persists_counters_across_chunks(tmp_path):
    assigner = StreamingMolIdxAssigner(tmp_path / "state.sqlite", tmp_path)
    try:
        first = assigner.assign(
            pd.DataFrame(
                {
                    "smiles": ["CCO", "CCN", "CCC"],
                    "model_name": ["model_a", "model_a", "model_b"],
                    "mol_idx": [pd.NA, pd.NA, pd.NA],
                }
            )
        )
        second = assigner.assign(
            pd.DataFrame(
                {
                    "smiles": ["CO", "CN"],
                    "model_name": ["model_a", "model_b"],
                    "mol_idx": [pd.NA, pd.NA],
                }
            )
        )
    finally:
        assigner.close()

    assert first["mol_idx"].tolist() == [
        "LP-0001-00001",
        "LP-0001-00002",
        "LP-0002-00001",
    ]
    assert second["mol_idx"].tolist() == ["LP-0001-00003", "LP-0002-00002"]

    model_map = json.loads((tmp_path / "configs" / "model_index_map.json").read_text())
    assert model_map == {"model_a": 1, "model_b": 2}


def test_data_checker_accepts_large_dataset_parts_output(tmp_path):
    output_csv = (
        tmp_path / DIR_DESCRIPTORS_INITIAL / "filtered" / FILE_FILTERED_MOLECULES
    )
    writer = ShardedCsvWriter(
        parts_dir_for_csv(output_csv), {"large_dataset_output_format": "csv.gz"}
    )
    writer.write(
        pd.DataFrame(
            {
                "smiles": ["CCO"],
                "model_name": ["model_a"],
                "mol_idx": ["LP-0001-00001"],
            }
        )
    )

    checker = DataChecker({"folder_to_save": str(tmp_path), "large_dataset_mode": True})

    assert checker.check_stage_data(DIR_DESCRIPTORS_INITIAL) is True
    assert checker.stage_has_molecules(DIR_DESCRIPTORS_INITIAL) is True


def test_large_final_output_does_not_use_stale_later_stage(tmp_path):
    struct_csv = tmp_path / DIR_STRUCT_FILTERS_POST / FILE_FILTERED_MOLECULES
    synth_csv = tmp_path / DIR_SYNTHESIS / FILE_FILTERED_MOLECULES
    struct_csv.parent.mkdir(parents=True)
    synth_csv.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "smiles": ["CCO"],
            "model_name": ["struct"],
            "mol_idx": ["LP-0001-00001"],
        }
    ).to_csv(struct_csv, index=False)
    pd.DataFrame(
        {
            "smiles": ["CCN"],
            "model_name": ["stale_synthesis"],
            "mol_idx": ["LP-0001-00002"],
        }
    ).to_csv(synth_csv, index=False)

    pipeline = object.__new__(MolecularAnalysisPipeline)
    pipeline.data_checker = DataChecker(
        {"folder_to_save": str(tmp_path), "large_dataset_mode": True}
    )

    assert pipeline._find_large_final_output(STAGE_STRUCT_FILTERS) == struct_csv
    assert pipeline._find_large_final_output(None) is None


def test_large_synthesis_validation_accepts_generated_input(tmp_path):
    input_csv = tmp_path / "mols.csv"
    input_csv.write_text("smiles\nCCO\n", encoding="utf-8")
    config = {
        "folder_to_save": str(tmp_path),
        "generated_mols_path": str(input_csv),
        "large_dataset_mode": True,
    }
    runner = PipelineStageRunner(config, DataChecker(config))

    assert runner._validate_synthesis_input() is True


def test_struct_filters_respect_empty_descriptor_output(tmp_path):
    descriptor_csv = (
        tmp_path / DIR_DESCRIPTORS_INITIAL / "filtered" / FILE_FILTERED_MOLECULES
    )
    descriptor_csv.parent.mkdir(parents=True)
    pd.DataFrame(columns=["smiles", "model_name", "mol_idx"]).to_csv(
        descriptor_csv, index=False
    )

    molprep_csv = tmp_path / "stages" / "01_mol_prep" / FILE_FILTERED_MOLECULES
    writer = ShardedCsvWriter(
        parts_dir_for_csv(molprep_csv), {"large_dataset_output_format": "csv.gz"}
    )
    writer.write(
        pd.DataFrame(
            {
                "smiles": ["CCO"],
                "model_name": ["model_a"],
                "mol_idx": ["LP-0001-00001"],
            }
        )
    )

    chunks = list(
        _struct_input_chunks(
            {"folder_to_save": str(tmp_path), "generated_mols_path": "unused.csv"},
            "stages/03_structural_filters_post",
            chunk_rows=10,
        )
    )

    assert len(chunks) == 1
    assert chunks[0].empty


def test_struct_filters_large_computes_without_filtering_outputs(tmp_path, monkeypatch):
    input_df = pd.DataFrame(
        {
            "smiles": ["CCO", "CCN"],
            "model_name": ["model_a", "model_a"],
            "mol_idx": ["LP-0001-00001", "LP-0001-00002"],
        }
    )

    monkeypatch.setattr(
        struct_large,
        "_struct_input_chunks",
        lambda *_args, **_kwargs: iter([input_df]),
    )
    monkeypatch.setattr(
        struct_large,
        "prepare_structfilters_input",
        lambda *_args, **_kwargs: {
            "mols": [object()],
            "smiles_model_mols": [("CCO", "model_a", object(), "LP-0001-00001")],
        },
    )
    monkeypatch.setattr(
        struct_large,
        "filter_function_applier",
        lambda _name: lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        struct_large,
        "process_prepared_payload",
        lambda *_args, **_kwargs: pd.DataFrame({"raw": [1, 2]}),
    )

    def _fake_basic_stats(*_args, **_kwargs):
        return pd.DataFrame(), pd.DataFrame(
            {
                "smiles": input_df["smiles"],
                "model_name": input_df["model_name"],
                "mol_idx": input_df["mol_idx"],
                "pass": [False, True],
            }
        )

    monkeypatch.setattr(struct_large, "get_basic_stats", _fake_basic_stats)

    out_dir = tmp_path / "stages" / "03_structural_filters_post"
    struct_large.run_large(
        {
            "folder_to_save": str(tmp_path),
            "generated_mols_path": "unused.csv",
            "large_dataset_filter_data": False,
            "large_dataset_enable_all_filters": True,
            "large_dataset_output_format": "csv.gz",
        },
        "stages/03_structural_filters_post",
        {"calculate_dummy": False, "parse_input_n_jobs": 1},
        out_dir,
    )

    filtered = pd.read_csv(out_dir / "filtered_molecules.csv")
    assert filtered["mol_idx"].tolist() == [
        "LP-0001-00001",
        "LP-0001-00002",
    ]
    assert (out_dir / "intermediate" / "dummy").exists()


def test_struct_filters_large_records_filter_failures_without_dropping(
    tmp_path, monkeypatch
):
    input_df = pd.DataFrame(
        {
            "smiles": ["CCO"],
            "model_name": ["model_a"],
            "mol_idx": ["LP-0001-00001"],
        }
    )

    monkeypatch.setattr(
        struct_large,
        "_struct_input_chunks",
        lambda *_args, **_kwargs: iter([input_df]),
    )
    monkeypatch.setattr(
        struct_large,
        "prepare_structfilters_input",
        lambda *_args, **_kwargs: {
            "mols": [object()],
            "smiles_model_mols": [("CCO", "model_a", object(), "LP-0001-00001")],
        },
    )

    def _failing_filter(*_args, **_kwargs):
        raise RuntimeError("missing optional backend")

    monkeypatch.setattr(
        struct_large,
        "filter_function_applier",
        lambda _name: _failing_filter,
    )

    out_dir = tmp_path / "stages" / "03_structural_filters_post"
    struct_large.run_large(
        {
            "folder_to_save": str(tmp_path),
            "generated_mols_path": "unused.csv",
            "large_dataset_filter_data": False,
            "large_dataset_enable_all_filters": True,
            "large_dataset_output_format": "csv.gz",
        },
        "stages/03_structural_filters_post",
        {"calculate_unavailable": False, "parse_input_n_jobs": 1},
        out_dir,
    )

    filtered = pd.read_csv(out_dir / "filtered_molecules.csv")
    failures = pd.read_csv(out_dir / "summary" / "filter_failures.tsv", sep="\t")

    assert filtered["mol_idx"].tolist() == ["LP-0001-00001"]
    assert failures[["filter", "error"]].to_dict("records") == [
        {"filter": "unavailable", "error": "missing optional backend"}
    ]


def test_synthesis_large_skips_retrosynthesis_and_keeps_score_failures(
    tmp_path, monkeypatch
):
    input_df = pd.DataFrame(
        {
            "smiles": ["CCO", "CCN"],
            "model_name": ["model_a", "model_a"],
            "mol_idx": ["LP-0001-00001", "LP-0001-00002"],
        }
    )

    monkeypatch.setattr(
        synthesis_large,
        "_synthesis_input_chunks",
        lambda *_args, **_kwargs: iter([input_df]),
    )

    def _fake_scores(df, *_args, **_kwargs):
        scored = df.copy()
        scored["sa_score"] = [2.0, 9.0]
        scored["ra_score"] = [0.9, 0.1]
        scored["syba_score"] = [1.0, -1.0]
        return scored

    monkeypatch.setattr(synthesis_large, "calculate_synthesis_scores", _fake_scores)

    out_dir = tmp_path / "stages" / "04_synthesis"
    synthesis_large.run_large(
        {
            "folder_to_save": str(tmp_path),
            "generated_mols_path": "unused.csv",
            "large_dataset_filter_data": False,
            "large_dataset_output_format": "csv.gz",
        },
        {
            "n_jobs": 1,
            "sa_score_min": 1,
            "sa_score_max": 4.5,
            "ra_score_min": 0.5,
            "ra_score_max": 1,
            "syba_score_min": 0,
            "syba_score_max": "inf",
        },
        out_dir,
    )

    filtered = pd.read_csv(out_dir / "filtered_molecules.csv")
    flags = pd.read_csv(out_dir / "pass_flags.csv")

    assert filtered["mol_idx"].tolist() == [
        "LP-0001-00001",
        "LP-0001-00002",
    ]
    assert flags["synthesis_score_pass"].tolist() == [True, False]
    assert not (out_dir / "retrosynthesis_results.json").exists()
