from __future__ import annotations

import csv
import glob
import gzip
import json
import shutil
import sqlite3
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pandas as pd

from hedgehog.utils.mol_index import (
    DEFAULT_MODEL_NAME,
    MODEL_INDEX_MAP_FILE,
    MODEL_NAME_COLUMN,
    MOL_IDX_COLUMN,
    MOL_IDX_FORMAT,
    RUN_CONFIGS_DIR,
)

LARGE_DATASET_MODE_KEY = "large_dataset_mode"
LARGE_DATASET_CHUNK_ROWS_KEY = "large_dataset_chunk_rows"
LARGE_DATASET_SINGLE_CSV_LIMIT_KEY = "large_dataset_single_csv_limit"
LARGE_DATASET_OUTPUT_FORMAT_KEY = "large_dataset_output_format"
LARGE_DATASET_FILTER_DATA_KEY = "large_dataset_filter_data"
LARGE_DATASET_ENABLE_ALL_FILTERS_KEY = "large_dataset_enable_all_filters"

DEFAULT_CHUNK_ROWS = 250_000
DEFAULT_SINGLE_CSV_LIMIT = 1_000_000
DEFAULT_OUTPUT_FORMAT = "csv.gz"
DEFAULT_FILTER_DATA = False
DEFAULT_ENABLE_ALL_FILTERS = True

SMILES_COLUMN = "smiles"
NAME_COLUMN = "name"
SUPPORTED_TABLE_SUFFIXES = {".csv", ".tsv"}
SMI_SUFFIXES = {".smi", ".ismi", ".cmi", ".smiles", ".txt"}
SUPPORTED_INPUT_SUFFIXES = SUPPORTED_TABLE_SUFFIXES | SMI_SUFFIXES


def apply_large_dataset_defaults(config: dict[str, Any]) -> None:
    """Populate large-dataset defaults without overwriting explicit config values."""
    config.setdefault(LARGE_DATASET_CHUNK_ROWS_KEY, DEFAULT_CHUNK_ROWS)
    config.setdefault(LARGE_DATASET_SINGLE_CSV_LIMIT_KEY, DEFAULT_SINGLE_CSV_LIMIT)
    config.setdefault(LARGE_DATASET_OUTPUT_FORMAT_KEY, DEFAULT_OUTPUT_FORMAT)
    config.setdefault(LARGE_DATASET_FILTER_DATA_KEY, DEFAULT_FILTER_DATA)
    config.setdefault(LARGE_DATASET_ENABLE_ALL_FILTERS_KEY, DEFAULT_ENABLE_ALL_FILTERS)


def is_large_dataset_mode(config: dict[str, Any] | None) -> bool:
    return bool((config or {}).get(LARGE_DATASET_MODE_KEY, False))


def _get_int_config(
    config: dict[str, Any] | None, key: str, default: int, *, allow_zero: bool
) -> int:
    raw = (config or {}).get(key, default)
    try:
        value = int(raw)
    except (TypeError, ValueError):
        return default
    if allow_zero:
        return value if value >= 0 else default
    return value if value > 0 else default


def get_chunk_rows(config: dict[str, Any] | None) -> int:
    return _get_int_config(
        config, LARGE_DATASET_CHUNK_ROWS_KEY, DEFAULT_CHUNK_ROWS, allow_zero=False
    )


def get_single_csv_limit(config: dict[str, Any] | None) -> int:
    return _get_int_config(
        config,
        LARGE_DATASET_SINGLE_CSV_LIMIT_KEY,
        DEFAULT_SINGLE_CSV_LIMIT,
        allow_zero=True,
    )


def get_output_format(config: dict[str, Any] | None) -> str:
    value = str((config or {}).get(LARGE_DATASET_OUTPUT_FORMAT_KEY, "")).strip()
    return value or DEFAULT_OUTPUT_FORMAT


def should_filter_large_outputs(config: dict[str, Any] | None) -> bool:
    return bool((config or {}).get(LARGE_DATASET_FILTER_DATA_KEY, DEFAULT_FILTER_DATA))


def should_enable_all_large_filters(config: dict[str, Any] | None) -> bool:
    return bool(
        (config or {}).get(
            LARGE_DATASET_ENABLE_ALL_FILTERS_KEY, DEFAULT_ENABLE_ALL_FILTERS
        )
    )


def parts_dir_for_csv(csv_path: Path) -> Path:
    """Return the sharded directory matching a CSV output path."""
    name = csv_path.name
    if name.endswith(".csv"):
        name = name[: -len(".csv")]
    return csv_path.with_name(f"{name}.parts")


def part_suffix(config: dict[str, Any] | None) -> str:
    fmt = get_output_format(config).lower()
    if fmt == "csv":
        return ".csv"
    return ".csv.gz"


def list_part_files(parts_dir: Path) -> list[Path]:
    if not parts_dir.exists() or not parts_dir.is_dir():
        return []
    return sorted(
        [
            p
            for p in parts_dir.iterdir()
            if p.is_file() and (p.name.endswith(".csv") or p.name.endswith(".csv.gz"))
        ]
    )


def _open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, mode="rt", encoding="utf-8", errors="ignore", newline="")
    return path.open(encoding="utf-8", errors="ignore", newline="")


def count_csv_rows(path: Path) -> int | None:
    if not path.exists() or not path.is_file():
        return None
    try:
        with _open_text(path) as handle:
            reader = csv.reader(handle)
            try:
                next(reader)
            except StopIteration:
                return 0
            return sum(1 for row in reader if any(str(cell).strip() for cell in row))
    except OSError:
        return None


def count_part_rows(parts_dir: Path) -> int | None:
    files = list_part_files(parts_dir)
    if not files:
        return None
    total = 0
    for path in files:
        counted = count_csv_rows(path)
        if counted is None:
            return None
        total += counted
    return total


def _resolve_csv_or_parts(csv_path: Path, *, allow_empty_csv: bool) -> Path | None:
    if csv_path.exists():
        if allow_empty_csv or csv_path.stat().st_size > 0:
            return csv_path
    parts_dir = parts_dir_for_csv(csv_path)
    if allow_empty_csv and parts_dir.exists():
        return parts_dir
    if not allow_empty_csv and list_part_files(parts_dir):
        return parts_dir
    return None


def output_path_or_parts(csv_path: Path) -> Path | None:
    return _resolve_csv_or_parts(csv_path, allow_empty_csv=False)


def stage_output_or_parts(csv_path: Path) -> Path | None:
    """Return a stage output even when it is an intentional empty CSV."""
    return _resolve_csv_or_parts(csv_path, allow_empty_csv=True)


def iter_csv_parts(path_or_parts: Path, chunk_rows: int | None = None) -> Iterator[pd.DataFrame]:
    """Yield DataFrames from one CSV/CSV.GZ file or a .parts directory."""
    if path_or_parts.is_dir():
        for part in list_part_files(path_or_parts):
            try:
                yield pd.read_csv(part)
            except pd.errors.EmptyDataError:
                continue
        return

    kwargs: dict[str, Any] = {}
    if chunk_rows is not None and chunk_rows > 0:
        kwargs["chunksize"] = chunk_rows
    try:
        reader = pd.read_csv(path_or_parts, **kwargs)
    except pd.errors.EmptyDataError:
        return
    if isinstance(reader, pd.DataFrame):
        yield reader
    else:
        yield from reader


def normalize_molecule_chunk(df: pd.DataFrame, source_path: Path) -> pd.DataFrame:
    """Normalize one molecule table chunk to smiles/model_name/mol_idx columns."""
    if df is None or len(df) == 0:
        return pd.DataFrame(columns=[SMILES_COLUMN, MODEL_NAME_COLUMN, MOL_IDX_COLUMN])

    out = df.copy()
    lower_cols = {str(c).lower(): c for c in out.columns}
    smiles_col = lower_cols.get(SMILES_COLUMN)
    if smiles_col is None:
        if len(out.columns) == 1:
            out = out.rename(columns={out.columns[0]: SMILES_COLUMN})
        else:
            msg = f"Input file {source_path} is missing required 'smiles' column."
            raise ValueError(msg)
    elif smiles_col != SMILES_COLUMN:
        out = out.rename(columns={smiles_col: SMILES_COLUMN})

    model_col = lower_cols.get(MODEL_NAME_COLUMN) or lower_cols.get(NAME_COLUMN)
    if model_col is not None and model_col != MODEL_NAME_COLUMN:
        out = out.rename(columns={model_col: MODEL_NAME_COLUMN})
    elif MODEL_NAME_COLUMN not in out.columns:
        out[MODEL_NAME_COLUMN] = source_path.stem

    if MOL_IDX_COLUMN not in out.columns:
        out[MOL_IDX_COLUMN] = pd.NA

    for col in (SMILES_COLUMN, MODEL_NAME_COLUMN, MOL_IDX_COLUMN):
        out[col] = out[col].astype("string").str.strip()
        out[col] = out[col].mask(out[col].eq(""), pd.NA)

    if out[SMILES_COLUMN].isna().any():
        invalid_rows = int(out[SMILES_COLUMN].isna().sum())
        msg = (
            f"Input file {source_path} contains {invalid_rows} row(s) with empty "
            "'smiles' values."
        )
        raise ValueError(msg)

    out[MODEL_NAME_COLUMN] = out[MODEL_NAME_COLUMN].fillna(DEFAULT_MODEL_NAME)
    keep_cols = [SMILES_COLUMN, MODEL_NAME_COLUMN, MOL_IDX_COLUMN]
    return out[keep_cols].copy()


def _expand_input_paths(input_path: str | Path) -> list[Path]:
    raw_path = str(input_path)
    if any(ch in raw_path for ch in "*?[]"):
        paths = [Path(p) for p in sorted(glob.glob(raw_path))]
    else:
        path = Path(input_path)
        if path.is_dir():
            paths = sorted(
                p
                for p in path.iterdir()
                if p.is_file() and p.suffix.lower() in SUPPORTED_INPUT_SUFFIXES
            )
        else:
            paths = [path]

    if not paths:
        raise FileNotFoundError(raw_path)
    return paths


def _iter_smi_chunks(path: Path, chunk_rows: int) -> Iterator[pd.DataFrame]:
    batch: list[dict[str, str]] = []
    with path.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if parts[0].lower() == SMILES_COLUMN:
                continue
            batch.append(
                {
                    SMILES_COLUMN: parts[0],
                    MODEL_NAME_COLUMN: parts[1] if len(parts) > 1 else path.stem,
                }
            )
            if len(batch) >= chunk_rows:
                yield normalize_molecule_chunk(pd.DataFrame(batch), path)
                batch = []
    if batch:
        yield normalize_molecule_chunk(pd.DataFrame(batch), path)


def _iter_table_chunks(path: Path, chunk_rows: int) -> Iterator[pd.DataFrame]:
    suffix = path.suffix.lower()
    sep = "\t" if suffix == ".tsv" else ","
    reader = pd.read_csv(path, sep=sep, chunksize=chunk_rows)
    for chunk in reader:
        yield normalize_molecule_chunk(chunk, path)


def iter_input_chunks(input_path: str | Path, chunk_rows: int) -> Iterator[pd.DataFrame]:
    """Stream molecule input chunks from CSV/TSV or headerless SMI-style files."""
    normalized_chunk_rows = max(int(chunk_rows), 1)
    for path in _expand_input_paths(input_path):
        suffix = path.suffix.lower()
        if suffix in SMI_SUFFIXES:
            yield from _iter_smi_chunks(path, normalized_chunk_rows)
            continue
        yield from _iter_table_chunks(path, normalized_chunk_rows)


class StreamingMolIdxAssigner:
    """Assign stable mol_idx values incrementally without loading all rows."""

    def __init__(self, state_db: Path, run_base: Path):
        self.state_db = Path(state_db)
        self.run_base = Path(run_base)
        self.state_db.parent.mkdir(parents=True, exist_ok=True)
        self.conn = sqlite3.connect(self.state_db)
        self.conn.execute(
            """
            CREATE TABLE IF NOT EXISTS model_state (
                model_name TEXT PRIMARY KEY,
                model_number INTEGER NOT NULL UNIQUE,
                next_counter INTEGER NOT NULL
            )
            """
        )
        self.conn.commit()

    def close(self) -> None:
        self.conn.close()

    def _get_or_create_model(self, model_name: str) -> tuple[int, int]:
        row = self.conn.execute(
            "SELECT model_number, next_counter FROM model_state WHERE model_name = ?",
            (model_name,),
        ).fetchone()
        if row is not None:
            return int(row[0]), int(row[1])

        max_row = self.conn.execute("SELECT MAX(model_number) FROM model_state").fetchone()
        next_model = int(max_row[0] or 0) + 1
        self.conn.execute(
            "INSERT INTO model_state(model_name, model_number, next_counter) VALUES (?, ?, ?)",
            (model_name, next_model, 1),
        )
        self.conn.commit()
        self._persist_model_map()
        return next_model, 1

    def _persist_model_map(self) -> None:
        rows = self.conn.execute(
            "SELECT model_name, model_number FROM model_state ORDER BY model_number"
        ).fetchall()
        model_map = {str(model): int(number) for model, number in rows}
        dest_dir = self.run_base / RUN_CONFIGS_DIR
        dest_dir.mkdir(parents=True, exist_ok=True)
        (dest_dir / MODEL_INDEX_MAP_FILE).write_text(
            json.dumps(model_map, indent=2, sort_keys=True),
            encoding="utf-8",
        )

    def assign(self, df: pd.DataFrame) -> pd.DataFrame:
        out = df.copy()
        if MOL_IDX_COLUMN not in out.columns:
            out[MOL_IDX_COLUMN] = pd.NA

        missing_mask = out[MOL_IDX_COLUMN].isna() | out[MOL_IDX_COLUMN].astype(
            "string"
        ).str.strip().eq("")
        if not missing_mask.any():
            return out

        for model_name, idxs in out.loc[missing_mask].groupby(MODEL_NAME_COLUMN).groups.items():
            model_number, next_counter = self._get_or_create_model(str(model_name))
            ordered_indices = list(idxs)
            values = [
                MOL_IDX_FORMAT.format(model_number, next_counter + offset)
                for offset in range(len(ordered_indices))
            ]
            out.loc[ordered_indices, MOL_IDX_COLUMN] = values
            self.conn.execute(
                "UPDATE model_state SET next_counter = ? WHERE model_name = ?",
                (next_counter + len(ordered_indices), str(model_name)),
            )
        self.conn.commit()
        return out


class ShardedCsvWriter:
    """Write row-level tables as part-NNNNNN CSV shards."""

    def __init__(self, parts_dir: Path, config: dict[str, Any] | None = None):
        self.parts_dir = Path(parts_dir)
        self.parts_dir.mkdir(parents=True, exist_ok=True)
        self.suffix = part_suffix(config)
        self.part_index = len(list_part_files(self.parts_dir))
        self.rows = count_part_rows(self.parts_dir) or 0
        self.columns: list[str] | None = None

    def write(self, df: pd.DataFrame) -> Path | None:
        if df is None or len(df) == 0:
            return None
        self.part_index += 1
        if self.columns is None:
            self.columns = [str(c) for c in df.columns]
        path = self.parts_dir / f"part-{self.part_index:06d}{self.suffix}"
        df.to_csv(path, index=False)
        self.rows += len(df)
        return path


def materialize_csv_if_small(
    parts_dir: Path,
    output_csv: Path,
    row_limit: int,
    columns: list[str] | None = None,
) -> int:
    """Write a compatibility CSV when the sharded table is small enough."""
    total = count_part_rows(parts_dir) or 0
    if total > row_limit:
        if output_csv.exists():
            output_csv.unlink()
        return total

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    first = True
    part_files = list_part_files(parts_dir)
    with output_csv.open("w", encoding="utf-8", newline="") as out_handle:
        for part in part_files:
            with _open_text(part) as in_handle:
                header = in_handle.readline()
                if first:
                    out_handle.write(header)
                    first = False
                shutil.copyfileobj(in_handle, out_handle)
    if first:
        pd.DataFrame(columns=columns or []).to_csv(output_csv, index=False)
    return total


def write_stage_counts(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def write_manifest(path: Path, rows: list[dict[str, Any]]) -> None:
    write_stage_counts(path, rows)


def copy_parts(src_parts: Path, dst_parts: Path) -> None:
    if dst_parts.exists():
        shutil.rmtree(dst_parts)
    dst_parts.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(src_parts, dst_parts)


def resolve_large_output(base_path: Path, *path_parts: str) -> Path | None:
    csv_path = base_path.joinpath(*path_parts)
    return output_path_or_parts(csv_path)
