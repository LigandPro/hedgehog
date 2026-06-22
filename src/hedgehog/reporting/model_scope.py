"""Generative model scoping helpers for report scoring and aggregation."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

MODEL_INDEX_MAP_REL_PATH = Path("configs") / "model_index_map.json"


def load_model_index_map(base_path: Path) -> dict[str, int]:
    """Load generative model names and numeric indices for a run."""
    map_path = base_path / MODEL_INDEX_MAP_REL_PATH
    if not map_path.exists():
        return {}
    try:
        data = json.loads(map_path.read_text(encoding="utf-8"))
        return {str(name): int(index) for name, index in data.items()}
    except (OSError, json.JSONDecodeError, TypeError, ValueError):
        return {}


def available_models_from_map(model_index_map: dict[str, int]) -> list[str]:
    if not model_index_map:
        return []
    return sorted(
        model_index_map.keys(), key=lambda name: (model_index_map[name], name)
    )


def load_input_molecules_df(base_path: Path) -> pd.DataFrame | None:
    input_path = base_path / "input" / "sampled_molecules.csv"
    if not input_path.exists():
        return None
    try:
        return pd.read_csv(input_path)
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
        return None


def get_available_models(base_path: Path) -> list[str]:
    """Return generative model names from model_index_map.json or input CSV."""
    models = available_models_from_map(load_model_index_map(base_path))
    if models:
        return models

    input_df = load_input_molecules_df(base_path)
    if input_df is not None and "model_name" in input_df.columns:
        return sorted(input_df["model_name"].dropna().astype(str).unique().tolist())
    return []


def _mol_idx_match_columns(df: pd.DataFrame) -> list[str]:
    return [col for col in ("mol_idx", "source_mol_idx") if col in df.columns]


def _filter_by_mol_idx_prefix(
    df: pd.DataFrame, model_number: int, columns: list[str]
) -> pd.DataFrame:
    prefix = f"LP-{model_number:04d}-"
    masks = [df[col].astype(str).str.startswith(prefix) for col in columns]
    if not masks:
        return df.iloc[0:0]
    combined = masks[0]
    for mask in masks[1:]:
        combined = combined | mask
    return df[combined]


def filter_df_by_model(
    df: pd.DataFrame | None, model: str, base_path: Path
) -> pd.DataFrame | None:
    """Filter rows for a generative model using model_name or mol_idx."""
    if df is None:
        return None
    if df.empty:
        return df

    if "model_name" in df.columns:
        by_name = df[df["model_name"].astype(str) == model]
        if not by_name.empty:
            return by_name

    mol_idx_columns = _mol_idx_match_columns(df)
    if not mol_idx_columns:
        return df.iloc[0:0]

    model_index_map = load_model_index_map(base_path)
    model_number = model_index_map.get(model)
    if model_number is not None:
        by_mol_idx = _filter_by_mol_idx_prefix(df, model_number, mol_idx_columns)
        if not by_mol_idx.empty:
            return by_mol_idx

    input_df = load_input_molecules_df(base_path)
    if (
        input_df is not None
        and "mol_idx" in input_df.columns
        and "model_name" in input_df.columns
    ):
        model_mol_idxs = set(
            input_df[input_df["model_name"].astype(str) == model]["mol_idx"]
            .astype(str)
            .tolist()
        )
        masks = [df[col].astype(str).isin(model_mol_idxs) for col in mol_idx_columns]
        combined = masks[0]
        for mask in masks[1:]:
            combined = combined | mask
        return df[combined]

    return df.iloc[0:0]
