"""Shared DataFrame utilities for column ordering and CSV I/O."""

from pathlib import Path

import pandas as pd

IDENTITY_COLUMNS = ["smiles", "model_name", "mol_idx"]


def order_identity_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Order dataframe columns with identity columns first."""
    existing_id_cols = [c for c in IDENTITY_COLUMNS if c in df.columns]
    ordered_cols = existing_id_cols + [
        c for c in df.columns if c not in IDENTITY_COLUMNS
    ]
    return df[ordered_cols]


def save_ordered_csv(df: pd.DataFrame, path: Path) -> None:
    """Save DataFrame with identity columns ordered first."""
    path.parent.mkdir(parents=True, exist_ok=True)
    order_identity_columns(df).to_csv(path, index=False)
