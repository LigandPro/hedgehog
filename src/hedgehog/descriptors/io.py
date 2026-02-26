"""File I/O operations for descriptors stage."""

from pathlib import Path

import pandas as pd

from hedgehog.struct_filters.utils import process_path  # noqa: F401 – re-export


def order_identity_columns(df):
    """Reorder dataframe columns with identity columns first."""
    id_cols = ["smiles", "model_name", "mol_idx"]
    existing_id_cols = [c for c in id_cols if c in df.columns]
    ordered = existing_id_cols + [c for c in df.columns if c not in id_cols]
    return df[ordered]


def order_descriptor_columns(df):
    """Order DataFrame columns with identity cols, then descriptors with their pass flags.

    Args:
        df: DataFrame to reorder

    Returns:
        list: Ordered column names
    """
    id_cols = ["smiles", "model_name", "mol_idx"]
    descriptor_cols = [
        col for col in df.columns if col not in id_cols and not col.endswith("_pass")
    ]
    pass_cols = [col for col in df.columns if col.endswith("_pass")]

    ordered_cols = id_cols.copy()
    for desc_col in sorted(descriptor_cols):
        if desc_col in df.columns:
            ordered_cols.append(desc_col)
            pass_col = f"{desc_col}_pass"
            if pass_col in pass_cols:
                ordered_cols.append(pass_col)

    for pass_col in sorted(pass_cols):
        if pass_col not in ordered_cols:
            ordered_cols.append(pass_col)

    return [col for col in ordered_cols if col in df.columns]


def merge_pass_flags(df, flags_path):
    """Merge pass flags from flags file into dataframe.

    Args:
        df: DataFrame to merge flags into
        flags_path: Path to flags CSV file

    Returns:
        pd.DataFrame: DataFrame with merged pass flags
    """
    if not Path(flags_path).exists():
        return df

    flags_df = pd.read_csv(flags_path)
    merge_cols = ["smiles", "model_name"]
    pass_cols = [
        col for col in flags_df.columns if col.endswith("_pass") or col == "pass"
    ]

    if not pass_cols:
        return df

    df = df.merge(
        flags_df[merge_cols + pass_cols],
        on=merge_cols,
        how="left",
        suffixes=("", "_flags"),
    )
    for col in pass_cols:
        if f"{col}_flags" in df.columns:
            df[col] = df[f"{col}_flags"].fillna(df.get(col, False))
            df = df.drop(columns=[f"{col}_flags"])

    return df


def save_failed_molecules(fail_filters, folder_to_save, flags_path):
    """Save failed molecules to CSV files.

    Args:
        fail_filters: DataFrame with failed molecules
        folder_to_save: Output folder path (Path object)
        flags_path: Path to pass flags CSV
    """
    fail_filters = merge_pass_flags(fail_filters, flags_path)
    fail_filters = order_identity_columns(fail_filters)

    ordered_cols = order_descriptor_columns(fail_filters)
    fail_filters[ordered_cols].to_csv(
        folder_to_save / "descriptors_failed.csv", index=False
    )

    id_cols = ["smiles", "model_name", "mol_idx"]
    fail_filters[id_cols].to_csv(folder_to_save / "failed_molecules.csv", index=False)
