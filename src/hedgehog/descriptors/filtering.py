"""Threshold-based filtering and pass/fail logic for descriptors."""

import ast
import json
import math
from pathlib import Path

import pandas as pd

from hedgehog.configs.logger import logger
from hedgehog.descriptors.constants import (
    _DESCRIPTOR_KEY_MAP,
    STRUCTURAL_ELEMENT_LIMIT_MAP,
)
from hedgehog.descriptors.io import (
    order_identity_columns,
    process_path,
    save_failed_molecules,
)


def _parse_literal_list(value, label):
    """Parse list-like values from strings safely using ast.literal_eval.

    This function deserializes lists that were saved as strings in CSV files
    (e.g., ['C', 'N', 'O'] for atom symbols or [5, 6] for ring sizes).

    Security: ast.literal_eval is safe for untrusted input. Unlike Python's
    built-in eval function, literal_eval only parses literal data structures
    (strings, numbers, tuples, lists, dicts, sets, booleans, None) and
    rejects any other Python code, preventing code injection attacks.

    Args:
        value: Value to parse (string representation of a list or list itself)
        label: Label for error messages (e.g., "ring sizes", "type limits")

    Returns:
        Parsed list, or empty list if parsing fails
    """
    if isinstance(value, list):
        return value
    if isinstance(value, str):
        try:
            parsed = ast.literal_eval(value)
            if isinstance(parsed, list):
                return parsed
        except Exception as e:
            logger.error("Error parsing %s: %s, %s", label, value, e)
            return []
    logger.error("Error parsing %s: %s, unsupported type", label, value)
    return []


def _parse_ring_size_column(series):
    """Parse ring size lists from series for plotting."""
    parsed = []
    for val in series.dropna():
        sizes = _parse_literal_list(val, "ring sizes")
        try:
            parsed.extend([float(size) for size in sizes])
        except Exception as e:
            logger.error("Error parsing ring sizes: %s, %s", val, e)
    return parsed


def drop_false_rows(df, borders):
    """Filter rows that passed all descriptor filters.

    Args:
        df: DataFrame with '_pass' columns
        borders: Configuration dict with filter settings

    Returns:
        pd.DataFrame: Filtered dataframe with only passed molecules
    """
    passed_cols = [col for col in df.columns if col.endswith("_pass") or col == "pass"]
    mask = df[passed_cols].all(axis=1)
    return df[mask].copy()


def _get_border_values(col, borders):
    """Get min and max border values for a column from borders config.

    Uses explicit key mapping to prevent case-insensitive confusion between
    similar descriptor names with different casing.

    Args:
        col: Column name
        borders: Dictionary with border configurations

    Returns:
        tuple: (min_border, max_border)
    """
    # Try exact match first (preserves case distinction)
    min_key = f"{col}_min"
    max_key = f"{col}_max"
    if min_key in borders or max_key in borders:
        return borders.get(min_key), borders.get(max_key)

    # Try canonical mapping for case-insensitive lookup
    col_lower = col.lower()
    canonical = _DESCRIPTOR_KEY_MAP.get(col_lower, col)
    min_key = f"{canonical}_min"
    max_key = f"{canonical}_max"
    if min_key in borders or max_key in borders:
        return borders.get(min_key), borders.get(max_key)

    # Fallback: search all keys (legacy behavior)
    min_border = None
    max_border = None
    for key, value in borders.items():
        if key.endswith("_min"):
            base = key[: -len("_min")]
            if base.lower() == col_lower:
                min_border = value
        elif key.endswith("_max"):
            base = key[: -len("_max")]
            if base.lower() == col_lower:
                max_border = value
    return min_border, max_border


def _apply_column_filter(df, col, borders):
    """Apply filter to a single column based on borders configuration.

    Args:
        df: DataFrame with data
        col: Column name to filter
        borders: Dictionary with border configurations

    Returns:
        pd.Series: Boolean series indicating pass/fail for each row
    """
    min_border, max_border = _get_border_values(col, borders)

    # Normalize common YAML/string sentinels (e.g., "inf") for numeric comparisons.
    def _normalize_border_value(val):
        if val is None:
            return None
        if isinstance(val, str) and val.lower() == "inf":
            return math.inf
        return val

    min_border = _normalize_border_value(min_border)
    max_border = _normalize_border_value(max_border)

    if col == "ring_size":
        # Missing bounds mean "no bound" on that side.
        if min_border is None:
            min_border = -math.inf
        if max_border is None:
            max_border = math.inf
        return df[col].apply(
            lambda x: all(
                min_border <= float(ring_size) <= max_border
                for ring_size in _parse_literal_list(x, "ring sizes")
            )
        )

    # Generic numeric range filter.
    #
    # Allow one-sided bounds:
    # - only min -> value >= min
    # - only max -> value <= max
    if min_border is None and max_border is None:
        return pd.Series(True, index=df.index)
    if min_border is None:
        return df[col] <= max_border
    if max_border is None or max_border == math.inf:
        return df[col] >= min_border
    return (df[col] >= min_border) & (df[col] <= max_border)


def _extract_borders_and_constraints(borders, structural_constraints):
    """Normalize descriptors config into plain borders + structural constraints."""
    if not isinstance(borders, dict):
        if isinstance(structural_constraints, dict):
            return {}, structural_constraints
        return {}, {}

    normalized_borders = borders
    constraints = structural_constraints

    if "borders" in borders and isinstance(borders.get("borders"), dict):
        normalized_borders = borders.get("borders", {}) or {}
        if constraints is None:
            constraints = borders.get("structural_constraints")
    elif constraints is None:
        constraints = borders.get("structural_constraints")

    normalized_borders = {
        key: value
        for key, value in normalized_borders.items()
        if key != "structural_constraints"
    }
    if not isinstance(constraints, dict):
        constraints = {}
    return normalized_borders, constraints


def _build_structural_constraint_borders(structural_constraints):
    """Convert structural constraints block to regular *_max border entries."""
    if not structural_constraints.get("enabled", False):
        return {}

    translated = {}

    type_limits = structural_constraints.get("type_limits") or {}
    if isinstance(type_limits, dict):
        for alias, max_value in type_limits.items():
            if max_value is None:
                continue
            alias_key = str(alias)
            canonical_alias = _DESCRIPTOR_KEY_MAP.get(alias_key.lower(), alias_key)
            translated[f"{canonical_alias}_max"] = max_value

    element_limits = structural_constraints.get("element_limits") or {}
    if isinstance(element_limits, dict):
        for key, max_value in element_limits.items():
            if max_value is None:
                continue
            key_str = str(key)
            element_col = STRUCTURAL_ELEMENT_LIMIT_MAP.get(
                key_str, STRUCTURAL_ELEMENT_LIMIT_MAP.get(key_str.lower())
            )
            if element_col is not None:
                translated[f"{element_col}_max"] = max_value

    direct_limits = {
        "max_n_or_o_atoms": "n_NO_atoms",
        "max_small_rings_3_4": "n_small_rings_3_4",
        "max_acyclic_chain_length": "max_acyclic_chain_length",
    }
    for config_key, column in direct_limits.items():
        max_value = structural_constraints.get(config_key)
        if max_value is not None:
            translated[f"{column}_max"] = max_value

    return translated


def filter_molecules(df, borders, folder_to_save, structural_constraints=None):
    """Filter molecules based on descriptor thresholds.

    Args:
        df: DataFrame with computed descriptors
        borders: Dictionary with min/max thresholds for each descriptor
        folder_to_save: Output folder path (should already include 'Descriptors' subfolder)
    """
    folder_to_save = Path(process_path(folder_to_save))
    id_cols = ["smiles", "model_name", "mol_idx"]
    normalized_borders, constraints = _extract_borders_and_constraints(
        borders, structural_constraints
    )
    effective_borders = dict(normalized_borders)
    effective_borders.update(_build_structural_constraint_borders(constraints))

    logger.info("[#B29EEE]Applied Descriptor Filters:[/#B29EEE]")
    for line in json.dumps(effective_borders, indent=2, ensure_ascii=False).split("\n"):
        logger.info("  %s", line)

    # Build filtered data with pass flags
    filtered_data = {}
    for col in df.columns.tolist():
        if col in id_cols:
            filtered_data[col] = df[col]
            continue

        col_in_borders = any(
            v is not None for v in _get_border_values(col, effective_borders)
        )
        if col_in_borders:
            filtered_data[col] = df[col]
            filtered_data[f"{col}_pass"] = _apply_column_filter(
                df, col, effective_borders
            )

    filtered_data_df = pd.DataFrame(filtered_data)
    pass_cols = [c for c in filtered_data_df.columns if c.endswith("_pass")]
    if pass_cols:
        filtered_data_df["stage_pass"] = filtered_data_df[pass_cols].all(axis=1)
    else:
        filtered_data_df["stage_pass"] = True
    filtered_data_df = order_identity_columns(filtered_data_df)
    filtered_data_df.to_csv(
        folder_to_save / "pass_flags.csv", index_label="SMILES", index=False
    )

    pass_filters = drop_false_rows(filtered_data_df, effective_borders)

    # Save passed molecules
    if len(pass_filters) > 0:
        pass_filters = order_identity_columns(pass_filters)
        descriptor_cols = [
            col
            for col in pass_filters.columns
            if col not in id_cols and not col.endswith("_pass")
        ]
        ordered_cols = [
            col
            for col in id_cols + sorted(descriptor_cols)
            if col in pass_filters.columns
        ]
        pass_filters[ordered_cols].to_csv(
            folder_to_save / "descriptors_passed.csv", index=False
        )
        pass_filters[id_cols].to_csv(
            folder_to_save / "filtered_molecules.csv", index=False
        )
    else:
        logger.warning("No molecules pass Descriptors Filters")
        descriptor_cols = [
            col
            for col in filtered_data_df.columns
            if col not in id_cols and not col.endswith("_pass")
        ]
        pd.DataFrame(
            columns=[
                col
                for col in id_cols + sorted(descriptor_cols)
                if col in filtered_data_df.columns or col in id_cols
            ]
        ).to_csv(folder_to_save / "descriptors_passed.csv", index=False)
        pd.DataFrame(columns=id_cols).to_csv(
            folder_to_save / "filtered_molecules.csv", index=False
        )

    # Save failed molecules
    all_computed_path = folder_to_save / "descriptors_all.csv"
    if not all_computed_path.exists():
        # descriptors_all.csv is written into the sibling "metrics/" folder
        # by compute_metrics(). Fall back to that location for provenance.
        sibling_metrics = folder_to_save.parent / "metrics" / "descriptors_all.csv"
        if sibling_metrics.exists():
            all_computed_path = sibling_metrics
    flags_path = folder_to_save / "pass_flags.csv"
    all_computed = None

    if all_computed_path.exists():
        all_computed = pd.read_csv(all_computed_path)

        if len(pass_filters) > 0:
            merge_cols = ["smiles", "model_name"]
            merged = all_computed.merge(
                pass_filters[merge_cols], on=merge_cols, how="left", indicator=True
            )
            fail_filters = (
                merged[merged["_merge"] == "left_only"].drop(columns=["_merge"]).copy()
            )
        else:
            fail_filters = all_computed.copy()

        if len(fail_filters) > 0:
            save_failed_molecules(fail_filters, folder_to_save, flags_path)

    # Re-order existing CSV files
    if (folder_to_save / "filtered_molecules.csv").exists():
        if all_computed is not None:
            order_identity_columns(all_computed).to_csv(all_computed_path, index=False)

        if flags_path.exists():
            flags = pd.read_csv(flags_path)
            order_identity_columns(flags).to_csv(flags_path, index=False)
