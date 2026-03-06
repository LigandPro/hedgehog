"""Internal helpers shared across struct_filters submodules."""

import importlib
import os

import numpy as np
import pandas as pd

from hedgehog.configs.logger import logger
from hedgehog.utils.datamol_import import import_datamol_quietly

_VALID_SCHEDULERS = {"threads", "processes"}

# Default columns for Lilly filter results
_LILLY_DEFAULT_COLUMNS = ["smiles", "status", "pass_filter", "demerit_score", "reasons"]


# ---------------------------------------------------------------------------
# Medchem / datamol lazy imports
# ---------------------------------------------------------------------------


def _import_medchem_quietly():
    """Import medchem while muting native stdout/stderr noise."""
    stdout_fd = os.dup(1)
    stderr_fd = os.dup(2)
    devnull_fd = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull_fd, 1)
        os.dup2(devnull_fd, 2)
        return importlib.import_module("medchem")
    finally:
        os.dup2(stdout_fd, 1)
        os.dup2(stderr_fd, 2)
        os.close(devnull_fd)
        os.close(stdout_fd)
        os.close(stderr_fd)


mc = _import_medchem_quietly()

try:
    from medchem.structural.lilly_demerits import LillyDemeritsFilters

    LILLY_AVAILABLE = True
except ImportError:
    LILLY_AVAILABLE = False
    LillyDemeritsFilters = None

dm = import_datamol_quietly()


# ---------------------------------------------------------------------------
# Lilly PATH setup
# ---------------------------------------------------------------------------

_LILLY_BIN_PATH = (
    __import__("pathlib").Path(__file__).parent
    / ".."
    / ".."
    / ".."
    / "modules"
    / "lilly_medchem_rules"
    / "bin"
).resolve()
if _LILLY_BIN_PATH.exists():
    current_path = os.environ.get("PATH", "")
    lilly_str = str(_LILLY_BIN_PATH)
    if lilly_str not in current_path:
        os.environ["PATH"] = f"{lilly_str}:{current_path}"


# ---------------------------------------------------------------------------
# Worker helpers
# ---------------------------------------------------------------------------


def _silence_worker_stdio() -> None:
    """Redirect worker stdout/stderr to /dev/null."""
    try:
        devnull_fd = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull_fd, 1)
        os.dup2(devnull_fd, 2)
        os.close(devnull_fd)
    except OSError:
        pass


def _resolve_scheduler(config_struct_filters, scheduler_key, default="processes"):
    """Resolve scheduler name for medchem parallel APIs."""
    scheduler_raw = config_struct_filters.get(scheduler_key)
    if scheduler_raw is None:
        scheduler_raw = config_struct_filters.get("parallel_scheduler", default)
    scheduler = str(scheduler_raw).strip().lower()
    if scheduler not in _VALID_SCHEDULERS:
        logger.warning(
            "Invalid scheduler '%s' for '%s'. Using '%s'.",
            scheduler_raw,
            scheduler_key,
            default,
        )
        return default
    return scheduler


def _split_indexed_mols(indexed_mols, n_chunks):
    """Split indexed molecules into near-equal chunks."""
    if not indexed_mols:
        return []
    n_chunks = max(1, min(n_chunks, len(indexed_mols)))
    chunk_size = (len(indexed_mols) + n_chunks - 1) // n_chunks
    return [
        indexed_mols[i : i + chunk_size]
        for i in range(0, len(indexed_mols), chunk_size)
    ]


# ---------------------------------------------------------------------------
# String / display helpers
# ---------------------------------------------------------------------------


def camelcase(any_str):
    """Convert underscore-separated string to CamelCase."""
    return "".join(word.capitalize() for word in any_str.split("_"))


def format_number(x, pos=None):
    """Format number for display. pos parameter is for matplotlib FuncFormatter compatibility."""
    if x >= 1e6:
        return f"{x / 1e6:.1f}M"
    elif x >= 1e3:
        return f"{x / 1e3:.1f}K"
    return f"{x:.0f}"


def clean_name(name):
    """Clean metric names for display."""
    for pattern in ["metrics", "_", ".csv"]:
        name = name.replace(pattern, "")
    return name.strip()


# ---------------------------------------------------------------------------
# Molecule helpers
# ---------------------------------------------------------------------------


def sdf_to_mols(sdf_file, subsample):
    """Read molecules from SDF file with subsampling."""
    molecules = dm.read_sdf(sdf_file)
    mols_list = []
    smiles_list = []

    for i, mol in enumerate(molecules):
        if i >= subsample:
            break
        if mol is not None:
            mols_list.append(mol)
            smiles_list.append(dm.to_smiles(mol))

    return mols_list, smiles_list


def dropna(mols, smiles):
    """Remove None values from molecule and SMILES lists."""
    df = pd.DataFrame({"mols": mols, "smiles": smiles}).dropna()
    return df["mols"].tolist(), df["smiles"].tolist()


def get_model_colors(model_names, cmap=None):
    """Generate color map for models."""
    import matplotlib.pyplot as plt

    if cmap is None:
        colors = plt.cm.YlOrRd(np.linspace(1, 0, len(model_names) + 1))
    else:
        colors = plt.colormaps.get_cmap(cmap)(np.linspace(1, 0, len(model_names) + 1))

    return dict(zip(model_names, colors, strict=False))


# ---------------------------------------------------------------------------
# DataFrame alignment / padding
# ---------------------------------------------------------------------------


def _create_failed_row(smiles=None, reason="unsupported_or_missing"):
    """Create a failed molecule row for Lilly filter results."""
    return {
        "smiles": smiles,
        "status": "exclude",
        "pass_filter": False,
        "demerit_score": None,
        "reasons": reason,
    }


def _create_failed_dataframe(count, smiles_list=None, reason="unsupported_or_missing"):
    """Create a DataFrame of failed molecule rows."""
    if smiles_list is None:
        smiles_list = [None] * count
    return pd.DataFrame([_create_failed_row(smi, reason) for smi in smiles_list])


def _pad_dataframe_to_length(df, target_length, template_row=None):
    """Pad a DataFrame to reach target length using failed rows."""
    current_length = len(df)
    if current_length >= target_length:
        return df

    missing_count = target_length - current_length
    if template_row is not None:
        missing_rows = []
        for _ in range(missing_count):
            row = template_row.copy()
            row["smiles"] = None
            row["status"] = "exclude"
            row["pass_filter"] = False
            row["demerit_score"] = None
            if "reasons" in row:
                row["reasons"] = "unsupported_or_missing"
            missing_rows.append(row)
        missing_df = pd.DataFrame(missing_rows)
    else:
        missing_df = _create_failed_dataframe(missing_count)

    return pd.concat([df, missing_df], ignore_index=True)


def _ensure_dataframe_length(df, expected_length, template_row=None):
    """Ensure DataFrame has exactly the expected length by padding or trimming."""
    if len(df) < expected_length:
        return _pad_dataframe_to_length(df, expected_length, template_row)
    if len(df) > expected_length:
        return df.iloc[:expected_length].reset_index(drop=True)
    return df


def _filter_valid_molecule_tuples(smiles):
    """Filter out invalid molecules from a list of (smiles, model, mol, ...) tuples."""
    cleaned = []
    for item in smiles:
        if len(item) < 3:
            continue
        smi, model, mol = item[0], item[1], item[2]
        if mol is None:
            continue
        if len(item) >= 4:
            cleaned.append((smi, model, mol, item[3]))
        else:
            cleaned.append((smi, model, mol))
    return cleaned


def _align_result_length(final_result, expected_len, smiles_with_model):
    """Pad or trim final_result to match expected_len."""
    actual_len = len(final_result)
    if actual_len == expected_len:
        return final_result

    logger.error(
        "Length mismatch in add_model_name_col: final_result has %d rows, "
        "but smiles_with_model has %d items.",
        actual_len,
        expected_len,
    )

    if actual_len < expected_len:
        logger.warning(
            "Padding final_result from %d to %d rows.", actual_len, expected_len
        )
        if len(final_result) > 0:
            template_row = final_result.iloc[-1].to_dict()
            final_result = _pad_dataframe_to_length(
                final_result, expected_len, template_row
            )
            logger.info("Padded to %d rows.", len(final_result))
        else:
            logger.error("final_result is empty. Creating all rows from scratch.")
            smiles_list = [
                item[0] if isinstance(item, tuple) and len(item) > 0 else None
                for item in smiles_with_model
            ]
            final_result = _create_failed_dataframe(expected_len, smiles_list)
    else:
        logger.warning(
            "Trimming final_result from %d to %d rows.", actual_len, expected_len
        )
        final_result = final_result.iloc[:expected_len].reset_index(drop=True)

    if len(final_result) != expected_len:
        msg = (
            f"CRITICAL: After padding/trimming, final_result length ({len(final_result)}) "
            f"still doesn't match expected ({expected_len})."
        )
        raise ValueError(msg)

    return final_result


def _realign_result_by_temp_mol_idx(final_result):
    """Restore original molecule order when worker code returns ``_mol_idx``."""
    if "_mol_idx" not in final_result.columns:
        return final_result
    return final_result.sort_values("_mol_idx").reset_index(drop=True)


def add_model_name_col(final_result, smiles_with_model):
    """Add smiles, model_name, and mol_idx columns from input data."""
    expected_len = len(smiles_with_model)
    final_result = _realign_result_by_temp_mol_idx(final_result)
    final_result = _align_result_length(final_result, expected_len, smiles_with_model)

    smiles_vals = [item[0] for item in smiles_with_model]
    model_vals = [
        item[1] if item[1] is not None else "single" for item in smiles_with_model
    ]
    mol_idx_vals = [item[3] if len(item) >= 4 else None for item in smiles_with_model]

    final_result["smiles"] = smiles_vals
    final_result["model_name"] = model_vals
    final_result["mol_idx"] = mol_idx_vals
    if "_mol_idx" in final_result.columns:
        final_result = final_result.drop(columns=["_mol_idx"])

    return final_result


def filter_alerts(config):
    """Filter structural alerts based on configuration."""
    df = pd.read_csv(config["alerts_data_path"])
    mask = df["rule_set_name"].isin(config["include_rulesets"])

    for ruleset in config["exclude_descriptions"]:
        is_other_ruleset = df["rule_set_name"] != ruleset
        is_not_excluded = ~df["description"].isin(
            config["exclude_descriptions"][ruleset]
        )
        mask &= is_other_ruleset | is_not_excluded

    return df[mask]


def common_postprocessing_statistics(filter_results, res_df, stat, extend):
    """Combine filter results with statistics."""
    if stat is not None:
        res_df = pd.concat([stat, res_df])

    filter_results = filter_results.drop(columns="mol")
    if extend is not None:
        filter_extended = pd.concat([extend, filter_results], ignore_index=True)
    else:
        filter_extended = filter_results.copy()

    return res_df, filter_extended


def _apply_subsample_to_dataframe(df, subsample):
    """Apply structFilters subsampling policy to input DataFrame."""
    if subsample is None or subsample <= 0:
        return df

    if "model_name" in df.columns and df["model_name"].nunique(dropna=True) > 1:
        return df.groupby("model_name", group_keys=False).head(subsample)

    return df.head(subsample) if len(df) > subsample else df
