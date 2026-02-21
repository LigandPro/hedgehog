import os
from pathlib import Path

import numpy as np
import pandas as pd

from hedgehog.configs.logger import load_config, logger
from hedgehog.stages.structFilters import filters as _filters
from hedgehog.stages.structFilters import plotting as _plotting
from hedgehog.stages.structFilters import stats as _stats
from hedgehog.utils.datamol_import import import_datamol_quietly
from hedgehog.utils.parallel import parallel_map, resolve_n_jobs
from hedgehog.utils.paths import process_path  # noqa: F401

# Mirror filter-module globals for backward compatibility and monkeypatch-based tests.
mc = _filters.mc
LILLY_AVAILABLE = _filters.LILLY_AVAILABLE
LillyDemeritsFilters = _filters.LillyDemeritsFilters
_LILLY_DEFAULT_COLUMNS = _filters._LILLY_DEFAULT_COLUMNS
_ORIG_FILTERS_ADD_MODEL_NAME_COL = _filters.add_model_name_col

dm = import_datamol_quietly()


def _sync_filters_module_state() -> None:
    _filters.load_config = load_config
    _filters.logger = logger
    _filters.parallel_map = parallel_map
    _filters.resolve_n_jobs = resolve_n_jobs
    _filters.dm = dm
    _filters.mc = mc
    _filters.LILLY_AVAILABLE = LILLY_AVAILABLE
    _filters.LillyDemeritsFilters = LillyDemeritsFilters
    _filters.filter_alerts = filter_alerts


def _sync_stats_module_state() -> None:
    _stats.mc = mc
    _stats.logger = logger
    _stats.process_path = process_path


def _sync_plotting_module_state() -> None:
    _plotting.load_config = load_config
    _plotting.logger = logger
    _plotting.process_path = process_path


def _silence_worker_stdio() -> None:
    """Redirect worker stdout/stderr to /dev/null."""
    try:
        devnull_fd = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull_fd, 1)
        os.dup2(devnull_fd, 2)
        os.close(devnull_fd)
    except OSError:
        pass


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


def camelcase(any_str):
    """Convert underscore-separated string to CamelCase."""
    return "".join(word.capitalize() for word in any_str.split("_"))


def build_identity_map_from_descriptors(config):
    """Build a map of (smiles, model_name) -> mol_idx from descriptors output."""
    base_folder = Path(process_path(config["folder_to_save"]))
    id_path = base_folder / "Descriptors" / "passDescriptorsSMILES.csv"

    try:
        if id_path.exists():
            id_df = pd.read_csv(id_path)
            identity_map = {
                (row["smiles"], row["model_name"]): row["mol_idx"]
                for _, row in id_df.iterrows()
            }
            return identity_map, id_df
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError, KeyError):
        pass

    return {}, None


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


def _apply_subsample_to_dataframe(df, subsample):
    """Apply structFilters subsampling policy to input DataFrame."""
    if subsample is None or subsample <= 0:
        return df

    if "model_name" in df.columns and df["model_name"].nunique(dropna=True) > 1:
        return df.groupby("model_name", group_keys=False).head(subsample)

    return df.head(subsample) if len(df) > subsample else df


def _parse_smiles_item(args):
    """Parse one SMILES string into RDKit Mol."""
    row_idx, smiles_raw, model_name, mol_idx = args
    try:
        mol = dm.to_mol(smiles_raw, sanitize=True)
    except (ValueError, TypeError):
        mol = None
    return row_idx, smiles_raw, model_name, mol_idx, mol


def prepare_structfilters_input(df, subsample, parse_n_jobs, progress_cb=None):
    """Prepare structFilters payload once: subsample + SMILES parsing.

    Returns:
        dict with keys:
          - ``mols``: list of valid RDKit mols
          - ``smiles_model_mols``: tuples (smiles, model_name_or_none, mol, mol_idx)
          - ``base_df``: DataFrame with only valid rows
    """
    if df is None or len(df) == 0:
        return {
            "mols": [],
            "smiles_model_mols": [],
            "base_df": pd.DataFrame(columns=["smiles", "model_name", "mol_idx"]),
        }

    if "smiles" not in df.columns:
        raise KeyError("Input DataFrame must contain a 'smiles' column")

    data = _apply_subsample_to_dataframe(df, subsample).copy()
    if "model_name" not in data.columns:
        data["model_name"] = "single"
    if "mol_idx" not in data.columns:
        data["mol_idx"] = range(len(data))

    if len(data) == 0:
        return {
            "mols": [],
            "smiles_model_mols": [],
            "base_df": data.iloc[0:0].copy(),
        }

    is_multi = data["model_name"].nunique(dropna=True) > 1
    model_vals = data["model_name"].tolist() if is_multi else [None] * len(data)
    items = list(
        enumerate(
            zip(
                data["smiles"].tolist(),
                model_vals,
                data["mol_idx"].tolist(),
                strict=False,
            )
        )
    )
    parse_payload = [(row_idx, item[0], item[1], item[2]) for row_idx, item in items]

    n_jobs = int(parse_n_jobs)
    if n_jobs <= 0:
        n_jobs = os.cpu_count() or 1
    parsed = parallel_map(
        _parse_smiles_item,
        parse_payload,
        n_jobs,
        progress=progress_cb,
        initializer=_silence_worker_stdio
        if n_jobs > 1 and len(parse_payload) > 1
        else None,
    )

    valid_row_indices = []
    smiles_model_mols = []
    mols = []
    for row_idx, smiles_raw, model_name, mol_idx, mol in parsed:
        if mol is None:
            continue
        valid_row_indices.append(row_idx)
        smiles_model_mols.append((smiles_raw, model_name, mol, mol_idx))
        mols.append(mol)

    if valid_row_indices:
        base_df = data.iloc[valid_row_indices].copy().reset_index(drop=True)
    else:
        base_df = data.iloc[0:0].copy()

    return {
        "mols": mols,
        "smiles_model_mols": smiles_model_mols,
        "base_df": base_df,
    }


def process_prepared_payload(config, prepared_payload, apply_filter, progress_cb=None):
    """Apply one structural filter on pre-parsed payload."""
    mols = prepared_payload["mols"]
    smiles = prepared_payload["smiles_model_mols"]

    if len(mols) == 0:
        return None

    if progress_cb is not None:
        try:
            return apply_filter(config, mols, smiles, progress_cb=progress_cb)
        except TypeError:
            pass
    return apply_filter(config, mols, smiles)


def process_one_dataframe(config, df, apply_filter, subsample, progress_cb=None):
    """Process a DataFrame directly without re-reading from file.

    Args:
        config: Configuration dictionary
        df: DataFrame with 'smiles', 'model_name', and 'mol_idx' columns
        apply_filter: Filter function to apply
        subsample: Maximum number of molecules to process

    Returns:
        DataFrame with filter results or None if no molecules to process
    """
    if df is None or len(df) == 0:
        return None

    if subsample is None or subsample <= 0:
        data = df
    elif "model_name" in df.columns and df["model_name"].nunique(dropna=True) > 1:
        # In multi-model mode, apply the subsample limit per model, not globally.
        # The global sampling is handled earlier by prepare_input_data().
        data = df.groupby("model_name", group_keys=False).head(subsample)
    else:
        data = df.head(subsample) if len(df) > subsample else df
    smiles_col = "smiles"

    is_multi = data["model_name"].nunique(dropna=True) > 1
    if is_multi:
        smiles_str = data[smiles_col].tolist()
        model_names = data["model_name"].tolist()
        mols = [dm.to_mol(x, sanitize=True) for x in smiles_str]
        mol_indices = data["mol_idx"].tolist()
        smiles = list(zip(smiles_str, model_names, mols, mol_indices, strict=False))
    else:
        smiles_list = data[smiles_col].tolist()
        mols = [dm.to_mol(x, sanitize=True) for x in smiles_list]
        mol_indices = data["mol_idx"].tolist()
        smiles = list(
            zip(smiles_list, [None] * len(smiles_list), mols, mol_indices, strict=False)
        )

    # Filter out invalid molecules
    cleaned = []
    for item in smiles:
        if len(item) >= 3:
            smi, model, mol = item[0], item[1], item[2]
            if mol is not None:
                if len(item) >= 4:
                    mol_idx = item[3]
                    cleaned.append((smi, model, mol, mol_idx))
                else:
                    cleaned.append((smi, model, mol))
    smiles = cleaned
    mols = [it[2] for it in smiles]

    if len(mols) == 0:
        return None

    if progress_cb is not None:
        try:
            return apply_filter(config, mols, smiles, progress_cb=progress_cb)
        except TypeError:
            pass
    return apply_filter(config, mols, smiles)


def process_one_file(config, input_path, apply_filter, subsample, progress_cb=None):
    """Process molecules from file through a filter function.

    For CSV files, consider using process_one_dataframe() directly if you
    already have the DataFrame loaded to avoid reading the file twice.

    Args:
        config: Configuration dictionary
        input_path: Path to input file (csv, smi, sdf, or txt)
        apply_filter: Filter function to apply
        subsample: Maximum number of molecules to process

    Returns:
        DataFrame with filter results or None if no molecules to process
    """
    input_type = input_path[input_path.rfind(".") + 1 :]
    assert input_type in {"csv", "smi", "sdf", "txt"}

    if input_type == "csv":
        data = pd.read_csv(input_path)
        return process_one_dataframe(config, data, apply_filter, subsample)

    elif input_type == "smi" or input_type == "txt":
        with open(input_path) as file:
            lines = [line.rstrip("\n") for line in file]

        if subsample <= len(lines):
            lines = np.random.permutation(lines)[:subsample].tolist()

        smiles = []
        model_names = []
        for line in lines:
            parts = line.split(",")
            if len(parts) == 2:
                smi, model = parts
                smiles.append(smi)
                model_names.append(model)
            else:
                smiles.append(parts[0])
        mols = [dm.to_mol(x) for x in smiles]
        if len(model_names) == len(smiles):
            smiles = list(zip(smiles, model_names, mols, strict=False))

    elif input_type == "sdf":
        mols, smiles = sdf_to_mols(input_path, subsample)

    if isinstance(smiles[0], tuple):
        cleaned = []
        for item in smiles:
            if len(item) >= 3:
                smi, model, mol = item[0], item[1], item[2]
                if mol is not None:
                    if len(item) >= 4:
                        mol_idx = item[3]
                        cleaned.append((smi, model, mol, mol_idx))
                    else:
                        cleaned.append((smi, model, mol))
        smiles = cleaned
        mols = [it[2] for it in smiles]
    else:
        mols, smiles = dropna(mols, smiles)

    assert len(mols) == len(smiles), f"{len(mols)}, {len(smiles)}"
    if not isinstance(smiles[0], tuple):
        assert len(mols) <= subsample

    for mol, smi in zip(mols, smiles, strict=False):
        smi_val = smi[0] if isinstance(smi, tuple) else smi
        assert mol is not None, f"{smi_val}"

    final_result = None
    if len(mols) > 0:
        if isinstance(smiles[0], tuple):
            if progress_cb is not None:
                try:
                    final_result = apply_filter(
                        config, mols, smiles, progress_cb=progress_cb
                    )
                except TypeError:
                    final_result = apply_filter(config, mols, smiles)
            else:
                final_result = apply_filter(config, mols, smiles)
        else:
            if progress_cb is not None:
                try:
                    final_result = apply_filter(config, mols, progress_cb=progress_cb)
                except TypeError:
                    final_result = apply_filter(config, mols)
            else:
                final_result = apply_filter(config, mols)

    return final_result


def inject_identity_columns_to_all_csvs(config, stage_dir):
    """Ensure identity columns are ordered consistently in all CSVs."""
    base_folder = Path(process_path(config["folder_to_save"]))
    target_folder = base_folder / stage_dir

    csv_paths = [str(p) for p in target_folder.glob("*.csv")]
    for path in csv_paths:
        try:
            df = pd.read_csv(path)
            if "smiles" not in df.columns:
                continue

            identity_order = ["smiles", "model_name", "mol_idx"]
            ordered = [c for c in identity_order if c in df.columns] + [
                c for c in df.columns if c not in identity_order
            ]
            df = df[ordered]
            df.to_csv(path, index=False)
        except (OSError, pd.errors.ParserError, KeyError):
            continue


def format_number(x, pos=None):
    return _plotting.format_number(x, pos)


def get_model_colors(model_names, cmap=None):
    return _plotting.get_model_colors(model_names, cmap)


def clean_name(name):
    return _plotting.clean_name(name)


def check_paths(config, paths):
    _sync_plotting_module_state()
    return _plotting.check_paths(config, paths)


def plot_calculated_stats(config, stage_dir):
    _sync_plotting_module_state()
    return _plotting.plot_calculated_stats(config, stage_dir)


def plot_restriction_ratios(config, stage_dir):
    _sync_plotting_module_state()
    return _plotting.plot_restriction_ratios(config, stage_dir)


def plot_filter_failures_analysis(config, stage_dir):
    _sync_plotting_module_state()
    return _plotting.plot_filter_failures_analysis(config, stage_dir)


def common_postprocessing_statistics(filter_results, res_df, stat, extend):
    return _stats.common_postprocessing_statistics(filter_results, res_df, stat, extend)


def get_basic_stats(
    config, filter_results, model_name, filter_name, stat=None, extend=None
):
    _sync_stats_module_state()
    return _stats.get_basic_stats(
        config, filter_results, model_name, filter_name, stat=stat, extend=extend
    )


def combine_filter_results_in_memory(output_dir, input_df, pass_mask_by_filter):
    return _stats.combine_filter_results_in_memory(
        output_dir, input_df, pass_mask_by_filter
    )


def filter_data(config, stage_dir):
    _sync_stats_module_state()
    return _stats.filter_data(config, stage_dir)


def filter_alerts(config):
    return _filters.filter_alerts(config)


def add_model_name_col(final_result, smiles_with_model):
    return _ORIG_FILTERS_ADD_MODEL_NAME_COL(final_result, smiles_with_model)


def apply_structural_alerts(config, mols, smiles_modelName_mols=None, progress_cb=None):
    _sync_filters_module_state()
    return _filters.apply_structural_alerts(
        config,
        mols,
        smiles_modelName_mols=smiles_modelName_mols,
        progress_cb=progress_cb,
    )


def apply_molgraph_stats(config, mols, smiles_modelName_mols=None):
    _sync_filters_module_state()
    return _filters.apply_molgraph_stats(
        config, mols, smiles_modelName_mols=smiles_modelName_mols
    )


def apply_molcomplexity_filters(config, mols, smiles_modelName_mols=None):
    _sync_filters_module_state()
    return _filters.apply_molcomplexity_filters(
        config, mols, smiles_modelName_mols=smiles_modelName_mols
    )


def apply_bredt_filter(config, mols, smiles_modelName_mols=None):
    _sync_filters_module_state()
    return _filters.apply_bredt_filter(
        config, mols, smiles_modelName_mols=smiles_modelName_mols
    )


def apply_protecting_groups(config, mols, smiles_modelName_mols=None):
    _sync_filters_module_state()
    return _filters.apply_protecting_groups(
        config, mols, smiles_modelName_mols=smiles_modelName_mols
    )


def apply_ring_infraction(config, mols, smiles_modelName_mols=None):
    _sync_filters_module_state()
    return _filters.apply_ring_infraction(
        config, mols, smiles_modelName_mols=smiles_modelName_mols
    )


def apply_stereo_center(config, mols, smiles_modelName_mols=None):
    _sync_filters_module_state()
    return _filters.apply_stereo_center(
        config, mols, smiles_modelName_mols=smiles_modelName_mols
    )


def apply_halogenicity(config, mols, smiles_modelName_mols=None):
    _sync_filters_module_state()
    return _filters.apply_halogenicity(
        config, mols, smiles_modelName_mols=smiles_modelName_mols
    )


def apply_symmetry(config, mols, smiles_modelName_mols=None):
    _sync_filters_module_state()
    return _filters.apply_symmetry(
        config, mols, smiles_modelName_mols=smiles_modelName_mols
    )


def apply_nibr_filter(config, mols, smiles_modelName_mols=None):
    _sync_filters_module_state()
    return _filters.apply_nibr_filter(
        config, mols, smiles_modelName_mols=smiles_modelName_mols
    )


def apply_lilly_filter(config, mols, smiles_modelName_mols=None):
    _sync_filters_module_state()
    return _filters.apply_lilly_filter(
        config, mols, smiles_modelName_mols=smiles_modelName_mols
    )


def filter_function_applier(filter_name):
    filters = {
        "common_alerts": apply_structural_alerts,
        "molgraph_stats": apply_molgraph_stats,
        "molcomplexity": apply_molcomplexity_filters,
        "NIBR": apply_nibr_filter,
        "bredt": apply_bredt_filter,
        "lilly": apply_lilly_filter,
        "protecting_groups": apply_protecting_groups,
        "ring_infraction": apply_ring_infraction,
        "stereo_center": apply_stereo_center,
        "halogenicity": apply_halogenicity,
        "symmetry": apply_symmetry,
    }
    if filter_name not in filters:
        raise ValueError(f"Filter {filter_name} not found")
    return filters[filter_name]
