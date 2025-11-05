"""Utility functions for structural filters stage."""

import os
import warnings
from pathlib import Path
from typing import Any

import datamol as dm
import matplotlib.pyplot as plt
import medchem as mc
import numpy as np
import pandas as pd
import seaborn as sns
from pandarallel import pandarallel
from rdkit import Chem

try:
    from medchem.structural.lilly_demerits import LillyDemeritsFilters

    LILLY_AVAILABLE = True
except ImportError:
    LILLY_AVAILABLE = False
    LillyDemeritsFilters = None

from hedge.configs.logger import load_config, logger

warnings.filterwarnings("ignore", category=FutureWarning)

def camelcase(any_str: str) -> str:
    """Convert snake_case to CamelCase.

    Args:
        any_str: String in snake_case format

    Returns
    -------
        String in CamelCase format
    """
    return "".join(word.capitalize() for word in any_str.split("_"))


def build_identity_map_from_descriptors(
    config: dict[str, Any],
) -> tuple[dict[tuple[str, str], Any], pd.DataFrame | None]:
    """Build a map of (smiles, model_name) -> mol_idx from descriptors output.

    Args:
        config: Configuration dictionary

    Returns
    -------
        Tuple of (identity_map, dataframe) or ({}, None) if not found
    """
    base_folder = process_path(config["folder_to_save"])
    id_path = Path(base_folder + "Descriptors/passDescriptorsSMILES.csv")

    try:
        if id_path.exists():
            id_df = pd.read_csv(id_path)
            identity_map = {
                (row["smiles"], row["model_name"]): row["mol_idx"]
                for _, row in id_df.iterrows()
            }
            return identity_map, id_df
    except Exception:
        pass

    return {}, None


def process_path(folder_to_save: str, key_word: str | None = None) -> str:
    """Ensure path ends with '/' and create directory if needed.

    Args:
        folder_to_save: Base folder path
        key_word: Optional subfolder name

    Returns
    -------
        Processed folder path ending with /
    """
    if not folder_to_save.endswith("/"):
        folder_to_save += "/"

    if key_word:
        folder_to_save += f"{key_word}/"

    os.makedirs(folder_to_save, exist_ok=True)
    return folder_to_save


def sdf_to_mols(sdf_file: str, subsample: int) -> tuple[list, list]:
    """Read molecules from SDF file with subsampling.

    Args:
        sdf_file: Path to SDF file
        subsample: Maximum number of molecules to read

    Returns
    -------
        Tuple of (molecules list, SMILES list)
    """
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


def dropna(mols: list, smiles: list) -> tuple[list, list]:
    """Remove None values from molecule and SMILES lists.

    Args:
        mols: List of molecule objects
        smiles: List of SMILES strings

    Returns
    -------
        Tuple of (cleaned mols, cleaned smiles)
    """
    df = pd.DataFrame({"mols": mols, "smiles": smiles}).dropna()
    return df["mols"].tolist(), df["smiles"].tolist()


def format_number(x: float, pos: Any = None) -> str:
    """Format number for display.

    Args:
        x: Number to format
        pos: Position (for matplotlib FuncFormatter compatibility)

    Returns
    -------
        Formatted number string
    """
    if x >= 1e6:
        return f"{x/1e6:.1f}M"
    if x >= 1e3:
        return f"{x/1e3:.1f}K"
    return f"{x:.0f}"


def get_model_colors(
    model_names: list[str], cmap: str | None = None
) -> dict[str, Any]:
    """Generate color map for models.

    Args:
        model_names: List of model names
        cmap: Optional colormap name

    Returns
    -------
        Dictionary mapping model names to colors
    """
    if cmap is None:
        colors = plt.cm.YlOrRd(np.linspace(1, 0, len(model_names) + 1))
    else:
        colors = plt.colormaps.get_cmap(cmap)(
            np.linspace(1, 0, len(model_names) + 1)
        )

    return dict(zip(model_names, colors, strict=False))


def clean_name(name: str) -> str:
    """Clean metric names for display.

    Args:
        name: Name to clean

    Returns
    -------
        Cleaned name
    """
    for pattern in ["metrics", "_", ".csv"]:
        name = name.replace(pattern, "")
    return name.strip()


def filter_alerts(config: dict[str, Any]) -> pd.DataFrame:
    """Filter structural alerts based on configuration.

    Args:
        config: Configuration dictionary

    Returns
    -------
        Filtered alerts dataframe
    """
    df = pd.read_csv(config["alerts_data_path"])
    mask = df["rule_set_name"].isin(config["include_rulesets"])

    for ruleset in config["exclude_descriptions"]:
        is_other_ruleset = df["rule_set_name"] != ruleset
        is_not_excluded = ~df["description"].isin(
            config["exclude_descriptions"][ruleset]
        )
        mask &= is_other_ruleset | is_not_excluded

    return df[mask]


def common_postprocessing_statistics(
    filter_results: pd.DataFrame,
    res_df: pd.DataFrame,
    stat: pd.DataFrame | None,
    extend: pd.DataFrame | None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Combine filter results with statistics.

    Args:
        filter_results: Filter results dataframe
        res_df: Results dataframe
        stat: Optional statistics dataframe
        extend: Optional extended dataframe

    Returns
    -------
        Tuple of (results dataframe, extended dataframe)
    """
    if stat is not None:
        res_df = pd.concat([stat, res_df])

    filter_results = filter_results.drop(columns="mol")
    if extend is not None:
        filter_extended = pd.concat([extend, filter_results], ignore_index=True)
    else:
        filter_extended = filter_results.copy()

    return res_df, filter_extended


def process_one_file(
    config: dict[str, Any],
    input_path: str,
    apply_filter: Any,
    subsample: int,
) -> pd.DataFrame | None:
    """Process a single input file through filters.

    Args:
        config: Configuration dictionary
        input_path: Path to input file
        apply_filter: Filter function to apply
        subsample: Maximum number of molecules to process

    Returns
    -------
        Filter results dataframe or None if no molecules
    """
    input_type = input_path[input_path.rfind(".")+1:]
    assert input_type in {"csv", "smi", "sdf", "txt"}

    if input_type == "csv":
        data = pd.read_csv(input_path)
        smiles_col = "smiles"

        is_multi = data["model_name"].nunique(dropna=True) > 1
        if is_multi:
            smiles_str = data[smiles_col].tolist()
            model_names = data["model_name"].tolist()
            mols = [dm.to_mol(x) for x in smiles_str]
            mol_indices = data["mol_idx"].tolist()
            smiles = list(
                zip(smiles_str, model_names, mols, mol_indices, strict=False)
            )
        else:
            smiles = data[smiles_col].tolist()
            mols = [dm.to_mol(x) for x in smiles]
            mol_indices = data["mol_idx"].tolist()
            smiles = list(
                zip(smiles, [None] * len(smiles), mols, mol_indices, strict=False)
            )

    elif input_type in {"smi", "txt"}:
        input_file = Path(input_path)
        with input_file.open() as file:
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
            final_result = apply_filter(config, mols, smiles)
        else:
            final_result = apply_filter(config, mols)

    return final_result


def add_model_name_col(
    final_result: pd.DataFrame, smiles_with_model: list[tuple]
) -> pd.DataFrame:
    """Add smiles, model_name, and mol_idx columns from input data.

    Args:
        final_result: Results dataframe to update
        smiles_with_model: List of tuples with (smiles, model_name, mol, mol_idx)

    Returns
    -------
        Updated dataframe with identity columns
    """
    smiles_vals = [item[0] for item in smiles_with_model]
    model_vals = [
        item[1] if item[1] is not None else "single"
        for item in smiles_with_model
    ]
    mol_idx_vals = [item[3] if len(item) >= 4 else None for item in smiles_with_model]

    final_result["smiles"] = smiles_vals
    final_result["model_name"] = model_vals
    final_result["mol_idx"] = mol_idx_vals

    return final_result


def filter_function_applier(filter_name: str) -> Any:
    """Get filter function by name.

    Args:
        filter_name: Name of the filter

    Returns
    -------
        Filter function

    Raises
    ------
        ValueError: If filter name is not recognized
    """
    filters = {
        "common_alerts": apply_structural_alerts,
        "molgraph_stats": apply_molgraph_stats,
        "molcomplexity": apply_molcomplexity_filters,
        "NIBR": apply_nibr_filter,
        "bredt": apply_bredt_filter,
        "lilly": apply_lilly_filter,
    }
    if filter_name not in filters:
        msg = f"Filter {filter_name} not found"
        raise ValueError(msg)
    return filters[filter_name]


def apply_structural_alerts(
    config: dict[str, Any], mols: list, smiles_modelName_mols: list | None = None
) -> pd.DataFrame:
    """Apply structural alerts filters.

    Args:
        config: Configuration dictionary
        mols: List of molecule objects
        smiles_modelName_mols: Optional list of (smiles, model_name, mol) tuples

    Returns
    -------
        Dataframe with filter results
    """
    logger.info("Calculating Common Alerts...")
    def _apply_alerts(row):
        mol = row["mol"]
        row["smiles"] = dm.to_smiles(mol) if mol is not None else None

        config_structFilters = load_config(config["config_structFilters"])
        alert_data = filter_alerts(config_structFilters)
        df = alert_data.copy()
        df["matches"] = df.smarts.apply(lambda x, y: y.GetSubstructMatches(Chem.MolFromSmarts(x)), args=(mol,))
        grouped = df.groupby("rule_set_name").apply(
            lambda group: pd.Series({"matches": [match for matches in group["matches"]
                                                       for match in matches]
                                                       if any(matches
                                                              for matches in group["matches"])
                                                       else (),
                                     "reasons": ";".join(group[group["matches"].apply(lambda x: len(x) > 0)]["description"].fillna("").tolist())
                                                                if any(matches
                                                                       for matches in group["matches"])
                                                                else "",
                                    }),
            include_groups=False
        ).reset_index()
        grouped["pass_filter"] = grouped["matches"].apply(lambda x: bool(not x))
        for _, g_row in grouped.iterrows():
            name = g_row["rule_set_name"]
            row[f"pass_{name}"] = g_row["pass_filter"]
            row[f"reasons_{name}"] = g_row["reasons"]
        return row


    def _get_full_any_pass(row):
        pass_val = True
        any_pass_val  = False
        for col in row.index:
            if col.startswith("pass_") and col != "pass_any":
                pass_val &= row[col]
                any_pass_val |= row[col]
        row["pass"] = pass_val
        row["pass_any"] = any_pass_val
        return row


    logger.info("Processing %s filtered molecules", len(mols))

    n_jobs = config["n_jobs"]
    pandarallel.initialize(progress_bar=False, nb_workers=n_jobs, verbose=0)
    logger.info("Pandarallel initialized with %s workers", n_jobs)

    mols_df = pd.DataFrame({"mol" : mols})
    results = mols_df.parallel_apply(_apply_alerts, axis=1)
    results = results.parallel_apply(_get_full_any_pass, axis=1)

    if smiles_modelName_mols is not None:
        results = add_model_name_col(results, smiles_modelName_mols)
    return results


def apply_molgraph_stats(config, mols, smiles_modelName_mols=None):
    logger.info("Calculating Molecular Graph statistics...")
    severities = list(range(1, 12))

    results = {"mol" : mols}
    for s in severities:
        out = mc.functional.molecular_graph_filter(mols=mols,
                                                   max_severity=s,
                                                   n_jobs=-1,
                                                   progress=False,
                                                   return_idx=False,
                                                  )
        results[f"pass_{s}"] = out
    results = pd.DataFrame(results)

    if smiles_modelName_mols is not None:
        results = add_model_name_col(results, smiles_modelName_mols)
    return results


def apply_molcomplexity_filters(config, mols, smiles_modelName_mols=None):
    logger.info("Calculating Complexity filters...")
    final_result = pd.DataFrame({"mol" : mols,
                                 "pass" : True,
                                 "pass_any" : False
                               })
    alert_names = mc.complexity.ComplexityFilter.list_default_available_filters()
    for name in alert_names:
        cfilter = mc.complexity.ComplexityFilter(complexity_metric=name)
        final_result[f"pass_{name}"] = final_result["mol"].apply(cfilter)
        final_result["pass"] = final_result["pass"] * final_result[f"pass_{name}"]
        final_result["pass_any"] = final_result["pass_any"] + final_result[f"pass_{name}"]

    if smiles_modelName_mols is not None:
        final_result = add_model_name_col(final_result, smiles_modelName_mols)
    return final_result


def apply_bredt_filter(config, mols, smiles_modelName_mols=None):
    logger.info("Calculating Bredt filter...")
    out = mc.functional.bredt_filter(mols=mols,
                                     n_jobs=-1,
                                     progress=False,
                                     return_idx=False,
                                    )
    results = pd.DataFrame({"mol" : mols,
                            "pass" : out
                            })
    if smiles_modelName_mols is not None:
        results = add_model_name_col(results, smiles_modelName_mols)
    return results


def apply_nibr_filter(config, mols, smiles_modelName_mols=None):
    logger.info("Calculating NIBR filter...")
    n_jobs = config["n_jobs"]

    config_structFilters = load_config(config["config_structFilters"])
    scheduler = config_structFilters["nibr_scheduler"]

    nibr_filters = mc.structural.NIBRFilters()
    results = nibr_filters(mols=mols,
                           n_jobs=n_jobs,
                           scheduler=scheduler,
                           keep_details=True,
                        )
    if smiles_modelName_mols is not None:
        results = add_model_name_col(results, smiles_modelName_mols)
    return results


def apply_lilly_filter(config, mols, smiles_modelName_mols=None):
    if not LILLY_AVAILABLE:
        msg = (
            "Lilly demerits filter is not available. "
            "This filter requires conda/mamba-installed binaries. "
            "Install with: conda install lilly-medchem-rules\n"
            "Or disable this filter by setting 'calculate_lilly: False' "
            "in config_structFilters.yml"
        )
        raise ImportError(
            msg
        )

    logger.info("Calculating Lilly filter...")
    n_jobs = config["n_jobs"]

    config_strcuFilters = load_config(config["config_structFilters"])
    scheduler = config_strcuFilters["lilly_scheduler"]

    dfilter = LillyDemeritsFilters()
    results = dfilter(mols=mols,
                      n_jobs=n_jobs,
                      scheduler=scheduler,
                      )

    if smiles_modelName_mols is not None:
        results = add_model_name_col(results, smiles_modelName_mols)
    return results


def get_basic_stats(config, filter_results, model_name, filter_name, stat=None, extend=None):
    is_multi = filter_results["model_name"].nunique(dropna=True) > 1
    if is_multi:
        all_res = []
        all_extended = []
        for model, group in filter_results.groupby("model_name"):
            res_df, filter_extended = get_basic_stats(config, group.copy(), model, filter_name, stat, extend)
            all_res.append(res_df)
            all_extended.append(filter_extended)

        res_df = pd.concat(all_res, ignore_index=True)
        filter_extended = pd.concat(all_extended, ignore_index=True)
        return res_df, filter_extended

    num_mol = len(filter_results)
    filter_results = filter_results.dropna(subset="mol")
    filter_results["model_name"] = model_name

    if filter_name == "common_alerts":
        any_banned_percent = (~filter_results["pass"]).mean()
        all_banned_percent = (~filter_results["pass_any"]).mean()

        res_df = pd.DataFrame({"model_name" : [model_name],
                               "num_mol" : [num_mol],
                               "all_banned_ratio" : [all_banned_percent],
                               "any_banned_ratio" : [any_banned_percent]
                             })

        for name in config["include_rulesets"]:
            res_df[f"{name}_banned_ratio"] = 1 - filter_results[f"pass_{name}"].mean()

        res_df, filter_extended = common_postprocessing_statistics(filter_results, res_df, stat, extend)
        return res_df, filter_extended

    if filter_name == "molgraph_stats":
        res_df = pd.DataFrame({"model_name" : [model_name],
                               "num_mol" : [num_mol],
                              })
        for i in range(1, 12):
            res_df[f"banned_ratio_s_{i}"] = 1 - filter_results[f"pass_{i}"].mean()

        res_df, filter_extended = common_postprocessing_statistics(filter_results, res_df, stat, extend)
        pass_cols = [col for col in filter_extended.columns if col.startswith("pass_") and col != "pass_any"]
        filter_extended["pass"] = filter_extended[pass_cols].all(axis=1)
        return res_df, filter_extended


    if filter_name == "molcomplexity":
        any_banned_percent = 1 - filter_results["pass"].mean()
        all_banned_percent = 1 - filter_results["pass_any"].mean()

        res_df = pd.DataFrame({"model_name" : [model_name],
                                "num_mol" : [num_mol],
                                "all_banned_ratio" : [all_banned_percent],
                                "any_banned_ratio" : [any_banned_percent]
                             })
        alert_names = mc.complexity.ComplexityFilter.list_default_available_filters()
        for name in alert_names:
            res_df[f"{name}_banned_ratio"] = 1 - filter_results[f"pass_{name}"].mean()

        res_df, filter_extended = common_postprocessing_statistics(filter_results, res_df, stat, extend)
        return res_df, filter_extended


    if filter_name == "bredt":
        res_df = pd.DataFrame({"model_name" : [model_name],
                               "num_mol" : [num_mol],
                              })
        res_df["banned_ratio"] = 1 - filter_results["pass"].mean()

        res_df, filter_extended = common_postprocessing_statistics(filter_results, res_df, stat, extend)
        return res_df, filter_extended


    if filter_name == "NIBR":
        mean_severity = filter_results.severity.mean()
        max_severity = filter_results.severity.max()

        mean_n_covalent_motif = filter_results.n_covalent_motif.mean()
        mean_nonzero_special_mol = (filter_results.special_mol > 0).mean()

        res_df = pd.DataFrame({"model_name" : [model_name],
                               "num_mol" : [num_mol],
                               "mean_severity" : [mean_severity],
                               "max_severity" : [max_severity],
                               "mean_n_covalent_motif" : [mean_n_covalent_motif],
                               "mean_nonzero_special_mol" : [mean_nonzero_special_mol]
                             })

        pass_col = None
        if "pass" in filter_results.columns:
            pass_col = "pass"
        elif "pass_filter" in filter_results.columns:
            pass_col = "pass_filter"
        else:
            filter_results["pass"] = (filter_results["severity"] == 0)
            pass_col = "pass"

        res_df["banned_ratio"] = 1 - filter_results[pass_col].mean()

        res_df, filter_extended = common_postprocessing_statistics(filter_results, res_df, stat, extend)
        if "pass" not in filter_extended.columns:
            if "pass_filter" in filter_extended.columns:
                filter_extended = filter_extended.rename(columns={"pass_filter" : "pass"})
            elif pass_col in filter_extended.columns and pass_col != "pass":
                filter_extended = filter_extended.rename(columns={pass_col : "pass"})
            elif "severity" in filter_extended.columns:
                filter_extended["pass"] = (filter_extended["severity"] == 0)
        return res_df, filter_extended


    if filter_name == "lilly":
        mean_noNA_demerit_score = filter_results.demerit_score.dropna().mean()

        res_df = pd.DataFrame({"model_name": [model_name],
                               "num_mol": [num_mol],
                               "mean_noNA_demerit_score": mean_noNA_demerit_score
                             })

        pass_col = None
        if "pass" in filter_results.columns:
            pass_col = "pass"
        elif "pass_filter" in filter_results.columns:
            pass_col = "pass_filter"
        else:
            filter_results["pass"] = (filter_results["demerit_score"] == 0)
            pass_col = "pass"

        res_df["banned_ratio"] = 1 - filter_results[pass_col].mean()

        res_df, filter_extended = common_postprocessing_statistics(filter_results, res_df, stat, extend)
        if "pass" not in filter_extended.columns:
            if pass_col in filter_extended.columns:
                if pass_col != "pass":
                    filter_extended = filter_extended.rename(columns={pass_col : "pass"})
            elif "pass_filter" in filter_extended.columns:
                filter_extended = filter_extended.rename(columns={"pass_filter" : "pass"})
            elif "demerit_score" in filter_extended.columns:
                filter_extended["pass"] = (filter_extended["demerit_score"] == 0)
            elif "pass" in filter_results.columns and len(filter_extended) == len(filter_results):
                filter_extended["pass"] = filter_results["pass"].values
        return res_df, filter_extended

    msg = f"Filter {filter_name} not found"
    raise ValueError(msg)


def check_paths(config, paths) -> bool:
    all_filters = {}
    for k, v in config.items():
        if "calculate_" in k:
            k = k.replace("calculate_", "")
            all_filters[k] = v

    required_patterns = ["".join(k.split("_")) for k, v in all_filters.items() if v]
    missing_patterns = [pattern
                        for pattern in required_patterns
                        if not any(pattern.lower() in path.lower()
                            for path in paths)]
    if len(missing_patterns) > 0:
        msg = f"Invalid filter name(s) missing: {', '.join(missing_patterns)}"
        raise AssertionError(msg)
    return True


def plot_calculated_stats(config, prefix) -> None:
    folder_to_save = process_path(config["folder_to_save"])

    config_structFilters = load_config(config["config_structFilters"])

    if prefix == "beforeDescriptors":
        folder_path = Path(folder_to_save + f"{prefix}_StructFilters/")
        paths = list(folder_path.glob("*metrics.csv"))
    else:
        folder_path = Path(folder_to_save + "StructFilters/")
        paths = list(folder_path.glob("*metrics.csv"))
    check_paths(config_structFilters, paths)

    datas = []
    filter_names = []
    all_model_names = set()

    for path in paths:
        data = pd.read_csv(path)

        all_model_names.update(data["model_name"].dropna().unique())
        data = data.set_index("model_name")

        banned_cols = [col for col in data.columns if "banned_ratio" in col]
        data_filtered = data[[*banned_cols, "num_mol"]].copy()
        for banned_col in banned_cols:
            data_filtered.loc[:, f"num_banned_{banned_col}"] = data_filtered[banned_col] * data_filtered["num_mol"]
        datas.append(data_filtered)

        filter_name = path.split("/")[-1].replace("_metrics.csv", "")
        filter_names.append(filter_name)

    model_name_set = sorted(all_model_names)

    filter_results = {}
    subdir = f"{prefix}_StructFilters" if prefix == "beforeDescriptors" else "StructFilters"
    subdir_path = Path(folder_to_save + f"{subdir}/")
    filters_to_find = list(subdir_path.glob("*filteredMols.csv"))

    for path in filters_to_find:
        try:
            filter_data = pd.read_csv(path)
            filter_name = path.split("/")[-1].split("filteredMols.csv")[0].strip("_")

            num_passed_by_model = None
            if "pass" in filter_data.columns:
                passed = filter_data[filter_data["pass"]]
                if len(passed) > 0:
                    num_passed_by_model = passed.groupby("model_name").size().to_dict()

            if num_passed_by_model is not None:
                filter_results[filter_name] = num_passed_by_model
            else:
                default_models = filter_data["model_name"].unique().tolist()
                filter_results[filter_name] = dict.fromkeys(default_models, 0)

        except (IndexError, FileNotFoundError) as e:
            logger.warning(f"Could not process {path}: {e}")
            filter_results[filter_name] = {}

    all_models = model_name_set

    for filter_name, values in filter_results.items():
        if len(values) != len(all_models):
            for model in all_models:
                if model not in values:
                    filter_results[filter_name][model] = 0
        filter_results[filter_name] = dict(sorted(filter_results[filter_name].items()))

    n_plots = len(datas)
    n_cols = 2
    n_rows = (n_plots + n_cols - 1) // n_cols

    plt.figure(figsize=(40, 5*n_rows))
    for idx, (data, filter_name) in enumerate(zip(datas, filter_names, strict=False)):
        ax = plt.subplot(n_rows, n_cols, idx + 1)
        models = data.index
        x = np.arange(len(models))
        width = 0.8
        total_mols = data["num_mol"].sum()
        total = ax.barh(x, data.loc[models, "num_mol"], width, label=f"Total Molecules ({format_number(total_mols)})", color="#E5E5E5", alpha=0.5)

        clean_filter_name = filter_name.split("/")[-1].lower()
        for known_filter in filter_results:
            if known_filter.lower() in clean_filter_name:
                for i, (model, passed) in enumerate(filter_results[known_filter].items()):
                    bar_center_x = data.loc[models, "num_mol"].values[0] / 2
                    bar_center_y = x[i]
                    model_total = data.loc[model, "num_mol"] if model in data.index else data["num_mol"].iloc[0]
                    if model_total != 0:
                        pct = (passed / model_total * 100)
                        text = f"Passed molecules: {passed} ({pct:.1f}%)"
                    else:
                        text = f"Passed molecules: {passed} (0%)"
                    ax.annotate(text, (bar_center_x, bar_center_y), ha="center", va="center", fontsize=12, color="black", fontweight="bold",
                                bbox={"facecolor": "white", "alpha": 0.7, "edgecolor": "none", "pad": 3}, zorder=1000)

        for i, model in enumerate(models):
            count = data.loc[model, "num_mol"]
            max_bar_width = data["num_mol"].max()
            text_x_position = max_bar_width
            ax.text(text_x_position, i, int(count),  va="center", ha="left", fontsize=12, color="black", fontweight="bold")

        banned_bars = []
        banned_percentages = []
        ratio_cols = [col for col in data.columns if "banned_ratio" in col and "num_banned" not in col]
        colors = get_model_colors(model_names=ratio_cols, cmap="Paired")
        for col, color in zip(ratio_cols, colors.values(), strict=False):
            num_banned_col = f"num_banned_{col}"
            ratio_name = col.replace("banned_ratio", "").strip("_")
            ratio_name = clean_name(ratio_name)

            banned_count = data[num_banned_col]
            total_banned = banned_count.sum()
            banned_percent = (total_banned / total_mols) * 100 if total_mols > 0 else 0
            banned_percentages.append(banned_percent)

            if banned_percent == 0.0:
                label = f"{ratio_name} (0%)"
            else:
                label = (
                    f"{ratio_name} ({format_number(total_banned)}, "
                    f"{banned_percent:.1f}%)"
                )

            bar = ax.barh(x, banned_count, width, label=label, color=color, alpha=0.8)
            banned_bars.append(bar)

        clean_filter_name = filter_name.split("/")[-1]
        ax.set_title(clean_name(clean_filter_name), fontsize=14, pad=20, fontweight="bold")
        ax.set_yticks(x)
        ax.set_yticklabels(models, fontsize=12)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(format_number))
        ax.set_xlim(left=0)

    ax.set_xlabel(f"Number of Molecules (Total: {format_number(total_mols)})", fontsize=12, labelpad=10)
    ax.set_ylabel("Models", fontsize=12, labelpad=10)

    ax.grid(True, axis="x", alpha=0.2, linestyle="--")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    sorted_indices = np.argsort(banned_percentages)[::-1]
    sorted_handles = [total] + [banned_bars[i] for i in sorted_indices]

    legend = ax.legend(handles=sorted_handles, loc="center left",  bbox_to_anchor=(1.02, 0.5), fontsize=11, ncol=1)
    legend.get_frame().set_alpha(0.9)
    legend.get_frame().set_edgecolor("lightgray")

    plt.subplots_adjust(right=0.85, hspace=0.6, wspace=0.5)

    subdir = f"{prefix}_StructFilters" if prefix == "beforeDescriptors" else "StructFilters"
    plt.savefig(folder_to_save + f"{subdir}/MoleculeCountsComparison.png", dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")
    plt.close()


def plot_restriction_ratios(config, prefix) -> None:
    folder_to_save = process_path(config["folder_to_save"])
    folder_name = config["folder_to_save"].split("/")[-1]

    config_structFilters = load_config(config["config_structFilters"])

    if prefix == "beforeDescriptors":
        folder_path = Path(folder_to_save + f"{prefix}_StructFilters/")
        paths = list(folder_path.glob("*metrics.csv"))
    else:
        folder_path = Path(folder_to_save + "StructFilters/")
        paths = list(folder_path.glob("*metrics.csv"))
    check_paths(config_structFilters, paths)

    if not paths:
        logger.error(f"No data files found in {folder_to_save}")
        return

    filter_data = {}
    model_names_filters = {}
    for path in paths:
        filter_name = path.split(f"{folder_name}/")[-1].split("_metrics.csv")[0]
        data = pd.read_csv(path)

        ratio_cols = [col for col in data.columns if "banned_ratio" in col]
        model_names_filters[filter_name] = dict(zip(data["model_name"].tolist(), data["num_mol"].tolist(), strict=False))

        if not ratio_cols:
            continue

        clean_cols = {col: col.replace("_banned_ratio", "").replace("banned_ratio", "").replace("_s", "s") for col in ratio_cols}
        ratios = data[ratio_cols].rename(columns=clean_cols)
        actual_model_names = data["model_name"].tolist()
        ratios.index = actual_model_names

        row = ratios.iloc[0]
        if row.isna().all():
            continue

        all_value = None
        if "all" in row.index:
            all_value = row["all"]
            row = row.drop("all")

        sorted_values = row.sort_values(ascending=False)

        if all_value is not None:
            if all_value >= sorted_values.iloc[0]:
                sorted_index = pd.Index(["all"]).append(sorted_values.index)
            else:
                sorted_index = sorted_values.index.append(pd.Index(["all"]))
            ratios = ratios[sorted_index]
        else:
            ratios = ratios[sorted_values.index]

        filter_data[filter_name] = ratios
    if not filter_data:
        logger.error("No valid data to plot")
        return

    model_names_filters = pd.DataFrame(model_names_filters).reset_index()
    model_names_filters = model_names_filters.rename(columns={"index": "model_name"})

    plt.style.use("default")
    sns.set_style("white")
    sns.set_context("talk")

    n_filters = len(filter_data)
    n_cols = min(2, n_filters)
    n_rows = (n_filters + n_cols - 1) // n_cols

    fig = plt.figure(figsize=(16, 7*n_rows))
    fig.suptitle("Comparison of Restriction Ratios Across Different Filters", fontsize=16, y=0.98, fontweight="bold")

    for idx, (filter_name, data) in enumerate(filter_data.items()):
        number_of_mols = np.array(model_names_filters[filter_name].tolist())
        ax = plt.subplot(n_rows, n_cols, idx + 1)
        for col in data.columns:
            if col not in ["any", "all", "model_name", "num_mol"]:
                data[col] = number_of_mols * (1 - np.array(data[col].tolist()))
        if not data.empty and data.notna().any().any():
            if "all" in data.columns:
                data = data.drop(columns=["all"])
            if "any" in data.columns:
                data = data.drop(columns=["any"])

            from matplotlib.colors import LinearSegmentedColormap
            custom_cmap = LinearSegmentedColormap.from_list("custom", ["white", "#B29EEE"])
            if idx == 1:

                sns.heatmap(data.T, cmap=custom_cmap, cbar_kws={"label": "Passed Molecules", "format": "%d" }, ax=ax, vmin=0, vmax=max(data.max()),
                            fmt=".0f", annot=True, annot_kws={"size": 12, "rotation": 0, "color": "black"}, cbar=True)
            else:
                sns.heatmap(data.T, cmap=custom_cmap, cbar_kws={"label": "Passed Molecules", "format": "%d" }, ax=ax, vmin=0, vmax=max(data.max()),
                            fmt=".0f", annot=True, annot_kws={"size": 12, "rotation": 0, "color": "black"}, cbar=False)

            ax.set_title(f"{clean_name(filter_name)} Filter", fontsize=12, fontweight="bold")
            plt.setp(ax.get_yticklabels(), rotation=0, ha="right", fontsize=12)
            plt.setp(ax.get_xticklabels(), rotation=0, ha="right", fontsize=12)
            ax.set_xlabel("Model")

            actual_model_names = data.index.tolist()
            if len(actual_model_names) == len(ax.get_xticklabels()):
                ax.set_xticklabels(actual_model_names)

        else:
            ax.text(0.5, 0.5, "No data available", horizontalalignment="center", verticalalignment="center", transform=ax.transAxes)
            ax.set_title(f"{clean_name(filter_name)} Filter", pad=10, fontsize=11, fontweight="bold")

    plt.tight_layout()
    subdir = f"{prefix}_StructFilters" if prefix == "beforeDescriptors" else "StructFilters"
    plt.savefig(folder_to_save + f"{subdir}/RestrictionRatiosComparison.png", dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")
    plt.close()


def filter_data(
    config: dict[str, Any], prefix: str
) -> pd.DataFrame:
    """Filter and merge molecular data from different filters.

    Args:
        config: Configuration dictionary
        prefix: Stage prefix

    Returns
    -------
        Filtered dataframe
    """
    if prefix == "beforeDescriptors":
        subdir = f"{prefix}_StructFilters"
    else:
        subdir = "StructFilters"
    base_folder = process_path(config["folder_to_save"])
    folder_to_save = process_path(config["folder_to_save"], key_word=subdir)

    folder_path = Path(folder_to_save)
    paths = list(folder_path.glob("*filteredMols.csv"))

    columns_to_drop = ["pass", "any_pass", "name", "pass_any"]
    datas = []
    for path in paths:
        data = pd.read_csv(path)
        for col in columns_to_drop:
            if col in data.columns:
                data = data.drop(columns=[col])
        datas.append(data)

    if len(datas) > 0:
        filtered_data = datas[0].copy()

        for df in datas[1:]:
            merge_cols = ["smiles", "model_name"]
            existing_cols = set(filtered_data.columns) - set(merge_cols)
            new_cols = [col for col in df.columns if col not in existing_cols and col not in merge_cols]

            if new_cols:
                cols_to_merge = merge_cols + new_cols
                filtered_data = filtered_data.merge(df[cols_to_merge], on=merge_cols, how="inner")
            else:
                filtered_data = filtered_data.merge(df[merge_cols], on=merge_cols, how="inner")
    else:
        filtered_data = pd.DataFrame(columns=["smiles", "model_name", "mol_idx"])

    if "mol_idx" not in filtered_data.columns:
        filtered_data["mol_idx"] = None

    cols = ["smiles", "model_name", "mol_idx"]
    out_df = filtered_data[cols].copy()
    out_df.to_csv(folder_to_save + "passStructFiltersSMILES.csv", index=False)

    if prefix != "beforeDescriptors":
        descriptors_path = Path(base_folder + "Descriptors/passDescriptorsSMILES.csv")
        if descriptors_path.exists():
            input_path = str(descriptors_path)
        else:
            sampled_path = Path(base_folder + "sampledMols.csv")
            if sampled_path.exists():
                input_path = str(sampled_path)
            else:
                try:
                    from hedge.stages.structFilters.main import (
                        _get_input_path,
                    )

                    input_path = _get_input_path(config, prefix, base_folder)
                except Exception:
                    logger.debug(
                        "Could not determine input path from main module"
                    )
                    input_path = None
    else:
        sampled_path = Path(base_folder + "sampledMols.csv")
        if sampled_path.exists():
            input_path = str(sampled_path)
        else:
            try:
                from hedge.stages.structFilters.main import (
                    _get_input_path,
                )

                input_path = _get_input_path(config, prefix, base_folder)
            except Exception:
                logger.debug("Could not determine input path from main module")
                input_path = None

    input_path_obj = Path(input_path) if input_path else None
    if input_path_obj and input_path_obj.exists():
        try:
            all_input = pd.read_csv(input_path)
            if len(out_df) > 0:
                merge_cols = ["smiles", "model_name"]
                merged = all_input.merge(
                    out_df[merge_cols],
                    on=merge_cols,
                    how="left",
                    indicator=True,
                )
                fail_molecules = merged[merged["_merge"] == "left_only"].drop(
                    columns=["_merge"]
                )
            else:
                fail_molecules = all_input.copy()

            if len(fail_molecules) > 0:
                extended_paths = list(folder_path.glob("*extended.csv"))
                all_extended = None
                for ext_path in extended_paths:
                    try:
                        ext_df = pd.read_csv(ext_path)
                        if all_extended is None:
                            all_extended = ext_df.copy()
                        else:
                            merge_cols = ["smiles", "model_name"]
                            pass_cols = [
                                col
                                for col in ext_df.columns
                                if col.startswith("pass_") or col == "pass"
                            ]
                            if pass_cols:
                                cols_to_merge = merge_cols + pass_cols
                                all_extended = all_extended.merge(
                                    ext_df[cols_to_merge],
                                    on=merge_cols,
                                    how="outer",
                                    suffixes=("", "_dup"),
                                )
                                for col in pass_cols:
                                    if f"{col}_dup" in all_extended.columns:
                                        all_extended[col] = all_extended[
                                            col
                                        ].fillna(all_extended[f"{col}_dup"])
                                        all_extended = all_extended.drop(
                                            columns=[f"{col}_dup"]
                                        )
                    except Exception:
                        logger.debug(
                            f"Could not process extended file: {ext_path}"
                        )
                        continue

                if all_extended is not None:
                    merge_cols = ["smiles", "model_name"]
                    pass_cols = [
                        col
                        for col in all_extended.columns
                        if col.startswith("pass_") or col == "pass"
                    ]
                    if pass_cols:
                        cols_to_merge = merge_cols + pass_cols
                        fail_molecules = fail_molecules.merge(
                            all_extended[cols_to_merge],
                            on=merge_cols,
                            how="left",
                        )
                        for col in pass_cols:
                            if col in fail_molecules.columns:
                                fail_molecules[col] = fail_molecules[col].fillna(
                                    False
                                )

                id_cols = ["smiles", "model_name", "mol_idx"]
                pass_cols_final = [
                    col
                    for col in fail_molecules.columns
                    if col.startswith("pass_") or col == "pass"
                ]
                fail_cols = [
                    c for c in id_cols if c in fail_molecules.columns
                ] + pass_cols_final
                fail_molecules[fail_cols].to_csv(
                    folder_to_save + "failStructFiltersSMILES.csv", index=False
                )
        except Exception as e:
            logger.warning(f"Could not create failStructFiltersSMILES.csv: {e}")

    return filtered_data


def inject_identity_columns_to_all_csvs(
    config: dict[str, Any], prefix: str
) -> None:
    """Ensure identity columns are ordered consistently in all CSVs.

    Args:
        config: Configuration dictionary
        prefix: Stage prefix
    """
    if prefix == "beforeDescriptors":
        subdir = f"{prefix}_StructFilters"
    else:
        subdir = "StructFilters"
    target_folder = process_path(config["folder_to_save"], key_word=subdir)

    folder_path = Path(target_folder)
    csv_paths = list(folder_path.glob("*.csv"))
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
        except Exception:
            logger.debug(f"Could not process CSV file: {path}")
            continue


def analyze_filter_failures(
    file_path: str,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]] | tuple[None, None, None]:
    """Analyze filter failures from extended CSV file and generate visualizations.

    Args:
        file_path: Path to extended CSV file

    Returns
    -------
        Tuple of (filter_failures, filter_reasons, all_detailed_reasons)
        or (None, None, None) if no filter columns found
    """
    logger.debug(f"Analyzing filter failures from: {file_path}")
    df = pd.read_csv(file_path, low_memory=False)

    filter_columns = [col for col in df.columns if col.startswith("pass_") and col not in {"pass", "pass_any"}]

    if not filter_columns:
        return None, None, None

    filter_failures = {}
    filter_reasons = {}
    all_detailed_reasons = {}

    for col in filter_columns:
        filter_name = col.replace("pass_", "")

        failures = (~df[col]).sum()
        total = len(df)
        failure_percentage = (failures / total) * 100

        filter_failures[filter_name] = {"failures": failures,
                                        "total": total,
                                        "percentage": failure_percentage
                                       }

        reasons_col = f"reasons_{filter_name}"
        if reasons_col in df.columns:
            failed_molecules = df[~df[col]]
            reasons_data = failed_molecules[reasons_col].dropna()

            reason_counts = {}
            for reasons_str in reasons_data:
                if pd.notna(reasons_str) and str(reasons_str).strip():
                    individual_reasons = [r.strip() for r in str(reasons_str).split(";") if r.strip()]
                    for reason in individual_reasons:
                        reason_counts[reason] = reason_counts.get(reason, 0) + 1

            sorted_reasons = sorted(reason_counts.items(), key=lambda x: x[1], reverse=True)
            filter_reasons[filter_name] = sorted_reasons
            all_detailed_reasons[filter_name] = reason_counts

    _create_main_filter_plot(filter_failures, file_path)
    _create_individual_filter_plots(filter_failures, filter_reasons, file_path)
    _create_multi_panel_filter_plot(filter_failures, filter_reasons, file_path)

    all_reasons = {}
    for filter_name, reasons in filter_reasons.items():
        for reason, count in reasons:
            if reason in all_reasons:
                all_reasons[reason] += count
            else:
                all_reasons[reason] = count

    top_reasons = sorted(all_reasons.items(), key=lambda x: x[1], reverse=True)[:5]

    if top_reasons:
        logger.info(
            "Top 5 most common filter failure reasons "
            "(molecules may have multiple reasons):"
        )
        for i, (reason, count) in enumerate(top_reasons, 1):
            logger.info("  %d. %s: %s failures", i, reason, count)

    _create_complete_reasons_breakdown(all_detailed_reasons, filter_failures, file_path)
    _create_comprehensive_overview(filter_reasons, filter_failures, file_path)
    _create_summary_table(filter_failures, filter_reasons, file_path)

    return filter_failures, filter_reasons, all_detailed_reasons


def _create_main_filter_plot(
    filter_failures: dict[str, Any], file_path: str
) -> None:
    """Create main filter failures bar chart.

    Args:
        filter_failures: Dictionary of filter failure statistics
        file_path: Path to save output
    """
    plot_data = []
    for filter_name, stats in filter_failures.items():
        plot_data.append(
            {
                "filter": filter_name,
                "failures": stats["failures"],
                "percentage": stats["percentage"],
            }
        )

    plot_df = pd.DataFrame(plot_data)
    plot_df = plot_df.sort_values("failures", ascending=False)

    plt.figure(figsize=(max(16, len(plot_df) * 0.6), 16))
    plt.bar(
        range(len(plot_df)),
        plot_df["failures"],
        color="steelblue",
        alpha=0.8,
        width=0.3,
    )

    plt.xlabel("Filters", fontsize=20)
    plt.ylabel("Number of Molecules Failed", fontsize=20)
    plt.title(
        "Number of Molecules Failed by Each Filter",
        fontsize=26,
        fontweight="bold",
    )
    plt.xticks(
        range(len(plot_df)),
        plot_df["filter"],
        rotation=45,
        ha="right",
        fontsize=16,
    )

    for i, (_, row) in enumerate(plot_df.iterrows()):
        plt.text(
            i,
            row["failures"] + max(plot_df["failures"]) * 0.01,
            f"{row['failures']}\n({row['percentage']:.1f}%)",
            ha="center",
            va="bottom",
            fontsize=14,
        )

    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()

    file_path_obj = Path(file_path)
    path_to_save = file_path_obj.parent / "CommonAlertsBreakdown"
    path_to_save.mkdir(parents=True, exist_ok=True)
    output_path = path_to_save / "filter_failures_plot.png"
    plt.savefig(output_path, dpi=600, bbox_inches="tight")
    plt.close()


def _create_individual_filter_plots(filter_failures, filter_reasons, file_path) -> None:
    """Create individual plots for each filter showing failure reasons."""
    for filter_name, stats in filter_failures.items():
        if stats["failures"] > 0:
            reasons_data = filter_reasons.get(filter_name, [])

            if not reasons_data:
                continue

            plot_data = []
            for reason, count in reasons_data:
                plot_data.append({"Reason": reason,
                                  "Count": count,
                                  "Percentage_of_Filter_Failures": (count / stats["failures"]) * 100 if stats["failures"] > 0 else 0
                                })

            if not plot_data:
                continue

            plot_df = pd.DataFrame(plot_data)
            plot_df = plot_df.sort_values("Count", ascending=False)

            plt.figure(figsize=(max(16, len(plot_df) * 0.6), 20))
            bars = plt.bar(range(len(plot_df)), plot_df["Count"], color="steelblue", alpha=0.8, width=0.3)

            plt.xlabel("Failure Reasons", fontsize=20)
            plt.ylabel("Number of Molecules Failed", fontsize=20)
            title_text = (
                f'{filter_name.upper()} - Failure Reasons '
                f'({len(plot_df)} reasons, {stats["failures"]} total failures)'
            )
            plt.title(title_text, fontsize=26, fontweight="bold")
            plt.xticks(range(len(plot_df)), plot_df["Reason"], rotation=45, ha="right", fontsize=max(10, min(16, 300 // len(plot_df))))

            for i, (_bar, count) in enumerate(zip(bars, plot_df["Count"], strict=False)):
                plt.text(i, count + max(plot_df["Count"]) * 0.01, f"{count}\n({plot_df.iloc[i]['Percentage_of_Filter_Failures']:.1f}%)", ha="center", va="bottom", fontsize=12)

            plt.grid(axis="y", alpha=0.3)
            plt.tight_layout()

            file_path_obj = Path(file_path)
            path_to_save = file_path_obj.parent / "CommonAlertsBreakdown"
            path_to_save.mkdir(parents=True, exist_ok=True)
            output_path = path_to_save / f"{filter_name}_reasons_plot.png"
            plt.savefig(output_path, dpi=600, bbox_inches="tight")
            plt.close()


def _create_multi_panel_filter_plot(filter_failures, filter_reasons, file_path) -> None:
    """Create multi-panel plot showing all filters with reasons."""
    sorted_filters = sorted(filter_failures.items(), key=lambda x: x[1]["failures"], reverse=True)
    sorted_filters = [(name, stats) for name, stats in sorted_filters if stats["failures"] > 0]

    num_filters = len(sorted_filters)
    if num_filters == 0:
        return

    if num_filters <= 3:
        rows = 1
        cols = num_filters
    elif num_filters <= 6:
        rows = 2
        cols = 3
    elif num_filters <= 9:
        rows = 3
        cols = 3
    elif num_filters <= 12:
        rows = 3
        cols = 4
    else:
        cols = 4
        rows = (num_filters + cols - 1) // cols

    plt.figure(figsize=(cols * 6, rows * 6))

    for i, (filter_name, _stats) in enumerate(sorted_filters):
        all_reasons_data = filter_reasons.get(filter_name, [])

        plt.subplot(rows, cols, i + 1)
        reason_names = [r[0] for r in all_reasons_data]
        reason_counts = [r[1] for r in all_reasons_data]

        if len(reason_names) > 10:
            reason_names = reason_names[:10]
            reason_counts = reason_counts[:10]
            title_suffix = f"(Top 10 of {len(all_reasons_data)} reasons)"
        else:
            title_suffix = f"({len(all_reasons_data)} reasons)"

        plt.bar(range(len(reason_names)), reason_counts, color="steelblue", alpha=0.8, width=0.3)

        plt.xlabel("Reasons", fontsize=14)
        plt.ylabel("Molecules Failed", fontsize=14)
        plt.title(f"{filter_name.upper()}\n{title_suffix}", fontsize=16, fontweight="bold")

        truncated_names = []
        for name in reason_names:
            if len(name) > 15:
                truncated_names.append(name[:12] + "...")
            else:
                truncated_names.append(name)

        plt.xticks(range(len(truncated_names)), truncated_names, rotation=45, ha="right", fontsize=12)

        for j, count in enumerate(reason_counts):
            plt.text(j, count + max(reason_counts) * 0.01 if reason_counts else 0, f"{count}", ha="center", va="bottom", fontsize=11)

        plt.grid(axis="y", alpha=0.3)

    plt.tight_layout(h_pad=1.5, w_pad=1.0)

    file_path_obj = Path(file_path)
    path_to_save = file_path_obj.parent / "CommonAlertsBreakdown"
    path_to_save.mkdir(parents=True, exist_ok=True)
    output_path = path_to_save / "all_filters_reasons_plot.png"
    plt.savefig(output_path, dpi=600, bbox_inches="tight")
    plt.close()


def _create_complete_reasons_breakdown(
    all_detailed_reasons: dict[str, Any],
    filter_failures: dict[str, Any],
    file_path: str,
) -> pd.DataFrame:
    """Create complete CSV breakdown of all reasons.

    Args:
        all_detailed_reasons: Dictionary of detailed failure reasons
        filter_failures: Dictionary of filter failure statistics
        file_path: Path to save output

    Returns
    -------
        Breakdown dataframe
    """
    breakdown_data = []

    for filter_name, reasons_dict in all_detailed_reasons.items():
        total_failures = filter_failures[filter_name]["failures"]

        for reason, count in sorted(
            reasons_dict.items(), key=lambda x: x[1], reverse=True
        ):
            percentage = (count / total_failures) * 100 if total_failures > 0 else 0
            breakdown_data.append(
                {
                    "Ruleset": filter_name,
                    "Reason": reason,
                    "Count": count,
                    "Percentage_of_Filter_Failures": percentage,
                    "Total_Filter_Failures": total_failures,
                }
            )

    breakdown_df = pd.DataFrame(breakdown_data)

    file_path_obj = Path(file_path)
    path_to_save = file_path_obj.parent / "CommonAlertsBreakdown"
    path_to_save.mkdir(parents=True, exist_ok=True)
    output_path = path_to_save / "complete_reasons_breakdown.csv"
    breakdown_df.to_csv(output_path, index=False)

    return breakdown_df


def _create_comprehensive_overview(filter_reasons, filter_failures, file_path) -> None:
    """Create comprehensive overview plot of most common failure reasons."""
    all_reasons = {}

    for reasons in filter_reasons.values():
        for reason, count in reasons:
            if reason in all_reasons:
                all_reasons[reason] += count
            else:
                all_reasons[reason] = count
    top_reasons = sorted(all_reasons.items(), key=lambda x: x[1], reverse=True)

    if not top_reasons:
        return

    display_count = min(30, len(top_reasons))
    plt.figure(figsize=(max(16, display_count * 0.6), 16))

    reason_names = []
    reason_counts = []

    for reason, count in top_reasons[:display_count]:
        reason_short = reason[:27] + "..." if len(reason) > 30 else reason
        reason_names.append(reason_short)
        reason_counts.append(count)

    bars = plt.bar(range(len(reason_names)), reason_counts, color="darkgreen", alpha=0.7, width=0.3)

    plt.xlabel("Failure Reasons", fontsize=20)
    plt.ylabel("Total Number of Molecules Failed", fontsize=20)
    title_text = (
        f"Most Common Molecular Filter Failure Reasons "
        f"(Top {display_count} of {len(top_reasons)})"
    )
    plt.title(title_text, fontsize=26, fontweight="bold")
    plt.xticks(range(len(reason_names)), reason_names, rotation=45, ha="right", fontsize=16)

    for i, (_bar, count) in enumerate(zip(bars, reason_counts, strict=False)):
        plt.text(i, count + max(reason_counts) * 0.01, f"{count}", ha="center", va="bottom", fontsize=14)

    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()

    file_path_obj = Path(file_path)
    path_to_save = file_path_obj.parent / "CommonAlertsBreakdown"
    path_to_save.mkdir(parents=True, exist_ok=True)
    output_path = path_to_save / "comprehensive_reasons_overview.png"
    plt.savefig(output_path, dpi=600, bbox_inches="tight")

    plt.close()

    all_reasons_df = pd.DataFrame(top_reasons, columns=["Reason", "Total_Count"])
    path_to_save = file_path_obj.parent / "CommonAlertsBreakdown"
    path_to_save.mkdir(parents=True, exist_ok=True)
    output_path = path_to_save / "all_reasons_summary.csv"
    all_reasons_df.to_csv(output_path, index=False)


def _create_summary_table(
    filter_failures: dict[str, Any],
    filter_reasons: dict[str, Any],
    file_path: str,
) -> pd.DataFrame:
    """Create summary table CSV with filter statistics.

    Args:
        filter_failures: Dictionary of filter failure statistics
        filter_reasons: Dictionary of filter reasons
        file_path: Path to save output

    Returns
    -------
        Summary dataframe
    """
    summary_data = []

    for filter_name, stats in filter_failures.items():
        row = {"Ruleset": filter_name,
               "Total_Failures": stats["failures"],
               "Failure_Percentage": stats["percentage"],
               "Total_Molecules": stats["total"],
               "Unique_Reasons_Count": len(filter_reasons.get(filter_name, []))
              }

        if filter_reasons.get(filter_name):
            for i, (reason, count) in enumerate(filter_reasons[filter_name][:5], 1):
                row[f"Top_Reason_{i}"] = reason
                row[f"Top_Reason_{i}_Count"] = count
                if stats["failures"] > 0:
                    row[f"Top_Reason_{i}_Percentage"] = (
                        count / stats["failures"]
                    ) * 100
                else:
                    row[f"Top_Reason_{i}_Percentage"] = 0

        summary_data.append(row)

    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df.sort_values("Total_Failures", ascending=False)

    file_path_obj = Path(file_path)
    path_to_save = file_path_obj.parent / "CommonAlertsBreakdown"
    path_to_save.mkdir(parents=True, exist_ok=True)
    output_path = path_to_save / "filter_summary_table.csv"
    summary_df.to_csv(output_path, index=False)
    logger.debug(f"Summary table saved to: {output_path}")

    return summary_df


def plot_filter_failures_analysis(
    config: dict[str, Any], prefix: str
) -> None:
    """Analyze and plot filter failures for extended CSV files.

    Args:
        config: Configuration dictionary
        prefix: Stage prefix
    """
    folder_to_save = process_path(config["folder_to_save"])
    if prefix == "beforeDescriptors":
        subfolder = f"{prefix}_StructFilters"
    else:
        subfolder = "StructFilters"
    struct_folder_path = Path(folder_to_save + f"{subfolder}/")

    if not struct_folder_path.exists():
        return

    extended_files = list(struct_folder_path.glob("*_extended.csv"))

    if not extended_files:
        logger.debug("No extended CSV files found for failure analysis")
        return

    for file_path in extended_files:
        try:
            analyze_filter_failures(str(file_path))
        except Exception as e:
            logger.debug(f"Error analyzing filter failures for {file_path}: {e}")
