from pathlib import Path

import pandas as pd

from hedgehog.configs.logger import logger
from hedgehog.stages.structFilters.filters import mc
from hedgehog.utils.paths import process_path


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


def _get_pass_column(df, fallback_column=None, fallback_condition=None):
    """Determine the pass column name, creating it if necessary.

    Args:
        df: DataFrame to check
        fallback_column: Column name to use for creating pass column if not found
        fallback_condition: Lambda function to apply to fallback_column

    Returns:
        Name of the pass column
    """
    if "pass" in df.columns:
        return "pass"
    if "pass_filter" in df.columns:
        return "pass_filter"
    if fallback_column and fallback_column in df.columns and fallback_condition:
        df["pass"] = fallback_condition(df[fallback_column])
        return "pass"
    return None


def _ensure_pass_column_in_extended(
    filter_extended, pass_col, filter_results, fallback_column=None
):
    """Ensure the 'pass' column exists in filter_extended DataFrame."""
    if "pass" in filter_extended.columns:
        return filter_extended

    if pass_col and pass_col in filter_extended.columns and pass_col != "pass":
        filter_extended.rename(columns={pass_col: "pass"}, inplace=True)
    elif "pass_filter" in filter_extended.columns:
        filter_extended.rename(columns={"pass_filter": "pass"}, inplace=True)
    elif fallback_column and fallback_column in filter_extended.columns:
        filter_extended["pass"] = filter_extended[fallback_column] == 0
    elif "pass" in filter_results.columns and len(filter_extended) == len(
        filter_results
    ):
        filter_extended["pass"] = filter_results["pass"].values

    return filter_extended


def _create_base_stats_df(model_name, num_mol, **extra_columns):
    """Create base statistics DataFrame with model_name and num_mol."""
    data = {"model_name": [model_name], "num_mol": [num_mol]}
    data.update({k: [v] for k, v in extra_columns.items()})
    return pd.DataFrame(data)


def get_basic_stats(
    config, filter_results, model_name, filter_name, stat=None, extend=None
):
    """Calculate basic statistics for filter results."""
    model_col = (
        filter_results["model_name"] if "model_name" in filter_results.columns else None
    )
    is_multi = (not isinstance(model_name, str)) or (
        model_col is not None and model_col.nunique(dropna=True) > 1
    )
    if is_multi:
        if model_col is None:
            raise ValueError(
                "Multi-model statistics requested but 'model_name' column is missing"
            )
        all_res = []
        all_extended = []
        for model, group in filter_results.groupby("model_name"):
            res_df, filter_extended = get_basic_stats(
                config, group.copy(), model, filter_name, stat, extend
            )
            all_res.append(res_df)
            all_extended.append(filter_extended)
        return pd.concat(all_res, ignore_index=True), pd.concat(
            all_extended, ignore_index=True
        )

    num_mol = len(filter_results)
    filter_results.dropna(subset="mol", inplace=True)
    if isinstance(model_name, str):
        filter_results["model_name"] = model_name

    if filter_name == "common_alerts":
        res_df = _create_base_stats_df(
            model_name,
            num_mol,
            all_banned_ratio=(~filter_results["pass"]).mean(),
            any_banned_ratio=(~filter_results["pass_any"]).mean(),
        )
        for name in config["include_rulesets"]:
            res_df[f"{name}_banned_ratio"] = 1 - filter_results[f"pass_{name}"].mean()
        return common_postprocessing_statistics(filter_results, res_df, stat, extend)

    if filter_name == "molgraph_stats":
        res_df = _create_base_stats_df(model_name, num_mol)
        for i in range(1, 12):
            res_df[f"banned_ratio_s_{i}"] = 1 - filter_results[f"pass_{i}"].mean()
        res_df, filter_extended = common_postprocessing_statistics(
            filter_results, res_df, stat, extend
        )
        pass_cols = [
            col
            for col in filter_extended.columns
            if col.startswith("pass_") and col != "pass_any"
        ]
        filter_extended["pass"] = filter_extended[pass_cols].all(axis=1)
        return res_df, filter_extended

    if filter_name == "molcomplexity":
        res_df = _create_base_stats_df(
            model_name,
            num_mol,
            all_banned_ratio=1 - filter_results["pass"].mean(),
            any_banned_ratio=1 - filter_results["pass_any"].mean(),
        )
        for name in mc.complexity.ComplexityFilter.list_default_available_filters():
            res_df[f"{name}_banned_ratio"] = 1 - filter_results[f"pass_{name}"].mean()
        return common_postprocessing_statistics(filter_results, res_df, stat, extend)

    if filter_name in (
        "bredt",
        "protecting_groups",
        "ring_infraction",
        "stereo_center",
        "halogenicity",
        "symmetry",
    ):
        res_df = _create_base_stats_df(model_name, num_mol)
        res_df["banned_ratio"] = 1 - filter_results["pass"].mean()
        return common_postprocessing_statistics(filter_results, res_df, stat, extend)

    if filter_name == "NIBR":
        res_df = _create_base_stats_df(
            model_name,
            num_mol,
            mean_severity=filter_results.severity.mean(),
            max_severity=filter_results.severity.max(),
            mean_n_covalent_motif=filter_results.n_covalent_motif.mean(),
            mean_nonzero_special_mol=(filter_results.special_mol > 0).mean(),
        )
        pass_col = _get_pass_column(filter_results, "severity", lambda x: x == 0)
        res_df["banned_ratio"] = 1 - filter_results[pass_col].mean()
        res_df, filter_extended = common_postprocessing_statistics(
            filter_results, res_df, stat, extend
        )
        filter_extended = _ensure_pass_column_in_extended(
            filter_extended, pass_col, filter_results, "severity"
        )
        return res_df, filter_extended

    if filter_name == "lilly":
        res_df = _create_base_stats_df(
            model_name,
            num_mol,
            mean_noNA_demerit_score=filter_results.demerit_score.dropna().mean(),
        )
        pass_col = _get_pass_column(filter_results, "demerit_score", lambda x: x == 0)
        res_df["banned_ratio"] = 1 - filter_results[pass_col].mean()
        res_df, filter_extended = common_postprocessing_statistics(
            filter_results, res_df, stat, extend
        )
        filter_extended = _ensure_pass_column_in_extended(
            filter_extended, pass_col, filter_results, "demerit_score"
        )
        return res_df, filter_extended

    raise ValueError(f"Filter {filter_name} not found")


def combine_filter_results_in_memory(output_dir, input_df, pass_mask_by_filter):
    """Combine per-filter pass masks in-memory and persist stage outputs."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    identity_cols = [
        c for c in ["smiles", "model_name", "mol_idx"] if c in input_df.columns
    ]
    combined = input_df[identity_cols].copy()

    pass_columns = []
    for filter_name, pass_df in pass_mask_by_filter.items():
        col_name = f"pass_{filter_name}"
        pass_columns.append(col_name)

        if pass_df is None or len(pass_df) == 0:
            combined[col_name] = False
            continue

        local_df = pass_df.copy()
        pass_col = None
        if "pass" in local_df.columns:
            pass_col = "pass"
        elif "pass_filter" in local_df.columns:
            pass_col = "pass_filter"

        if pass_col is None:
            combined[col_name] = False
            continue

        merge_cols = [c for c in identity_cols if c in local_df.columns]
        if not merge_cols:
            combined[col_name] = False
            continue

        cols_to_take = merge_cols + [pass_col]
        local_df = local_df[cols_to_take].drop_duplicates(
            subset=merge_cols, keep="last"
        )
        local_df = local_df.rename(columns={pass_col: col_name})
        combined = combined.merge(local_df, on=merge_cols, how="left")
        combined[col_name] = combined[col_name].fillna(False).astype(bool)

    if pass_columns:
        combined["_pass_all"] = combined[pass_columns].all(axis=1)
    else:
        combined["_pass_all"] = False

    filtered_df = combined[combined["_pass_all"]][identity_cols].copy()
    filtered_df.to_csv(output_path / "filtered_molecules.csv", index=False)

    failed_cols = identity_cols + pass_columns
    failed_df = combined[~combined["_pass_all"]][failed_cols].copy()
    failed_df.to_csv(output_path / "failed_molecules.csv", index=False)
    return filtered_df, failed_df


def filter_data(config, stage_dir):
    """Filter and combine data from all structural filters.

    Args:
        config: Configuration dictionary
        stage_dir: Stage directory path
    """
    base_folder = Path(process_path(config["folder_to_save"]))
    folder_to_save = base_folder / stage_dir

    paths = list(folder_to_save.glob("*/filtered_molecules.csv"))
    if not paths:
        paths = list(folder_to_save.glob("*filteredMols.csv"))
    paths = [str(p) for p in paths]

    columns_to_drop = ["pass", "any_pass", "name", "pass_any"]
    datas = []
    for path in paths:
        data = pd.read_csv(path)
        for col in columns_to_drop:
            if col in data.columns:
                data.drop(columns=[col], inplace=True)
        datas.append(data)

    if len(datas) > 0:
        filtered_data = datas[0].copy()

        for df in datas[1:]:
            merge_cols = ["smiles", "model_name"]
            existing_cols = set(filtered_data.columns) - set(merge_cols)
            new_cols = [
                col
                for col in df.columns
                if col not in existing_cols and col not in merge_cols
            ]

            if new_cols:
                cols_to_merge = merge_cols + new_cols
                filtered_data = filtered_data.merge(
                    df[cols_to_merge], on=merge_cols, how="inner"
                )
            else:
                filtered_data = filtered_data.merge(
                    df[merge_cols], on=merge_cols, how="inner"
                )
    else:
        filtered_data = pd.DataFrame(columns=["smiles", "model_name", "mol_idx"])

    if "mol_idx" not in filtered_data.columns:
        filtered_data["mol_idx"] = None

    cols = ["smiles", "model_name", "mol_idx"]
    out_df = filtered_data[cols].copy()
    out_df.to_csv(folder_to_save / "filtered_molecules.csv", index=False)

    is_post_descriptors = (
        "03_structural_filters_post" in stage_dir or stage_dir == "StructFilters"
    )
    if is_post_descriptors:
        descriptors_path = base_folder / "Descriptors" / "passDescriptorsSMILES.csv"
        if descriptors_path.exists():
            input_path = str(descriptors_path)
        else:
            sampled_path = base_folder / "sampled_molecules.csv"
            if sampled_path.exists():
                input_path = str(sampled_path)
            else:
                try:
                    from hedgehog.stages.structFilters.main import _get_input_path

                    input_path = _get_input_path(config, stage_dir, str(base_folder))
                except (ImportError, OSError, KeyError):
                    input_path = None
    else:
        sampled_path = base_folder / "sampled_molecules.csv"
        if sampled_path.exists():
            input_path = str(sampled_path)
        else:
            try:
                from hedgehog.stages.structFilters.main import _get_input_path

                input_path = _get_input_path(config, stage_dir, str(base_folder))
            except (ImportError, OSError, KeyError):
                input_path = None

    if input_path and Path(input_path).exists():
        try:
            all_input = pd.read_csv(input_path)
            if len(out_df) > 0:
                merge_cols = ["smiles", "model_name"]
                merged = all_input.merge(
                    out_df[merge_cols], on=merge_cols, how="left", indicator=True
                )
                fail_molecules = merged[merged["_merge"] == "left_only"].drop(
                    columns=["_merge"]
                )
            else:
                fail_molecules = all_input.copy()

            if len(fail_molecules) > 0:
                extended_paths = [str(p) for p in folder_to_save.glob("*extended.csv")]
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
                                        all_extended[col] = all_extended[col].fillna(
                                            all_extended[f"{col}_dup"]
                                        )
                                        all_extended = all_extended.drop(
                                            columns=[f"{col}_dup"]
                                        )
                    except (KeyError, ValueError, TypeError):
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
                            all_extended[cols_to_merge], on=merge_cols, how="left"
                        )
                        for col in pass_cols:
                            if col in fail_molecules.columns:
                                fail_molecules[col] = fail_molecules[col].fillna(False)

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
                    folder_to_save / "failed_molecules.csv", index=False
                )
        except (OSError, pd.errors.ParserError, KeyError) as e:
            logger.warning("Could not create failStructFiltersSMILES.csv: %s", e)

    return filtered_data
