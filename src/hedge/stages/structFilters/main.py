import glob
import os

import pandas as pd

from hedge.configs.logger import load_config, logger
from hedge.stages.structFilters.utils import (
    camelcase,
    filter_data,
    filter_function_applier,
    get_basic_stats,
    inject_identity_columns_to_all_csvs,
    plot_calculated_stats,
    plot_filter_failures_analysis,
    plot_restriction_ratios,
    process_one_file,
    process_path,
)


def _order_identity_columns(df):
    """Order dataframe columns with identity columns first."""
    id_cols = ["smiles", "model_name", "mol_idx"]
    existing_id_cols = [c for c in id_cols if c in df.columns]
    ordered_cols = existing_id_cols + [c for c in df.columns if c not in id_cols]
    return df[ordered_cols]


def _get_input_path(config, prefix, folder_to_save):
    """Determine input path based on prefix.
    
    If running a stage independently and previous stage outputs don't exist,
    falls back to original molecules from config.
    """
    if prefix != "beforeDescriptors":
        descriptors_path = folder_to_save + "Descriptors/passDescriptorsSMILES.csv"
        if os.path.exists(descriptors_path):
            return descriptors_path
        sampled_path = folder_to_save + "sampledMols.csv"
        if os.path.exists(sampled_path):
            logger.info("Descriptors output not found, using sampledMols.csv")
            return sampled_path
        logger.info("Descriptors output not found, using molecules from config")
        return config["generated_mols_path"]

    matched = glob.glob(config["generated_mols_path"])
    if len(matched) > 1:
        return folder_to_save + "sampledMols.csv"

    single_path = matched[0]
    try:
        df_check = pd.read_csv(single_path)
        lower_cols = {c.lower(): c for c in df_check.columns}
        candidate = lower_cols.get("model_name") or lower_cols.get("name")
        if candidate and df_check[candidate].nunique(dropna=True) > 1:
            return folder_to_save + "sampledMols.csv"
    except Exception:
        pass

    return config["generated_mols_path"]


def main(config, prefix):
    """Main entry point for structural filters stage."""
    sample_size = config["sample_size"]
    folder_to_save = process_path(config["folder_to_save"])

    subfolder = f"{prefix}_StructFilters" if prefix == "beforeDescriptors" else "StructFilters"
    os.makedirs(folder_to_save + f"{subfolder}/", exist_ok=True)
    input_path = _get_input_path(config, prefix, folder_to_save)

    try:
        input_df = pd.read_csv(input_path)
        model_names = sorted(input_df["model_name"].dropna().unique().tolist())
        model_name = model_names[0] if len(model_names) == 1 else model_names
    except Exception as e:
        logger.error(f"Could not load input data: {e}")
        raise

    config_structFilters = load_config(config["config_structFilters"])
    filters_to_calculate = {k.replace("calculate_", ""): v
                            for k, v in config_structFilters.items()
                            if "calculate_" in k and v
                          }

    for filter_name in filters_to_calculate:
        apply_func = filter_function_applier(filter_name)
        filter_results = process_one_file(config, input_path, apply_func, sample_size)

        if filter_results is None:
            logger.warning("No molecules to process")
            continue

        final_res, final_extended = get_basic_stats(config_structFilters, filter_results, model_name, filter_name=filter_name)
        path_to_save = folder_to_save + f"{subfolder}/{camelcase(filter_name)}"

        final_res = _order_identity_columns(final_res)
        final_res.to_csv(f"{path_to_save}_metrics.csv", index=False)

        final_extended = _order_identity_columns(final_extended)
        final_extended.to_csv(f"{path_to_save}_extended.csv", index=False)

        if "pass" not in final_extended.columns:
            if "pass_filter" in final_extended.columns:
                final_extended["pass"] = final_extended["pass_filter"]
            else:
                logger.warning(f"'pass' column not found in {filter_name} extended results. Assuming all molecules pass.")
                final_extended["pass"] = True
        filtered_mols = final_extended[final_extended["pass"]].copy()

        filtered_mols = _order_identity_columns(filtered_mols)
        filtered_mols.to_csv(f"{path_to_save}_filteredMols.csv", index=False)

    plot_calculated_stats(config, prefix)
    plot_restriction_ratios(config, prefix)

    if prefix != "beforeDescriptors":
        plot_filter_failures_analysis(config, prefix)

    is_single_stage = config.get("_run_single_stage_override") == "struct_filters"
    if config_structFilters.get("filter_data", False) or is_single_stage:
        filter_data(config, prefix)

    inject_identity_columns_to_all_csvs(config, prefix)
