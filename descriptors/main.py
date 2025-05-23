import os
import pandas as pd
from descriptors.utils import *
from logger_config import logger

from syba import syba


def main(generated_mols_path, path_to_save):
    syba_model = syba.SybaClassifier()
    try:
        syba_model.fitDefaultScore()
    except Exception as e:
        logger.error(f"Failed to load SYBA model: {str(e)}")
        syba_model = None

    path_to_save = path_to_save + '/descriptors/'
    if not os.path.exists(path_to_save):
        os.makedirs(path_to_save)
    generated_mols_df = pd.read_csv(generated_mols_path, names=['smiles'])

    df_inhibitors_clean_seen, df_inhibitors_clean_unseen, df_inhibitors_total = load_inhibitors(path_to_save)
    dict_generated = loading_generated_mols(generated_mols_path)
    dict_generated = dict(sorted(dict_generated.items(), key=lambda x: int(''.join(filter(str.isdigit, x[0])))))

    check_intersection(dict_generated, df_inhibitors_clean_seen, df_inhibitors_clean_unseen, df_inhibitors_total)
    df_to_compare_with = df_inhibitors_total

    tanimoto_similarity_claculation(dict_generated, df_to_compare_with, path_to_save)
    descriptors_dict = compute_descriptors(dict_generated)
    dict_generated = compute_mce18_score(dict_generated)
    if syba_model is not None:
        dict_generated = compute_syba_score(dict_generated, syba_model)
        target_descriptors_df = computing_target_descriptors(df_to_compare_with, syba_model)
    else:
        target_descriptors_df = computing_target_descriptors(df_to_compare_with, None)
    dict_generated["Target"] = target_descriptors_df

    metrics_dict = collect_metrics_dict(dict_generated, descriptors_dict, target_descriptors_df)

    save_plots(metrics_dict, dict_generated, path_to_save, df_to_compare_with) 

    calculate_metrics(generated_mols_df, syba_model, path_to_save)
    calculate_cycle_metrics(generated_mols_df, path_to_save)

