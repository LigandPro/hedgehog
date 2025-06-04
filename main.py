import os
import argparse
import glob 
import yaml
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns 

from calculate_metrics import calculate_metrics
from logger_config import logger
from config_utils import load_config



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config',              type=str,                 default='config.yml')
    return parser.parse_args()


def main():
    args = parse_args()
    config = load_config(args.config)
    generated_mols_path = config['main']['generated_mols_path']
    folder_to_save = config['main']['folder_to_save']
    smiles_col_name = config['main']['smiles_col_name']
    save_mols_list = config['main']['save_mols_list']
    number_of_mols_to_save = config['main']['number_of_mols_to_save'] if save_mols_list else None

    logger.info(f'Loading generated mols from {generated_mols_path}...')
    try:
        data = pd.read_csv(generated_mols_path, names=[smiles_col_name])
        one_df = True

    except Exception as e:
        paths = glob.glob(generated_mols_path)
        data = [pd.read_csv(path) for path in paths]
        one_df = False

    if save_mols_list and one_df:
        data = data.sample(number_of_mols_to_save)
        os.makedirs(folder_to_save, exist_ok=True)
        data.to_csv(folder_to_save + f'{number_of_mols_to_save}mols.csv', index=False)
        logger.info(f'Sampled {len(data)} generated mols.')
    if not one_df:
        logger.info(f'Loaded {len(data)} dataframes.')

    logger.info(f'Start calculating metrics...\n')
        
    if one_df: data = calculate_metrics(config)
    # else:
    #     for df in data:
    #         calculate_metrics(df, config)



if __name__ == '__main__':
    main()