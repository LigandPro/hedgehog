import os
import argparse
import glob 
import yaml
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns 

from descriptors.utils import get_model_name    
from calculate_metrics import calculate_metrics
from logger_config import logger
from configs.config_utils import load_config

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=str, default='./configs/config.yml')
    return parser.parse_args()


def main():
    args = parse_args()
    config = load_config(args.config)

    generated_mols_path = config['generated_mols_path']
    folder_to_save = config['folder_to_save']
    save_sampled_mols = config['save_sampled_mols']
    sample_size = config['sample_size'] if save_sampled_mols else None

    logger.info(f'Loading generated mols from {generated_mols_path}...')
    try:
        data = pd.read_csv(generated_mols_path)
        
        smiles_col = next((col for col in data.columns if col.lower() == 'smiles'), None)
        if smiles_col:
            data = data[[smiles_col]]
            data.columns = ['smiles']
        else:
            data = pd.read_csv(generated_mols_path, header=None)
            data.columns = ['smiles'] + [f'col_{i}' for i in range(1, len(data.columns))]
            data = data[['smiles']]
        one_df = True
        mode = 'single_comparison'
        data = data.drop_duplicates(subset='smiles').reset_index(drop=True)
        if len(data) < sample_size:
            logger.warning(f'Sample size is not equal to {sample_size}. Chosen {len(data)} mols.')
            sample_size = len(data)

    except Exception as e:
        paths = glob.glob(generated_mols_path)
        data = []
        for path in paths:
            df = pd.read_csv(path)
            smiles_col = next((col for col in df.columns if col.lower() == 'smiles'), None)
            if smiles_col:
                df = df[[smiles_col]]
                col_name = get_model_name(path=path)
                df.columns = [col_name]
            else:
                df = pd.read_csv(generated_mols_path, header=None)
                df.columns = ['smiles'] + [f'col_{i}' for i in range(1, len(data.columns))]
                df = df[['smiles']]
                col_name = get_model_name(path=path)
                df.columns = [col_name]
            
            df = df.drop_duplicates(subset=df.columns[0]).reset_index(drop=True)
            data.append(df)
            
        one_df = False
        mode = 'multi_comparison'
        logger.info(f'Loaded {len(data)} dataframes.')


    if save_sampled_mols and one_df:
        data = data.sample(sample_size)
        os.makedirs(folder_to_save, exist_ok=True)
        data.to_csv(folder_to_save + f'{sample_size}mols.csv', index=False)
        logger.info(f'Sampled {len(data)} generated mols.')

    logger.info(f'Start calculating metrics...\n')
        
    calculate_metrics(data, config, mode=mode)


if __name__ == '__main__':
    main()