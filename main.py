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


def load_multi_comparison_data(paths, sample_size):
    dataframes = []
    model_names = []

    for path in paths: 
        df = pd.read_csv(path)
        smiles_col = next((col for col in df.columns if col.lower() == 'smiles'), None)
        if smiles_col:
            df = df[[smiles_col]].drop_duplicates().reset_index(drop=True)
            df.columns = ['smiles']
        else:
            df = pd.read_csv(path, header=None)
            df.columns = ['smiles'] + [f'col_{i}' for i in range(1, len(df.columns))]
            df = df[['smiles']].drop_duplicates().reset_index(drop=True)
            
        model_name = get_model_name(path=path)
        if len(df) < sample_size:
            logger.warning(f'Sample size is not equal to {sample_size}. Chosen {len(df)} for {model_name} mols.')
        else:
            df = df.sample(sample_size)

        df['model_name'] = model_name
        dataframes.append(df)

    merged = pd.concat(dataframes, axis=0)
    return merged
            

def main():
    args = parse_args()
    config = load_config(args.config)

    mode = config['mode']
    generated_mols_path = config['generated_mols_path']
    folder_to_save = config['folder_to_save']
    save_sampled_mols = config['save_sampled_mols']
    sample_size = config['sample_size'] if save_sampled_mols else None

    logger.info(f'Loading generated mols from {generated_mols_path}...')
    if mode == 'single_comparison':
        data = pd.read_csv(generated_mols_path)
        
        smiles_col = next((col for col in data.columns if col.lower() == 'smiles'), None)
        if smiles_col:
            data = data[[smiles_col]]
            data.columns = ['smiles']
        else:
            data = pd.read_csv(generated_mols_path, header=None)
            data.columns = ['smiles'] + [f'col_{i}' for i in range(1, len(data.columns))]
            data = data[['smiles']]

        data = data.drop_duplicates(subset='smiles').reset_index(drop=True)
        if len(data) < sample_size:
            logger.warning(f'Sample size is not equal to {sample_size}. Chosen {len(data)} mols.')
        else:
            data = data.sample(sample_size)


    elif mode == 'multi_comparison':
        generated_mols_path = glob.glob(generated_mols_path)
        logger.info(f'Loading multi comparison analysis for {len(generated_mols_path)} files...')
        data = load_multi_comparison_data(generated_mols_path, sample_size)

    else:
        raise ValueError(f'Invalid mode: {mode}')


    if save_sampled_mols:
        os.makedirs(folder_to_save, exist_ok=True)
        data.to_csv(folder_to_save + f'sampledMols.csv', index=False)
        logger.info(f'Sampled {len(data)} generated mols.')
    

    logger.info(f'Start calculating metrics...')
    calculate_metrics(data, config, mode)


if __name__ == '__main__':
    main()