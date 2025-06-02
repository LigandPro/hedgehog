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
    parser.add_argument('--generated_mols_path', type=str, required=True)
    parser.add_argument('--config',              type=str,                 default='config.yml')

    parser.add_argument('--path_to_save',        type=str, required=False, default='./descriptors/results/')
    parser.add_argument('--smiles_col_name',     type=str, required=False, default='smiles')
    return parser.parse_args()


def main():
    args = parse_args()
    config = load_config(args.config)
    
    logger.info(f'Loading generated mols from {args.generated_mols_path}...')
    try:
        data = pd.read_csv(args.generated_mols_path)
        one_df = True
        logger.info(f'Loaded {len(data)} generated mols.')
        
    except Exception as e:
        paths = glob.glob(args.generated_mols_path)
        data = [pd.read_csv(path) for path in paths]
        one_df = False
        logger.info(f'Loaded {len(data)} dataframes.')

    logger.info(f'Start calculating metrics...\n')
        
    if one_df: data = calculate_metrics(args.generated_mols_path, args.path_to_save, config)
    else:
        for df in data:
            calculate_metrics(df, args.path_to_save, config)



if __name__ == '__main__':
    main()