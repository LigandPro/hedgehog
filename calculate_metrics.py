import logging

from metrics_utils import *
from logger_config import logger

from descriptors_main import main as descriptors_main
from structural_filters_main import main as structural_filters_main
from retrosynthesis_main import main as retrosynthesis_main
from smina_score_main import main as smina_score_main


def run_descriptors(generated_mols_path, path_to_save):
    logger.info(f'Calculating Descriptors...')
    descriptors_main(generated_mols_path, path_to_save)


def run_structural_filters():
    logger.info(f'Running Structural Filters...')
    structural_filters_main()

def run_rethrosynth():
    logger.info(f'Running Retrosynthesis...')
    retrosynthesis_main()

def run_smina_score():
    logger.info(f'Running SMINA Score...')
    smina_score_main()


def calculate_metrics(generated_mols_path, path_to_save, config):
    run_descriptors_flag = config['descriptors']['run']
    run_structural_filters_flag = config['structure_filters']['run']
    run_rethrosynth_flag = config['retrosynthesis']['run']
    run_smina_score_flag = config['smina_score']['run']
    
    if run_descriptors_flag:
        run_descriptors(generated_mols_path, path_to_save)

    if run_structural_filters_flag:
        run_structural_filters()

    if run_rethrosynth_flag:
        run_rethrosynth()
        
    if run_smina_score_flag:
        run_smina_score()
    
    
