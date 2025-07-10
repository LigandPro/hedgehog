import os 

from configs.config_utils import load_config
from logger_config import logger

from descriptors.main import main as descriptors_main
from structFilters.main import main as structural_filters_main
from molScore.main import main as molscore_main
# from retroSynth.main import main as retrosynthesis_main
from docking.main import main as docking_main


def check_left_any_data(config, key_word):
    if 'after_descriptors' in key_word:
        key_word = key_word.replace('after_descriptors', '')
        path = config['folder_to_save'] + f'after_descriptors{key_word}/' + f'leftMolsAfter{key_word}SMILES.csv'
    else:
        path = config['folder_to_save'] + f'{key_word}/' + f'leftMolsAfter{key_word}SMILES.csv'
    if os.path.exists(path):
        return True
    else:
        return False


def run_descriptors(data, config):
    config_descriptors = load_config(config['config_descriptors'])
    if config_descriptors['run']:
        logger.info(f'-----> Calculating Descriptors...')
        descriptors_main(data, config)


def run_structural_filters(config, prefix):
    if prefix == 'beforeDescriptors':
        config_structFilters = load_config(config['config_structFilters'])
        if config_structFilters['run']:
            logger.info(f'-----> Running Structural Filters...')
            structural_filters_main(config, prefix)
    else:
        if check_left_any_data(config, prefix):
            config_structFilters = load_config(config['config_structFilters'])
            if config_structFilters['run']:
                logger.info(f'-----> Running Structural Filters...')
                structural_filters_main(config, prefix)
        else:
            logger.warning(f'No data to process for Structural Filters')    


def run_molscore_evaluation(config):
    if check_left_any_data(config, 'StructFilters'):
        config_molscore = load_config(config['config_molScore'])
        if config_molscore['run']:
            logger.info(f'-----> Running MolScore Evaluation...')
            molscore_main(config)
    else:
        logger.warning(f'No data to process for MolScore Evaluation')

# def run_rethrosynth(config):
#     config_retroSynth = load_config(config['config_retroSynth'])
#     if config_retroSynth['run']:
#         logger.info(f'-----> Running Retrosynthesis...')
#         retrosynthesis_main()


def run_docking(config):
    if check_left_any_data(config, 'StructFilters'):
        config_docking = load_config(config['config_docking'])
        if config_docking['run']:
            logger.info(f'-----> Running docking...')
            docking_main(config)
    else:
        logger.warning(f'No data to process for docking')



def calculate_metrics(data, config):
    run_structural_filters(config, prefix='beforeDescriptors')
    run_descriptors(data, config)
    run_structural_filters(config, prefix='Descriptors')
    run_molscore_evaluation(config)
    # run_rethrosynth(config)  
    # run_docking(config)
    
    
