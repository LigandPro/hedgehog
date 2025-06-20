import os 

from configs.config_utils import load_config
from logger_config import logger

from descriptors.main import main as descriptors_main
from structFilters.main import main as structural_filters_main
from molScore_metrics.main import main as molscore_main
# from retroSynth.main import main as retrosynthesis_main
# from sminaScore.main import main as smina_score_main


def check_left_any_data(config, key_word):
    path = config['folder_to_save'] + f'{key_word}/' + f'leftMolsAfter{key_word}SMILES.csv'
    if os.path.exists(path):
        return True
    else:
        return False


def run_descriptors(data, config, mode):
    config_descriptors = load_config(config['config_descriptors'])
    if config_descriptors['run']:
        logger.info(f'Calculating Descriptors...')
        descriptors_main(data, config, mode)


def run_structural_filters(config, mode, prefix):
    if prefix == 'before_descriptors':
        config_structFilters = load_config(config['config_structFilters'])
        if config_structFilters['run']:
            logger.info(f'Running Structural Filters...')
            structural_filters_main(config, mode, prefix)
    else:
        if check_left_any_data(config, 'Descriptors'):
            config_structFilters = load_config(config['config_structFilters'])
            if config_structFilters['run']:
                logger.info(f'Running Structural Filters...')
                structural_filters_main(config, mode, prefix)
        else:
            logger.warning(f'No data to process for Structural Filters')    


def run_molscore_evaluation(config, mode):
    if check_left_any_data(config, 'StructFilters'):
        config_molscore = load_config(config['config_molscore'])
        if config_molscore['run']:
            logger.info(f'Running MolScore Evaluation...')
            molscore_main(config, mode)
    else:
        logger.warning(f'No data to process for MolScore Evaluation')

# def run_molscore_evaluation(config):
#     config_molscore = load_config(config['config_molscore'])
#     if config_molscore['run']:
#         logger.info(f'Running MolScore Evaluation...')
#         molscore_main()


# def run_rethrosynth(config):
#     config_retroSynth = load_config(config['config_retroSynth'])
#     if config_retroSynth['run']:
#         logger.info(f'Running Retrosynthesis...')
#         retrosynthesis_main()


# def run_smina_score():
#     if config_dockingScore['run']:
    #     logger.info(f'Running SMINA Score...')
    #     smina_score_main()


def calculate_metrics(data, config, mode):
    run_structural_filters(config, mode, prefix='before_descriptors')
    run_descriptors(data, config, mode)
    run_structural_filters(config, mode, prefix='after_descriptors')
    # run_molscore_evaluation(config, mode)
    # run_molscore_evaluation(config, mode)
    # run_rethrosynth(config)  
    # run_smina_score(config)
    
    
