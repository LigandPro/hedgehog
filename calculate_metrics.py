from configs.config_utils import load_config
from logger_config import logger

from descriptors.main import main as descriptors_main
# from structFilters.main import main as structural_filters_main
# from retroSynth.main import main as retrosynthesis_main
# from sminaScore.main import main as smina_score_main


def run_descriptors(data, config):
    config_descriptors = load_config(config['config_descriptors'])
    if config_descriptors['run']:
        logger.info(f'Calculating Descriptors...')
        descriptors_main(data, config)


# def run_structural_filters():
#     if config_structFilters['run']:
    #     logger.info(f'Running Structural Filters...')
    #     structural_filters_main()

# def run_rethrosynth():
#     if config_retroSynth['run']:
    #     logger.info(f'Running Retrosynthesis...')
    #     retrosynthesis_main()

# def run_smina_score():
#     if config_dockingScore['run']:
    #     logger.info(f'Running SMINA Score...')
    #     smina_score_main()


def calculate_metrics(data, config):
    run_descriptors(data, config)
    # run_structural_filters(config)
    # run_rethrosynth(config)  
    # run_smina_score(config)
    
    
