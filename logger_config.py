import logging
import yaml

def setup_logger(config_path='config.yml'):
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
            log_config = config.get('logging', {})
    except:
        log_config = {'level': 'INFO', 'file_output': False}

    logger = logging.getLogger('metrics_comparison')
    
    level = getattr(logging, log_config.get('level', 'INFO'))
    logger.setLevel(level)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)

    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    console_handler.setFormatter(formatter)

    if log_config.get('file_output', False):
        file_handler = logging.FileHandler('metrics_comparison.log')
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    logger.addHandler(console_handler)

    return logger

logger = setup_logger() 