import yaml
from os import path

import logging

logger = logging.getLogger(__file__)


def yaml_config_provider(file_path: str, cmd_name: str):
    if not path.exists(file_path):
        logger.warning(f'No configuration file found in [{file_path}]')
        return {}
    else:
        with open(file_path, 'r') as config_data:
            logger.debug(f'Reading config file {file_path}')
            data = yaml.safe_load(config_data)
            logger.debug(data)
            return data