from typing import Union
from pathlib import Path

import yaml
import logging


logger = logging.getLogger(__name__)


class ConfigFileManager:

    def __init__(self, config_file: Union[Path, str]):
        self._config_file = config_file

    def read_config(self):
        if not self._config_file.exists():
            raise Exception(f'Config file {self._config_file} does not exist')

        logger.info(f'Reading configuration from {self._config_file}')

        with open(self._config_file) as file:
            config = yaml.load(file, Loader=yaml.FullLoader)
            if config is None:
                config = {}

            if 'database_connection' not in config:
                logger.warning(f'Missing database_connection in config file {self._config_file}')
            if 'database_dir' not in config:
                logger.warning(f'Missing database_dir in config file {self._config_file}')

            return config