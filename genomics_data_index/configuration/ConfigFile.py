from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Dict

import yaml

logger = logging.getLogger(__name__)


class ConfigFile:

    def __init__(self, config_dict: Optional[Dict[str, str]] = None):
        if config_dict is None:
            config_dict = {}
        self._config = config_dict

    @property
    def database_connection(self) -> str:
        return self._config['database_connection']

    @database_connection.setter
    def database_connection(self, connection: str) -> None:
        self._config['database_connection'] = connection

    @property
    def database_dir(self) -> Path:
        return Path(self._config['database_dir'])

    @database_dir.setter
    def database_dir(self, dir: Path) -> None:
        self._config['database_dir'] = str(dir)

    @property
    def sqlite_database(self) -> Path:
        return Path(self._config['sqlite_database'])

    @sqlite_database.setter
    def sqlite_database(self, db_path: Path):
        self._config['sqlite_database'] = str(db_path)

    def write(self, file: Path) -> None:
        """
        Writes The configuration to the passed file.
        :param file: The file where configuration should be written.
        :return: None.
        """
        if file.exists():
            raise Exception(f'File [{file}] already exists')

        with open(file, 'w') as fh:
            yaml.dump(self._config, fh)

        logger.debug(f'Wrote configuration to [{file}]')

    @classmethod
    def read_config(self, file: Path) -> ConfigFile:
        if not file.exists():
            raise Exception(f'Config file {file} does not exist')

        logger.debug(f'Reading configuration from {file}')

        with open(file) as file:
            config = yaml.load(file, Loader=yaml.FullLoader)
            if config is None:
                config = {}

            if 'database_connection' not in config and 'sqlite_database' not in config:
                logger.warning(f'Missing database_connection and sqlite_database in config file {file}')

            if 'database_connection' not in config:
                config['database_connection'] = None
            if 'sqlite_database' not in config:
                config['sqlite_database'] = None
            if 'database_dir' not in config:
                logger.warning(f'Missing database_dir in config file {file}')
                config['database_dir'] = None

            return ConfigFile(config_dict=config)
