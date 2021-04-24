from __future__ import annotations
import logging
from pathlib import Path
import shutil
import os

from storage.configuration.connector.DataIndexConnection import DataIndexConnection
from storage.configuration.ConfigFile import ConfigFile

logger = logging.getLogger(__name__)


class Project:
    DATA_DIR = 'data'
    CONFIG_FILE = 'config.yaml'

    def __init__(self, root_dir: Path):
        self._root_dir = root_dir
        self._config_file = root_dir / self.CONFIG_FILE

        if not self._config_file.exists():
            raise Exception(f'Config file [{self._config_file}] does not exist. '
                            f'Please verify that project directory [{root_dir}] is valid.')

        self._config = ConfigFile.read_config(self._config_file)
        self._database_connection = self.get_database_connection_str()
        self._database_dir = self.get_database_dir()

    def get_database_connection_str(self) -> str:
        if self._config.database_connection is None:
            if self._config.sqlite_database is None:
                raise Exception(f'No valid configured database in config file {self._config_file}')
            else:
                # Provide support for defining sqlite file relative to project directory or as an absolute string
                if self._config.sqlite_database.is_absolute():
                    sqlite_path = self._config.sqlite_database
                else:
                    sqlite_path = self._root_dir / self._config.sqlite_database

            database_connection = f'sqlite:///{sqlite_path}'
        elif self._config.database_connection is None:
            raise Exception(f'No valid configured database in config file {self._config_file}')
        else:
            database_connection = self._config.database_connection

        return database_connection

    def get_database_dir(self) -> Path:
        if self._config.database_dir is None:
            raise Exception(f'database_dir is not configured properly in config file {self._config_file}')

        if self._config.database_dir.is_absolute():
            database_dir = self._config.database_dir
        else:
            database_dir = self._root_dir / self._config.database_dir

        if not database_dir.exists():
            raise Exception(f'database_dir=[{database_dir}] does not exist. '
                            f'Please verify that project directory [{self._root_dir}] is valid.')
        else:
            return database_dir

    def create_connection(self) -> DataIndexConnection:
        return DataIndexConnection.connect(database_connection=self._database_connection,
                                           database_dir=self._database_dir)

    @classmethod
    def create_new_project(cls, project_dir: Path, force_new: bool = False) -> Project:
        if project_dir.exists():
            if force_new:
                logger.warning(f'Project directory [{project_dir}] exists but force_new={force_new} so overwriting.')
                shutil.rmtree(project_dir)
            else:
                raise Exception(f'Project directory [{project_dir}] already exists')

        data_dir = project_dir / cls.DATA_DIR

        os.mkdir(project_dir)
        os.mkdir(data_dir)

        config = ConfigFile()
        config.sqlite_database = 'db.sqlite'
        config.database_dir = cls.DATA_DIR

        config.write(project_dir / cls.CONFIG_FILE)

        return Project(project_dir)
