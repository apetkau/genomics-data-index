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

        self._config_file_manager = ConfigFile.read_config(self._config_file)
        cfm = self._config_file_manager
        fs_dir = cfm.database_dir

        if not fs_dir.exists():
            raise Exception(f'Data directory [{fs_dir}] does not exist. '
                            f'Please verify that project directory [{root_dir}] is valid.')

        if cfm.database_connection is None and cfm.sqlite_database is None:
            raise Exception(f'No valid configured database in config file {self._config_file}')

    def create_connection(self) -> DataIndexConnection:
        cfm = self._config_file_manager
        database_connection = None
        if cfm.database_connection is None:
            if cfm.sqlite_database is None:
                raise Exception(f'No valid configured database in config file {self._config_file}')
            else:
                # Provide support for defining sqlite file relative to project directory or as an absolute string
                sqlite_path = None
                if cfm.sqlite_database.is_absolute():
                    sqlite_path = cfm.sqlite_database
                else:
                    sqlite_path = self._root_dir / cfm.sqlite_database

            database_connection = f'sqlite:///{sqlite_path}'
        elif cfm.database_connection is None:
            raise Exception(f'No valid configured database in config file {self._config_file}')
        else:
            database_connection = cfm.database_connection

        return DataIndexConnection.connect(database_connection=database_connection,
                                           database_dir=self._config_file_manager.database_dir)

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
        config.database_dir = data_dir

        config.write(project_dir / cls.CONFIG_FILE)

        return Project(project_dir)
