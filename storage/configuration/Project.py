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
        fs_dir = root_dir / self.DATA_DIR
        config_file = root_dir / self.CONFIG_FILE

        if not fs_dir.exists():
            raise Exception(f'Data directory [{fs_dir}] does not exist. '
                            f'Please verify that project directory [{root_dir}] is valid.')

        if not config_file.exists():
            raise Exception(f'Config file [{config_file}] does not exist. '
                            f'Please verify that project directory [{root_dir}] is valid.')

        self._config_file_manager = ConfigFile.read_config(config_file)

    def create_connection(self) -> DataIndexConnection:
        return DataIndexConnection.connect(database_connection=self._config_file_manager.database_connection,
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
        config.database_connection = 'sqlite:///database.sqlite'
        config.database_dir = data_dir

        config.write(project_dir / cls.CONFIG_FILE)

        return Project(project_dir)
