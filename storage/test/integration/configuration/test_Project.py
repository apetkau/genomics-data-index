import os
import shutil
import warnings
from pathlib import Path
from tempfile import TemporaryDirectory

import pytest

warnings.filterwarnings("ignore", category=DeprecationWarning)

from storage.test.integration import test_project_dir
from storage.configuration.Project import Project
from storage.configuration.ConfigFile import ConfigFile


def test_create_new_project():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)

        Project.initialize_project(tmp_dir / 'project')

        assert (tmp_dir / 'project').exists()
        assert (tmp_dir / 'project' / '.genomics-data').exists()
        assert (tmp_dir / 'project' / 'data-config.yaml').exists()

        config = ConfigFile.read_config(tmp_dir / 'project' / 'data-config.yaml')
        assert config.database_connection is None
        assert config.database_dir == Path('.genomics-data')
        assert config.sqlite_database == Path('.genomics-data') / 'genomics-db.sqlite'


def test_from_existing_project():
    project = Project(test_project_dir)

    assert f'sqlite:///{test_project_dir}/.genomics-data/db.sqlite' == project.get_database_connection_str()
    assert test_project_dir / '.genomics-data' == project.get_database_dir()

    connection = project.create_connection()
    assert connection.filesystem_storage.variation_dir.parent == test_project_dir / '.genomics-data'


def test_initialize_from_existing_project():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        project_dir = tmp_dir / 'project'
        shutil.copytree(test_project_dir, project_dir)
        print(os.listdir(project_dir))

        with pytest.raises(Exception) as execinfo:
            Project.initialize_project(project_dir)
        assert 'has already been initialized' in str(execinfo.value)
