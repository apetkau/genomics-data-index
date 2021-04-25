import shutil
import warnings
from pathlib import Path
from tempfile import TemporaryDirectory

import pytest

warnings.filterwarnings("ignore", category=DeprecationWarning)

from genomics_data_index.test.integration import test_project_dir
from genomics_data_index.configuration.Project import Project
from genomics_data_index.configuration.ConfigFile import ConfigFile


def test_create_new_project():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)

        Project.initialize_project(tmp_dir / 'project')

        assert (tmp_dir / 'project').exists()
        assert (tmp_dir / 'project' / '.gdi-data').exists()
        assert (tmp_dir / 'project' / 'gdi-config.yaml').exists()

        config = ConfigFile.read_config(tmp_dir / 'project' / 'gdi-config.yaml')
        assert config.database_connection is None
        assert config.database_dir == Path('.gdi-data')
        assert config.sqlite_database == Path('.gdi-data') / 'gdi-db.sqlite'


def test_from_existing_project():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        project_dir = tmp_dir / 'project'
        shutil.copytree(test_project_dir, project_dir)

        project = Project(project_dir)

        assert f'sqlite:///{project_dir}/.gdi-data/gdi-db.sqlite' == project.get_database_connection_str()
        assert project_dir / '.gdi-data' == project.get_database_dir()

        connection = project.create_connection()
        assert connection.filesystem_storage.variation_dir.parent == project_dir / '.gdi-data'


def test_initialize_from_existing_project():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)
        project_dir = tmp_dir / 'project'
        shutil.copytree(test_project_dir, project_dir)

        with pytest.raises(Exception) as execinfo:
            Project.initialize_project(project_dir)
        assert 'has already been initialized' in str(execinfo.value)
