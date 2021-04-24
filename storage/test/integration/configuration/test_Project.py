from pathlib import Path
from tempfile import TemporaryDirectory

import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)

from storage.test.integration import test_project_dir
from storage.configuration.Project import Project
from storage.configuration.ConfigFile import ConfigFile


def test_create_new_project():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)

        Project.create_new_project(tmp_dir / 'project')

        assert (tmp_dir / 'project').exists()
        assert (tmp_dir / 'project' / 'data').exists()
        assert (tmp_dir / 'project' / 'config.yaml').exists()

        config = ConfigFile.read_config(tmp_dir / 'project' / 'config.yaml')
        assert config.database_connection is None
        assert config.database_dir == Path('data')
        assert config.sqlite_database == Path('db.sqlite')


def test_from_existing_project():
    project = Project(test_project_dir)

    assert f'sqlite:///{test_project_dir}/db.sqlite' == project.get_database_connection_str()
    assert test_project_dir / 'data' == project.get_database_dir()

    connection = project.create_connection()
    assert connection.filesystem_storage.variation_dir.parent == test_project_dir / 'data'
