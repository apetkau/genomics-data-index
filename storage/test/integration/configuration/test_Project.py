from pathlib import Path
from tempfile import TemporaryDirectory

import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)

from storage.configuration.Project import Project


def test_create_new_project():
    with TemporaryDirectory() as tmp_dir_str:
        tmp_dir = Path(tmp_dir_str)

        project = Project.create_new_project(tmp_dir / 'project')

        assert (tmp_dir / 'project').exists()
        assert (tmp_dir / 'project' / 'data').exists()
        assert (tmp_dir / 'project' / 'config.yaml').exists()
