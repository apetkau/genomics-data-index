from os import path, listdir
from pathlib import Path
from typing import List

import pytest

from genomics_data_index.test.integration import data_dir
from genomics_data_index.test.integration import data_dir_empty


@pytest.fixture
def sample_dirs() -> List[Path]:
    return [data_dir / d for d in listdir(data_dir) if path.isdir(data_dir / d)]


@pytest.fixture
def sample_dirs_empty() -> List[Path]:
    return [data_dir_empty / d for d in listdir(data_dir_empty) if path.isdir(data_dir_empty / d)]
