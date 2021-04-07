import pytest

from storage.variant.io.BasicMLSTFeaturesReader import BasicMLSTFeaturesReader
from storage.test.integration.variant import basic_mlst_file


@pytest.fixture
def mlst_reader() -> BasicMLSTFeaturesReader:
    return BasicMLSTFeaturesReader(mlst_file=basic_mlst_file)


def test_get_features_table(mlst_reader):
    mlst_df = mlst_reader.get_features_table()
    print(mlst_df)
    assert False