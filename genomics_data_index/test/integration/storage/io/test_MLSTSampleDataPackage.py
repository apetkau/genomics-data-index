from typing import cast

import pandas as pd
import pytest

from genomics_data_index.storage.io.mlst.MLSTFeaturesReader import MLSTFeaturesReader
from genomics_data_index.storage.io.mlst.MLSTSampleData import MLSTSampleData
from genomics_data_index.storage.io.mlst.MLSTSampleDataPackage import MLSTSampleDataPackage
from genomics_data_index.storage.io.mlst.MLSTTSeemannFeaturesReader import MLSTTSeemannFeaturesReader
from genomics_data_index.test.integration import basic_mlst_file


@pytest.fixture
def mlst_reader() -> MLSTTSeemannFeaturesReader:
    return MLSTTSeemannFeaturesReader(mlst_file=basic_mlst_file)


@pytest.fixture
def sample_data_package(mlst_reader: MLSTFeaturesReader) -> MLSTSampleDataPackage:
    return MLSTSampleDataPackage(mlst_reader=mlst_reader)


def test_get_sample_names(sample_data_package: MLSTSampleDataPackage):
    assert {'2014C-3598', '2014C-3599', '2014D-0067', '2014D-0068',
            'CFSAN002349', 'CFSAN023463'} == sample_data_package.sample_names()


def test_iter_samples(sample_data_package: MLSTSampleDataPackage):
    sample_data_map = {}
    for sample_data in sample_data_package.iter_sample_data():
        sample_data_map[sample_data.sample_name] = cast(MLSTSampleData, sample_data)

    assert 6 == len(sample_data_map)
    assert {'2014C-3598', '2014C-3599', '2014D-0067', '2014D-0068',
            'CFSAN002349', 'CFSAN023463'} == set(sample_data_map.keys())

    allele_data = pd.Series({'adk': '100', 'fumC': '23', 'gyrB': '68', 'icd': '45',
                             'mdh': '1', 'purA': '35', 'recA': '7'}, name='Allele')
    allele_data.index.name = 'Locus'
    assert allele_data.equals(sample_data_map['2014C-3598'].get_alleles())
