import pandas as pd
import pytest

from genomics_data_index.storage.io.mlst.MLSTSampleData import MLSTSampleData


@pytest.fixture
def mlst_sample_data() -> MLSTSampleData:
    allele_data = pd.Series({'l1': '1', 'l2': '2'}, name='Allele', dtype=str)
    return MLSTSampleData(sample_name='SampleA', scheme='lmonocytogenes',
                          mlst_allele_data=allele_data)


def test_single_sample_data(mlst_sample_data: MLSTSampleData):
    assert 'lmonocytogenes' == mlst_sample_data.get_scheme()
    assert 'SampleA' == mlst_sample_data.sample_name
    assert mlst_sample_data.is_preprocessed()
    assert 2 == len(mlst_sample_data.get_alleles())
