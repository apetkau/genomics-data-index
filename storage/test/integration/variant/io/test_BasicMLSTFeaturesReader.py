from typing import List
import pytest

from storage.variant.io.BasicMLSTFeaturesReader import BasicMLSTFeaturesReader
from storage.test.integration.variant import basic_mlst_file


@pytest.fixture
def mlst_reader() -> BasicMLSTFeaturesReader:
    return BasicMLSTFeaturesReader(mlst_file=basic_mlst_file)


def expand_list_by(list_in: List[str], number: int) -> List[str]:
    new_list = []
    for value in list_in:
        new_list.extend([value] * number)

    return new_list


def test_get_features_table(mlst_reader):
    num_samples = 6
    num_loci = 7

    mlst_df = mlst_reader.get_features_table()

    assert ['File', 'Sample', 'Scheme', 'Locus', 'Allele', 'Sequence Type'] == list(mlst_df.columns)

    assert num_samples * num_loci == len(mlst_df)
    assert expand_list_by(['2014C-3598.fasta', '2014C-3599.fasta', '2014D-0067.fasta', '2014D-0068.fasta',
            'CFSAN002349.fasta', 'CFSAN023463.fasta'], num_loci) == list(mlst_df['File'].tolist())
    assert expand_list_by(['2014C-3598', '2014C-3599', '2014D-0067', '2014D-0068',
            'CFSAN002349', 'CFSAN023463'], num_loci) == list(mlst_df['Sample'].tolist())

    assert ['abcZ', 'bglA', 'cat', 'dapE', 'dat', 'ldh', 'lhkA'] == list(mlst_df.loc[mlst_df['Sample'] == 'CFSAN002349',
                                                                               'Locus'].tolist())
    assert ['1', '51', '11', '13', '2', '5', '5'] == list(mlst_df.loc[mlst_df['Sample'] == 'CFSAN002349',
                                                                               'Allele'].tolist())
    assert ['adk', 'fumC', 'gyrB', 'icd', 'mdh', 'purA', 'recA'] == list(mlst_df.loc[mlst_df['Sample'] == '2014C-3598',
                                                                               'Locus'].tolist())
    assert ['100', '23', '68', '45', '1', '35', '7'] == list(mlst_df.loc[mlst_df['Sample'] == '2014C-3598',
                                                                               'Allele'].tolist())


def test_samples_list(mlst_reader):
    assert {'2014C-3598', '2014C-3599', '2014D-0067', '2014D-0068',
            'CFSAN002349', 'CFSAN023463'} == set(mlst_reader.samples_list())
