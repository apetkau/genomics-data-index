import pandas as pd
import pytest

from genomics_data_index.storage.io.mlst.MLSTTSeemannFeaturesReader import MLSTTSeemannFeaturesReader
from genomics_data_index.test.integration import basic_mlst_file, mlst_file_unknown, expand_list_by


@pytest.fixture
def mlst_reader() -> MLSTTSeemannFeaturesReader:
    return MLSTTSeemannFeaturesReader(mlst_file=basic_mlst_file)


@pytest.fixture
def mlst_reader_unknown() -> MLSTTSeemannFeaturesReader:
    return MLSTTSeemannFeaturesReader(mlst_file=mlst_file_unknown)


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
    assert ['1', '51', '11', '13', '2', '5', '4'] == list(mlst_df.loc[mlst_df['Sample'] == 'CFSAN002349',
                                                                      'Allele'].tolist())
    assert ['adk', 'fumC', 'gyrB', 'icd', 'mdh', 'purA', 'recA'] == list(mlst_df.loc[mlst_df['Sample'] == '2014C-3598',
                                                                                     'Locus'].tolist())
    assert ['100', '23', '68', '45', '1', '35', '7'] == list(mlst_df.loc[mlst_df['Sample'] == '2014C-3598',
                                                                         'Allele'].tolist())


def test_get_features_table_with_unknown(mlst_reader_unknown):
    num_samples = 4
    num_loci = 7

    mlst_df = mlst_reader_unknown.get_features_table()

    assert ['File', 'Sample', 'Scheme', 'Locus', 'Allele', 'Sequence Type'] == list(mlst_df.columns)

    assert num_samples * num_loci == len(mlst_df)

    assert ['abcZ', 'bglA', 'cat', 'dapE', 'dat', 'ldh', 'lhkA'] == list(mlst_df.loc[mlst_df['Sample'] == 'CFSAN002349',
                                                                                     'Locus'].tolist())
    assert ['1', '?', '11', '13', '2', '5', '5'] == list(mlst_df.loc[mlst_df['Sample'] == 'CFSAN002349',
                                                                     'Allele'].tolist())
    assert ['abcZ', 'bglA', 'cat', 'dapE', 'dat', 'ldh', 'lhkA'] == list(mlst_df.loc[mlst_df['Sample'] == 'CFSAN023463',
                                                                                     'Locus'].tolist())
    assert ['?', '?', '11', '13', '2', '5', '5'] == list(mlst_df.loc[mlst_df['Sample'] == 'CFSAN023463',
                                                                     'Allele'].tolist())
    assert ['adk', 'fumC', 'gyrB', 'icd', 'mdh', 'purA', 'recA'] == list(mlst_df.loc[mlst_df['Sample'] == '2014C-3598',
                                                                                     'Locus'].tolist())
    assert ['100', '?', '?', '45', '1', '35', '7'] == list(mlst_df.loc[mlst_df['Sample'] == '2014C-3598',
                                                                       'Allele'].tolist())


def test_iter_data(mlst_reader: MLSTTSeemannFeaturesReader):
    sample_data_map = {}
    for sample_name, scheme, allele_data in mlst_reader.iter_sample_data():
        if sample_name in sample_data_map:
            raise Exception(f'Duplicate sample {sample_name} in sample_data_map')
        sample_data_map[sample_name] = {'scheme': scheme, 'allele_data': allele_data}

    assert 6 == len(sample_data_map)

    sample_schemes = {(n, sample_data_map[n]['scheme']) for n in sample_data_map}
    assert {('2014C-3598', 'ecoli'), ('2014C-3599', 'ecoli'),
            ('2014D-0067', 'campylobacter'), ('2014D-0068', 'campylobacter'),
            ('CFSAN002349', 'lmonocytogenes'), ('CFSAN023463', 'lmonocytogenes')} == sample_schemes

    allele_data = pd.Series({'adk': '100', 'fumC': '23', 'gyrB': '68', 'icd': '45',
                             'mdh': '1', 'purA': '35', 'recA': '7'}, name='Allele')
    allele_data.index.name = 'Locus'
    assert allele_data.equals(sample_data_map['2014C-3598']['allele_data'])

    allele_data = pd.Series({'abcZ': '1', 'bglA': '51', 'cat': '11', 'dapE': '13',
                             'dat': '2', 'ldh': '5', 'lhkA': '5'}, name='Allele')
    allele_data.index.name = 'Locus'
    assert allele_data.equals(sample_data_map['CFSAN023463']['allele_data'])


def test_samples_list(mlst_reader):
    assert {'2014C-3598', '2014C-3599', '2014D-0067', '2014D-0068',
            'CFSAN002349', 'CFSAN023463'} == set(mlst_reader.samples_list())
