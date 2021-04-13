from typing import List

import pytest

from storage.test.integration.variant import sistr_mlst_file
from storage.variant.io.mlst.MLSTSistrReader import MLSTSistrReader



@pytest.fixture
def mlst_reader() -> MLSTSistrReader:
    return MLSTSistrReader(mlst_file=sistr_mlst_file)


def expand_list_by(list_in: List[str], number: int) -> List[str]:
    new_list = []
    for value in list_in:
        new_list.extend([value] * number)

    return new_list


def test_get_features_table(mlst_reader):
    num_samples = 2
    num_loci = 330

    mlst_df = mlst_reader.get_features_table()

    assert ['File', 'Sample', 'Scheme', 'Locus', 'Allele'] == list(mlst_df.columns)

    assert num_samples * num_loci == len(mlst_df)
    assert expand_list_by(['GCF_000006945.2_ASM694v2_genomic', 'GCF_000007545.1_ASM754v1_genomic'],
                          num_loci) == list(mlst_df['File'].tolist())
    assert expand_list_by(['GCF_000006945', 'GCF_000007545'],
                          num_loci) == list(mlst_df['Sample'].tolist())

    assert {'sistr_330'} == set(mlst_df['Scheme'].tolist())

    assert ['3371366009'] == list(mlst_df.loc[(mlst_df['Sample'] == 'GCF_000006945') &
                                       (mlst_df['Locus'] == 'NZ_AOXE01000059.1_60'), 'Allele'].tolist())

    assert ['3848099890'] == list(mlst_df.loc[(mlst_df['Sample'] == 'GCF_000007545') &
                                       (mlst_df['Locus'] == 'NZ_AOXE01000019.1_14'), 'Allele'].tolist())
