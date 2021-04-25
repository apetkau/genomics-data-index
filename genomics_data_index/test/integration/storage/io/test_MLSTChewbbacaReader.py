import pytest

from genomics_data_index.storage.io.mlst.MLSTChewbbacaReader import MLSTChewbbacaReader
from genomics_data_index.test.integration import chewbbaca_mlst_file, expand_list_by


@pytest.fixture
def mlst_reader() -> MLSTChewbbacaReader:
    return MLSTChewbbacaReader(mlst_file=chewbbaca_mlst_file, scheme='test-scheme')


def test_get_features_table(mlst_reader):
    num_samples = 2
    num_loci = 11

    mlst_df = mlst_reader.get_features_table()

    assert ['File', 'Sample', 'Scheme', 'Locus', 'Allele'] == list(mlst_df.columns)
    assert {'test-scheme'} == set(mlst_df['Scheme'].tolist())

    assert num_samples * num_loci == len(mlst_df)
    assert expand_list_by(['GCF_000006945.2_ASM694v2_genomic.fna', 'GCF_000007545.1_ASM754v1_genomic.fna'],
                          num_loci) == list(mlst_df['File'].tolist())
    assert expand_list_by(['GCF_000006945', 'GCF_000007545'],
                          num_loci) == list(mlst_df['Sample'].tolist())

    assert ['1'] == list(mlst_df.loc[(mlst_df['Sample'] == 'GCF_000006945') &
                                     (mlst_df['Locus'] == 'GCF-000006945-protein10.fasta'), 'Allele'].tolist())

    assert ['2'] == list(mlst_df.loc[(mlst_df['Sample'] == 'GCF_000007545') &
                                     (mlst_df['Locus'] == 'GCF-000006945-protein1004.fasta'), 'Allele'].tolist())


def test_get_features_table_unknown(mlst_reader):
    mlst_df = mlst_reader.get_features_table()

    assert ['File', 'Sample', 'Scheme', 'Locus', 'Allele'] == list(mlst_df.columns)
    assert {'test-scheme'} == set(mlst_df['Scheme'].tolist())

    assert ['?'] == list(mlst_df.loc[(mlst_df['Sample'] == 'GCF_000007545') &
                                     (mlst_df['Locus'] == 'GCF-000006945-protein1.fasta'), 'Allele'].tolist())

    assert ['?'] == list(mlst_df.loc[(mlst_df['Sample'] == 'GCF_000007545') &
                                     (mlst_df['Locus'] == 'GCF-000006945-protein1008.fasta'), 'Allele'].tolist())
