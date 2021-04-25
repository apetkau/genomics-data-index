import tempfile
from pathlib import Path

from sourmash import load_file_as_signatures

from genomics_data_index.storage.index.KmerIndexer import KmerIndexManager
from genomics_data_index.storage.index.KmerIndexer import KmerIndexerSourmash
from genomics_data_index.test.integration import data_dir
from genomics_data_index.test.integration.storage.index import get_values_from_signatures


def test_index_multiple_files():
    sampleA = data_dir / 'SampleA' / 'snps.aligned.nogap.fa.gz'
    sampleB = data_dir / 'SampleB' / 'snps.aligned.nogap.fa.gz'

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        genomes_files = [
            ('SampleA', [sampleA]),
            ('SampleB', [sampleB])
        ]

        kmer_indexer = KmerIndexerSourmash(k=31, scaled=1000, compress=True)
        kmer_index_manager = KmerIndexManager(tmp_path, kmer_indexer=kmer_indexer)
        indexed_files = kmer_index_manager.index_all_genomes(genomes_files)

        assert {'SampleA', 'SampleB'} == set(indexed_files.keys())

        sigsA = load_file_as_signatures(str(indexed_files['SampleA']))
        valuesA = get_values_from_signatures(sigsA)
        assert {31} == valuesA['ksize']
        assert {1000} == valuesA['scaled']

        sigsB = load_file_as_signatures(str(indexed_files['SampleB']))
        valuesB = get_values_from_signatures(sigsB)
        assert {31} == valuesB['ksize']
        assert {1000} == valuesB['scaled']
