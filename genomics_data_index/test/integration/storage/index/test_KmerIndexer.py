import tempfile
from pathlib import Path

from sourmash import load_file_as_signatures

from genomics_data_index.storage.index.KmerIndexer import KmerIndexerSourmash
from genomics_data_index.test.integration import data_dir
from genomics_data_index.test.integration.storage.index import get_values_from_signatures


def test_index_single_file_uncompressed():
    sampleA = data_dir / 'SampleA' / 'snps.aligned.nogap.fa.gz'

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        index_out = tmp_path / 'SampleA.sig'

        kmer_indexer = KmerIndexerSourmash(k=31, scaled=1000, compress=False)
        indexed_path = kmer_indexer.index('SampleA', index_out, [sampleA])

        assert indexed_path == index_out
        assert indexed_path.exists()
        sigs = load_file_as_signatures(str(indexed_path))
        values = get_values_from_signatures(sigs)

        assert {31} == values['ksize']
        assert {1000} == values['scaled']


def test_index_single_file_compressed():
    sampleA = data_dir / 'SampleA' / 'snps.aligned.nogap.fa.gz'

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        index_out = tmp_path / 'SampleA.sig'

        kmer_indexer = KmerIndexerSourmash(k=31, scaled=1000, compress=True)
        indexed_path = kmer_indexer.index('SampleA', index_out, [sampleA])

        assert '.gz' == indexed_path.suffix
        assert indexed_path.exists()
        sigs = load_file_as_signatures(str(indexed_path))
        values = get_values_from_signatures(sigs)

        assert {31} == values['ksize']
        assert {1000} == values['scaled']


def test_index_single_file_multiple_k():
    sampleA = data_dir / 'SampleA' / 'snps.aligned.nogap.fa.gz'

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir)
        index_out = tmp_path / 'SampleA.sig.gz'

        kmer_indexer = KmerIndexerSourmash(k=[21, 31], scaled=1000, compress=True)
        indexed_path = kmer_indexer.index('SampleA', index_out, [sampleA])

        assert indexed_path.exists()
        sigs = load_file_as_signatures(str(indexed_path))
        values = get_values_from_signatures(sigs)

        assert {21, 31} == values['ksize']
        assert {1000} == values['scaled']
