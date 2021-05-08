import math

from genomics_data_index.storage.index.KmerSearchManager import KmerSearchManagerSourmash
from genomics_data_index.test.integration import sourmash_dir

sigs = {
    'SampleA': sourmash_dir / 'SampleA.sig.gz',
    'SampleB': sourmash_dir / 'SampleB.sig.gz',
    'SampleC': sourmash_dir / 'SampleC.sig.gz',
}


def test_search_single_file():
    search_manager = KmerSearchManagerSourmash()

    results_df = search_manager.search(kmer_size=31, similarity_threshold=0.0,
                                       query_file=sigs['SampleA'], search_files=[sigs['SampleB'], sigs['SampleC']])
    print(results_df)

    assert 2 == len(results_df)
    assert ['SampleC', 'SampleB'] == list(results_df['name'].tolist())
    assert math.isclose(0.5, results_df['similarity'].tolist()[0], rel_tol=1e-3)
    assert math.isclose(0.478, results_df['similarity'].tolist()[1], rel_tol=1e-3)


def test_search_single_file_higher_threshold():
    search_manager = KmerSearchManagerSourmash()

    results_df = search_manager.search(kmer_size=31, similarity_threshold=0.49,
                                       query_file=sigs['SampleA'], search_files=[sigs['SampleB'], sigs['SampleC']])
    print(results_df)

    assert 1 == len(results_df)
    assert ['SampleC'] == list(results_df['name'].tolist())
    assert math.isclose(0.5, results_df['similarity'].tolist()[0], rel_tol=1e-3)


def test_search_single_file_higher_threshold_no_results():
    search_manager = KmerSearchManagerSourmash()

    results_df = search_manager.search(kmer_size=31, similarity_threshold=0.9,
                                       query_file=sigs['SampleA'], search_files=[sigs['SampleB'], sigs['SampleC']])
    print(results_df)

    assert 0 == len(results_df)


def test_search_single_file_2():
    search_manager = KmerSearchManagerSourmash()

    results_df = search_manager.search(kmer_size=31, similarity_threshold=0.0,
                                       query_file=sigs['SampleB'], search_files=[sigs['SampleA'], sigs['SampleC']])
    print(results_df)

    assert 2 == len(results_df)
    assert ['SampleC', 'SampleA'] == list(results_df['name'].tolist())
    assert math.isclose(0.681, results_df['similarity'].tolist()[0], rel_tol=1e-3)
    assert math.isclose(0.478, results_df['similarity'].tolist()[1], rel_tol=1e-3)
