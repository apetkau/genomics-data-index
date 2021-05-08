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
    assert 2 == len(results_df)
    assert ['SampleC', 'SampleB'] == list(results_df['name'].tolist())
    assert math.isclose(0.5, results_df['similarity'].tolist()[0], rel_tol=1e-3)
    assert math.isclose(0.478, results_df['similarity'].tolist()[1], rel_tol=1e-3)


def test_search_single_file_higher_threshold():
    search_manager = KmerSearchManagerSourmash()

    results_df = search_manager.search(kmer_size=31, similarity_threshold=0.49,
                                       query_file=sigs['SampleA'], search_files=[sigs['SampleB'], sigs['SampleC']])
    assert 1 == len(results_df)
    assert ['SampleC'] == list(results_df['name'].tolist())
    assert math.isclose(0.5, results_df['similarity'].tolist()[0], rel_tol=1e-3)


def test_search_single_file_higher_threshold_no_results():
    search_manager = KmerSearchManagerSourmash()

    results_df = search_manager.search(kmer_size=31, similarity_threshold=0.9,
                                       query_file=sigs['SampleA'], search_files=[sigs['SampleB'], sigs['SampleC']])
    assert 0 == len(results_df)


def test_search_single_file_2():
    search_manager = KmerSearchManagerSourmash()

    results_df = search_manager.search(kmer_size=31, similarity_threshold=0.0,
                                       query_file=sigs['SampleB'], search_files=[sigs['SampleA'], sigs['SampleC']])
    assert 2 == len(results_df)
    assert ['SampleC', 'SampleA'] == list(results_df['name'].tolist())
    assert math.isclose(0.681, results_df['similarity'].tolist()[0], rel_tol=1e-3)
    assert math.isclose(0.478, results_df['similarity'].tolist()[1], rel_tol=1e-3)


def test_distances_all():
    search_manager = KmerSearchManagerSourmash()

    results_d, labels = search_manager.distances(kmer_size=31,
                                                 signature_files=[sigs['SampleA'], sigs['SampleB'], sigs['SampleC']])
    assert (3, 3) == results_d.shape
    assert ['SampleA', 'SampleB', 'SampleC'] == labels

    l = {element: idx for idx, element in enumerate(labels)}

    assert math.isclose(results_d[l['SampleA']][l['SampleA']], 0, rel_tol=1e-3)
    assert math.isclose(results_d[l['SampleA']][l['SampleB']], 0.522, rel_tol=1e-3)
    assert math.isclose(results_d[l['SampleA']][l['SampleC']], 0.5, rel_tol=1e-3)

    assert math.isclose(results_d[l['SampleB']][l['SampleA']], 0.522, rel_tol=1e-3)
    assert math.isclose(results_d[l['SampleB']][l['SampleB']], 0, rel_tol=1e-3)
    assert math.isclose(results_d[l['SampleB']][l['SampleC']], 0.3186, rel_tol=1e-3)

    assert math.isclose(results_d[l['SampleC']][l['SampleA']], 0.5, rel_tol=1e-3)
    assert math.isclose(results_d[l['SampleC']][l['SampleB']], 0.3186, rel_tol=1e-3)
    assert math.isclose(results_d[l['SampleC']][l['SampleC']], 0, rel_tol=1e-3)


def test_distances_pair():
    search_manager = KmerSearchManagerSourmash()

    results_d, labels = search_manager.distances(kmer_size=31,
                                                 signature_files=[sigs['SampleA'], sigs['SampleC']])
    assert (2, 2) == results_d.shape
    assert ['SampleA', 'SampleC'] == labels

    l = {element: idx for idx, element in enumerate(labels)}

    assert math.isclose(results_d[l['SampleA']][l['SampleA']], 0, rel_tol=1e-3)
    assert math.isclose(results_d[l['SampleA']][l['SampleC']], 0.5, rel_tol=1e-3)

    assert math.isclose(results_d[l['SampleC']][l['SampleA']], 0.5, rel_tol=1e-3)
    assert math.isclose(results_d[l['SampleC']][l['SampleC']], 0, rel_tol=1e-3)


def test_distances_single():
    search_manager = KmerSearchManagerSourmash()

    results_d, labels = search_manager.distances(kmer_size=31,
                                                 signature_files=[sigs['SampleA']])
    assert (1, 1) == results_d.shape
    assert ['SampleA'] == labels

    l = {element: idx for idx, element in enumerate(labels)}

    assert math.isclose(results_d[l['SampleA']][l['SampleA']], 0, rel_tol=1e-3)
