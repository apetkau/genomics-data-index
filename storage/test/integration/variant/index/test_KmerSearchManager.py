import math

from storage.test.integration.variant import sourmash_dir
from storage.variant.index.KmerSearchManager import KmerSearchManagerSourmash

sigs = {
    'SampleA': sourmash_dir / 'SampleA.sig.gz',
    'SampleB': sourmash_dir / 'SampleB.sig.gz',
    'SampleC': sourmash_dir / 'SampleC.sig.gz',
}


def test_search_single_file():
    search_manager = KmerSearchManagerSourmash()

    results_df = search_manager.search(31, sigs['SampleA'], [sigs['SampleB'], sigs['SampleC']])
    print(results_df)

    assert 2 == len(results_df)
    assert ['SampleC', 'SampleB'] == list(results_df['name'].tolist())
    assert math.isclose(0.5, results_df['similarity'].tolist()[0], rel_tol=1e-3)
    assert math.isclose(0.478, results_df['similarity'].tolist()[1], rel_tol=1e-3)
