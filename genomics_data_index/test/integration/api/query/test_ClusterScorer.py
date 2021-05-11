import math

from genomics_data_index.api.query.GenomicsDataIndex import GenomicsDataIndex

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.impl.ClusterScorer import ClusterScorer
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection


# wrapper methods to simplify writing tests
def query(connection: DataIndexConnection, **kwargs) -> SamplesQuery:
    return GenomicsDataIndex(connection=connection).samples_query(**kwargs)


def test_score_samples_with_mutation_tree(loaded_database_connection: DataIndexConnection):
    query_tree = query(loaded_database_connection).build_tree(kind='mutation',
                                                              scope='genome', include_reference=True,
                                                              extra_params='--seed 42 -m GTR')
    assert 3 == len(query_tree)
    assert 9 == len(query_tree.universe_set)

    # Test single sample
    query_result = query_tree.isa('SampleA')
    assert 1 == len(query_result)
    cluster_scorer = ClusterScorer.create(query_tree)
    assert math.isclose(1/1, cluster_scorer.score_samples(query_result), rel_tol=1e-3)

    # Test single sample set
    query_result = query_tree.isa('SampleA')
    assert 1 == len(query_result)
    cluster_scorer = ClusterScorer.create(query_tree)
    assert math.isclose(1/1, cluster_scorer.score_samples(query_result.sample_set), rel_tol=1e-3)

    # Test multiple samples
    query_result = query_tree.isin(['SampleB', 'SampleC'])
    assert 2 == len(query_result)
    cluster_scorer = ClusterScorer.create(query_tree)
    assert math.isclose(2/2, cluster_scorer.score_samples(query_result), rel_tol=1e-3)

    # Test multiple samples set
    query_result = query_tree.isin(['SampleB', 'SampleC'])
    assert 2 == len(query_result)
    cluster_scorer = ClusterScorer.create(query_tree)
    assert math.isclose(2/2, cluster_scorer.score_samples(query_result.sample_set), rel_tol=1e-3)

    # Test multiple samples non-perfect score
    query_result = query_tree.isin(['SampleA', 'SampleB'])
    assert 2 == len(query_result)
    cluster_scorer = ClusterScorer.create(query_tree)
    assert math.isclose(2/3, cluster_scorer.score_samples(query_result), rel_tol=1e-3)

    # Test multiple samples set non-perfect score
    query_result = query_tree.isin(['SampleA', 'SampleB'])
    assert 2 == len(query_result)
    cluster_scorer = ClusterScorer.create(query_tree)
    assert math.isclose(2/3, cluster_scorer.score_samples(query_result.sample_set), rel_tol=1e-3)

    # Test multiple samples non-perfect score 2
    query_result = query_tree.isin(['SampleA', 'SampleC'])
    assert 2 == len(query_result)
    cluster_scorer = ClusterScorer.create(query_tree)
    assert math.isclose(2/3, cluster_scorer.score_samples(query_result), rel_tol=1e-3)

    # Test all samples
    query_result = query_tree.isin(['SampleA', 'SampleB', 'SampleC'])
    assert 3 == len(query_result)
    cluster_scorer = ClusterScorer.create(query_tree)
    assert math.isclose(3/3, cluster_scorer.score_samples(query_result), rel_tol=1e-3)

    # # Test no samples
    query_result = query_tree.isa('not exists')
    assert 0 == len(query_result)
    cluster_scorer = ClusterScorer.create(query_tree)
    assert math.isnan(cluster_scorer.score_samples(query_result))

    # # Test no samples set
    query_result = query_tree.isa('not exists')
    assert 0 == len(query_result)
    cluster_scorer = ClusterScorer.create(query_tree)
    assert math.isnan(cluster_scorer.score_samples(query_result.sample_set))
