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
                                                              scope='genome', include_reference=True)
    assert 3 == len(query_tree)
    assert 9 == len(query_tree.universe_set)

    # Test single sample
    query_result = query_tree.isa('SampleA')
    assert 1 == len(query_result)
    assert 9 == len(query_result.universe_set)

    cluster_scorer = ClusterScorer.create(query_tree)
    assert math.isclose(1.0, cluster_scorer.score_samples(query_result), rel_tol=1e-3)

    # # Test no samples
    # query_result = query_tree.isa('not exists')
    # assert 0 == len(query_result)
    # assert 9 == len(query_result.universe_set)
    #
    # cluster_scorer = ClusterScorer.create(query_tree)
    # assert math.isclose(0.0, cluster_scorer.score_samples(query_result), rel_tol=1e-3)
