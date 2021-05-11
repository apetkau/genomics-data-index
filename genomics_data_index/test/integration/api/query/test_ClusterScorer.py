import math
import pandas as pd

from genomics_data_index.api.query.GenomicsDataIndex import GenomicsDataIndex
from genomics_data_index.storage.model.db import Sample
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


def test_score_groupby_with_mutation_tree(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    metadata_df1 = pd.DataFrame([
        [sampleA.id, 'red'],
        [sampleB.id, 'blue'],
        [sampleC.id, 'blue']
    ], columns=['Sample ID', 'Color'])

    query_tree = query(loaded_database_connection).build_tree(kind='mutation',
                                                              scope='genome', include_reference=True,
                                                              extra_params='--seed 42 -m GTR')
    query_tree_df1 = query_tree.join(metadata_df1, sample_ids_column='Sample ID')

    # Group by on color
    cluster_scorer = ClusterScorer.create(query_tree_df1)
    scores_series = cluster_scorer.score_groupby('Color')
    assert 2 == len(scores_series)
    assert {'red', 'blue'} == set(scores_series.index)
    assert math.isclose(1/1, scores_series.loc['red'], abs_tol=1e-3)
    assert math.isclose(2/2, scores_series.loc['blue'], abs_tol=1e-3)
