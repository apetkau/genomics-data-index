import math
import re

import pandas as pd
from ete3 import Tree

from genomics_data_index.api.query.GenomicsDataIndex import GenomicsDataIndex
from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.impl.ClusterScorer import ClusterScorer
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.storage.model.db import Sample


# wrapper methods to simplify writing tests
def query(connection: DataIndexConnection, **kwargs) -> SamplesQuery:
    return GenomicsDataIndex(connection=connection).samples_query(**kwargs)


def test_score_samples_with_mutation_tree(prebuilt_tree: Tree, loaded_database_connection: DataIndexConnection):
    query_tree = query(loaded_database_connection).join_tree(tree=prebuilt_tree, kind='mutation',
                                                             reference_name='genome',
                                                             alignment_length=5180)
    assert 3 == len(query_tree)
    assert 9 == len(query_tree.universe_set)

    # Test single sample
    query_result = query_tree.isa('SampleA')
    assert 1 == len(query_result)
    cluster_scorer = ClusterScorer.create(query_tree)
    assert math.isclose(1 / 1, cluster_scorer.score_samples(query_result), rel_tol=1e-3)

    # Test single sample set
    query_result = query_tree.isa('SampleA')
    assert 1 == len(query_result)
    cluster_scorer = ClusterScorer.create(query_tree)
    assert math.isclose(1 / 1, cluster_scorer.score_samples(query_result.sample_set), rel_tol=1e-3)

    # Test multiple samples
    query_result = query_tree.isin(['SampleB', 'SampleC'])
    assert 2 == len(query_result)
    cluster_scorer = ClusterScorer.create(query_tree)
    assert math.isclose(2 / 2, cluster_scorer.score_samples(query_result), rel_tol=1e-3)

    # Test multiple samples set
    query_result = query_tree.isin(['SampleB', 'SampleC'])
    assert 2 == len(query_result)
    cluster_scorer = ClusterScorer.create(query_tree)
    assert math.isclose(2 / 2, cluster_scorer.score_samples(query_result.sample_set), rel_tol=1e-3)

    # Test multiple samples non-perfect score
    query_result = query_tree.isin(['SampleA', 'SampleB'])
    assert 2 == len(query_result)
    cluster_scorer = ClusterScorer.create(query_tree)
    assert math.isclose(2 / 3, cluster_scorer.score_samples(query_result), rel_tol=1e-3)

    # Test multiple samples set non-perfect score
    query_result = query_tree.isin(['SampleA', 'SampleB'])
    assert 2 == len(query_result)
    cluster_scorer = ClusterScorer.create(query_tree)
    assert math.isclose(2 / 3, cluster_scorer.score_samples(query_result.sample_set), rel_tol=1e-3)

    # Test multiple samples non-perfect score 2
    query_result = query_tree.isin(['SampleA', 'SampleC'])
    assert 2 == len(query_result)
    cluster_scorer = ClusterScorer.create(query_tree)
    assert math.isclose(2 / 3, cluster_scorer.score_samples(query_result), rel_tol=1e-3)

    # Test all samples
    query_result = query_tree.isin(['SampleA', 'SampleB', 'SampleC'])
    assert 3 == len(query_result)
    cluster_scorer = ClusterScorer.create(query_tree)
    assert math.isclose(3 / 3, cluster_scorer.score_samples(query_result), rel_tol=1e-3)

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


def test_score_groupby_with_mutation_tree(prebuilt_tree: Tree, loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    metadata_df = pd.DataFrame([
        [sampleA.id, 'red', 'up', 'A', '1', 'red', 'A'],
        [sampleB.id, 'blue', 'up', 'B', '1', 'red', 'A.1'],
        [sampleC.id, 'blue', 'down', 'C', '1', pd.NA, 'A.1.1']
    ], columns=['Sample ID', 'Color', 'Direction', 'Letter', 'Number', 'Color_NA', 'Type'])

    query_tree = query(loaded_database_connection).join_tree(tree=prebuilt_tree, kind='mutation',
                                                             reference_name='genome',
                                                             alignment_length=5180)
    query_tree_df = query_tree.join(metadata_df, sample_ids_column='Sample ID')
    cluster_scorer = ClusterScorer.create(query_tree_df)

    # Group by 2 groups perfect score
    scores_df = cluster_scorer.score_groupby('Color')
    assert ['Score', 'Sample Count'] == scores_df.columns.tolist()
    assert 2 == len(scores_df)
    assert {'red', 'blue'} == set(scores_df.index)
    assert math.isclose(1 / 1, scores_df.loc['red', 'Score'], abs_tol=1e-3)
    assert math.isclose(2 / 2, scores_df.loc['blue', 'Score'], abs_tol=1e-3)

    # Group by 2 groups imperfect score
    scores_df = cluster_scorer.score_groupby('Direction')
    assert ['Score', 'Sample Count'] == scores_df.columns.tolist()
    assert 2 == len(scores_df)
    assert {'up', 'down'} == set(scores_df.index)
    assert math.isclose(2 / 3, scores_df.loc['up', 'Score'], abs_tol=1e-3)
    assert math.isclose(1 / 1, scores_df.loc['down', 'Score'], abs_tol=1e-3)

    # Group by 3 groups
    scores_df = cluster_scorer.score_groupby('Letter')
    assert ['Score', 'Sample Count'] == scores_df.columns.tolist()
    assert 3 == len(scores_df)
    assert {'A', 'B', 'C'} == set(scores_df.index)
    assert math.isclose(1 / 1, scores_df.loc['A', 'Score'], abs_tol=1e-3)
    assert math.isclose(1 / 1, scores_df.loc['B', 'Score'], abs_tol=1e-3)
    assert math.isclose(1 / 1, scores_df.loc['C', 'Score'], abs_tol=1e-3)

    # Group by 1 group
    scores_df = cluster_scorer.score_groupby('Number')
    assert ['Score', 'Sample Count'] == scores_df.columns.tolist()
    assert 1 == len(scores_df)
    assert {'1'} == set(scores_df.index)
    assert math.isclose(1 / 1, scores_df.loc['1', 'Score'], abs_tol=1e-3)

    # Group by with NA exclude NA
    scores_df = cluster_scorer.score_groupby('Color_NA')
    assert ['Score', 'Sample Count'] == scores_df.columns.tolist()
    assert 1 == len(scores_df)
    assert {'red'} == set(scores_df.index)
    assert math.isclose(2 / 3, scores_df.loc['red', 'Score'], abs_tol=1e-3)

    # Group by with NA, fill na
    scores_df = cluster_scorer.score_groupby('Color_NA', na_value='NA')
    assert ['Score', 'Sample Count'] == scores_df.columns.tolist()
    assert 2 == len(scores_df)
    assert {'red', 'NA'} == set(scores_df.index)
    assert math.isclose(2 / 3, scores_df.loc['red', 'Score'], abs_tol=1e-3)
    assert math.isclose(1 / 1, scores_df.loc['NA', 'Score'], abs_tol=1e-3)

    # Group by 2 groups exclude 1 due to minimum samples
    scores_df = cluster_scorer.score_groupby('Color', min_samples_count=2)
    assert ['Score', 'Sample Count'] == scores_df.columns.tolist()
    assert 1 == len(scores_df)
    assert {'blue'} == set(scores_df.index)
    assert math.isclose(2 / 2, scores_df.loc['blue', 'Score'], abs_tol=1e-3)

    # Group by 2 groups exclude 1 due to maximum samples
    scores_df = cluster_scorer.score_groupby('Color', max_samples_count=1)
    assert ['Score', 'Sample Count'] == scores_df.columns.tolist()
    assert 1 == len(scores_df)
    assert {'red'} == set(scores_df.index)
    assert math.isclose(1 / 1, scores_df.loc['red', 'Score'], abs_tol=1e-3)

    # Group by 3 groups using Type
    scores_df = cluster_scorer.score_groupby('Type')
    assert ['Score', 'Sample Count'] == scores_df.columns.tolist()
    assert 3 == len(scores_df)
    assert {'A', 'A.1', 'A.1.1'} == set(scores_df.index)
    assert math.isclose(1 / 1, scores_df.loc['A', 'Score'], abs_tol=1e-3)
    assert math.isclose(1 / 1, scores_df.loc['A.1', 'Score'], abs_tol=1e-3)
    assert math.isclose(1 / 1, scores_df.loc['A.1.1', 'Score'], abs_tol=1e-3)

    # Group by A.1 group using Type and groupby_func
    scores_df = cluster_scorer.score_groupby('Type', groupby_func=
    lambda x: m.group(1) if (m := re.search(r'^([^.]+\.[^.]+)', x)) else pd.NA)
    assert ['Score', 'Sample Count'] == scores_df.columns.tolist()
    assert 1 == len(scores_df)
    assert {'A.1'} == set(scores_df.index)
    assert math.isclose(2 / 2, scores_df.loc['A.1', 'Score'], abs_tol=1e-3)

    # Group A.1 group using Type and groupby_func
    scores_df = cluster_scorer.score_groupby('Type', groupby_func=
    lambda x: m.group(1) if (m := re.search(r'^([^.]+)', x)) else pd.NA)
    assert ['Score', 'Sample Count'] == scores_df.columns.tolist()
    assert 1 == len(scores_df)
    assert {'A'} == set(scores_df.index)
    assert math.isclose(3 / 3, scores_df.loc['A', 'Score'], abs_tol=1e-3)

    # Group by 2 groups using Type and groupby_func
    scores_df = cluster_scorer.score_groupby('Type', groupby_func=
    lambda x: m.group(1) if (m := re.search(r'^([^.]+\.[^.]+)', x)) else x)
    assert ['Score', 'Sample Count'] == scores_df.columns.tolist()
    assert 2 == len(scores_df)
    assert {'A', 'A.1'} == set(scores_df.index)
    assert math.isclose(1 / 1, scores_df.loc['A', 'Score'], abs_tol=1e-3)
    assert math.isclose(2 / 2, scores_df.loc['A.1', 'Score'], abs_tol=1e-3)
