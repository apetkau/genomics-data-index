from typing import cast

from genomics_data_index.api.query.GenomicsDataIndex import GenomicsDataIndex
from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.impl import ExperimentalTreeSamplesQuery
from genomics_data_index.configuration.connector import DataIndexConnection


# wrapper methods to simplify writing tests
def query(connection: DataIndexConnection, **kwargs) -> SamplesQuery:
    return GenomicsDataIndex(connection=connection).samples_query(**kwargs)


def test_query_pruning_tree(loaded_database_connection: DataIndexConnection):
    query_result = query(loaded_database_connection).build_tree(
        kind='mutation_experimental', scope='genome',
        include_reference=True, extra_params='--seed 42 -m GTR').reset_universe()
    assert 3 == len(query_result)
    assert 3 == len(query_result.universe_set)

    query_result = cast(ExperimentalTreeSamplesQuery, query_result)

    assert {'SampleA', 'SampleB', 'SampleC', 'genome'} == set(
        query_result.tree.get_leaf_names())
    assert {'SampleA', 'genome'} == set(
        query_result.isa('SampleA').tree.get_leaf_names())
    assert {'SampleA', 'SampleB', 'genome'} == set(
        query_result.isin(['SampleA', 'SampleB']).tree.get_leaf_names())
    assert {'SampleB', 'SampleC', 'genome'} == set(
        query_result.isin(['SampleB', 'SampleC']).tree.get_leaf_names())

    # Test complement
    assert {'SampleB', 'SampleC', 'genome'} == set(
        query_result.isin('SampleA').complement().tree.get_leaf_names())

    # Double complement
    assert {'SampleB', 'genome'} == set(
        query_result.isin('SampleB').complement().complement().tree.get_leaf_names())

    # Test reset universe
    assert {'SampleA', 'SampleB', 'SampleC', 'genome'} == set(
        query_result.isa('SampleA').universe_tree.get_leaf_names())
    assert {'SampleA', 'genome'} == set(
        query_result.isa('SampleA').reset_universe().universe_tree.get_leaf_names())

    # Test reset universe 2
    assert {'SampleA', 'SampleB', 'SampleC', 'genome'} == set(
        query_result.isa(['SampleB', 'SampleC']).universe_tree.get_leaf_names())
    assert {'SampleB', 'SampleC', 'genome'} == set(
        query_result.isa(['SampleB', 'SampleC']).reset_universe().universe_tree.get_leaf_names())

    # Test reset universe and complement
    assert {'genome'} == set(
        query_result.isin('SampleA').reset_universe().complement().tree.get_leaf_names())
