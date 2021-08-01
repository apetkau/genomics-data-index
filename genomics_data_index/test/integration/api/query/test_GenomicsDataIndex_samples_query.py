import math
import warnings
from typing import cast

import pandas as pd
import pytest
from ete3 import Tree, ClusterTree

from genomics_data_index.api.query.GenomicsDataIndex import GenomicsDataIndex
from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.impl.DataFrameSamplesQuery import DataFrameSamplesQuery
from genomics_data_index.api.query.impl.ExperimentalTreeSamplesQuery import ExperimentalTreeSamplesQuery
from genomics_data_index.api.query.impl.MutationTreeSamplesQuery import MutationTreeSamplesQuery
from genomics_data_index.api.query.impl.TreeSamplesQuery import TreeSamplesQuery
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.io.mutation.NucleotideSampleDataPackage import NucleotideSampleDataPackage
from genomics_data_index.storage.model.QueryFeatureHGVS import QueryFeatureHGVS
from genomics_data_index.storage.model.QueryFeatureMLST import QueryFeatureMLST
from genomics_data_index.storage.model.QueryFeatureMutationSPDI import QueryFeatureMutationSPDI
from genomics_data_index.storage.model.db import Sample
from genomics_data_index.test.integration import snippy_all_dataframes, data_dir, snpeff_tree_file


# wrapper methods to simplify writing tests
def query(connection: DataIndexConnection, **kwargs) -> SamplesQuery:
    return GenomicsDataIndex(connection=connection).samples_query(**kwargs)


def test_initialized_query_default(loaded_database_connection: DataIndexConnection):
    initial_query = query(loaded_database_connection)

    assert len(initial_query) == 9
    assert len(initial_query.universe_set) == 9
    assert len(initial_query.sample_set) == 9


def test_initialized_query_mutations(loaded_database_connection_with_built_tree: DataIndexConnection):
    initial_query = query(loaded_database_connection_with_built_tree)
    assert len(initial_query) == 9
    assert len(initial_query.sample_set) == 9
    assert len(initial_query.universe_set) == 9

    initial_query = query(loaded_database_connection_with_built_tree,
                          universe='mutations', reference_name='genome')
    assert len(initial_query) == 3
    assert len(initial_query.sample_set) == 3
    assert len(initial_query.universe_set) == 3
    assert isinstance(initial_query, TreeSamplesQuery)
    assert initial_query.tree is not None


def test_empty_universe(loaded_database_connection: DataIndexConnection):
    query_result = query(loaded_database_connection).isin('empty')
    assert 0 == len(query_result)
    assert 0 == len(query_result.unknown_set)
    assert 9 == len(query_result.absent_set)
    assert 9 == len(query_result.universe_set)

    query_result = query_result.reset_universe()
    assert 0 == len(query_result)
    assert 0 == len(query_result.unknown_set)
    assert 0 == len(query_result.absent_set)
    assert 0 == len(query_result.universe_set)

    assert '<SamplesQueryIndex[' in str(query_result)
    assert 'selected=NA%' in str(query_result)
    assert 'unknown=NA%' in str(query_result)


def test_query_isin_samples(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    query_result = query(loaded_database_connection).isin('SampleB')
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 0 == len(query_result.unknown_set)
    assert 9 == len(query_result.universe_set)


def test_query_isin_samples_no_exist(loaded_database_connection: DataIndexConnection):
    query_result = query(loaded_database_connection).isin('no_exist')
    assert 0 == len(query_result)
    assert 0 == len(query_result.unknown_set)
    assert 9 == len(query_result.universe_set)


def test_query_isin_sample_set_single(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleBSet = SampleSet([sampleB.id])

    query_result = query(loaded_database_connection).isin(sampleBSet)
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 0 == len(query_result.unknown_set)
    assert 9 == len(query_result.universe_set)


def test_query_isin_samples_query_single(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    query_result_B = query(loaded_database_connection).isin('SampleB')

    query_result = query(loaded_database_connection).isin(query_result_B)
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 0 == len(query_result.unknown_set)
    assert 9 == len(query_result.universe_set)


def test_query_isin_samples_query_no_matches(loaded_database_connection: DataIndexConnection):
    query_result_empty = query(loaded_database_connection).isin('no_exist')

    query_result = query(loaded_database_connection).isin(query_result_empty)
    assert 0 == len(query_result)
    assert 0 == len(query_result.unknown_set)
    assert 9 == len(query_result.universe_set)


def test_query_isin_samples_with_unknown(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    query_result = query(loaded_database_connection).hasa('reference:5061:G:A')

    assert {sampleB.id} == set(query_result.sample_set)
    assert {sampleA.id} == set(query_result.unknown_set)

    # case isin with present and unknown
    query_result_isin = query_result.isin(['SampleA', 'SampleB'])
    assert 1 == len(query_result_isin)
    assert {sampleB.id} == set(query_result_isin.sample_set)
    assert {sampleA.id} == set(query_result_isin.unknown_set)
    assert 9 == len(query_result_isin.universe_set)

    # isin with only unknown
    query_result_isin = query_result.isin('SampleA')
    assert 0 == len(query_result_isin)
    assert set() == set(query_result_isin.sample_set)
    assert {sampleA.id} == set(query_result_isin.unknown_set)
    assert 9 == len(query_result_isin.universe_set)

    # isin with only present
    query_result_isin = query_result.isin('SampleB')
    assert 1 == len(query_result_isin)
    assert {sampleB.id} == set(query_result_isin.sample_set)
    assert set() == set(query_result_isin.unknown_set)
    assert 9 == len(query_result_isin.universe_set)


def test_query_and_or(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    query_resultA = query(loaded_database_connection).isa('SampleA')
    assert 1 == len(query_resultA)
    assert {sampleA.id} == set(query_resultA.sample_set)
    assert 0 == len(query_resultA.unknown_set)
    assert 9 == len(query_resultA.universe_set)

    query_result_onlyA = query_resultA.reset_universe()
    assert 1 == len(query_result_onlyA)
    assert {sampleA.id} == set(query_result_onlyA.sample_set)
    assert 0 == len(query_result_onlyA.unknown_set)
    assert 1 == len(query_result_onlyA.universe_set)

    query_resultB = query(loaded_database_connection).isa('SampleB')
    assert 1 == len(query_resultB)
    assert {sampleB.id} == set(query_resultB.sample_set)
    assert 0 == len(query_resultB.unknown_set)
    assert 9 == len(query_resultB.universe_set)

    query_result_onlyB = query_resultB.reset_universe()
    assert 1 == len(query_result_onlyB)
    assert {sampleB.id} == set(query_result_onlyB.sample_set)
    assert 0 == len(query_result_onlyB.unknown_set)
    assert 1 == len(query_result_onlyB.universe_set)

    query_resultB_with_unknownA = query(loaded_database_connection).hasa('reference:5061:G:A')
    assert 1 == len(query_resultB_with_unknownA)
    assert {sampleB.id} == set(query_resultB_with_unknownA.sample_set)
    assert 1 == len(query_resultB_with_unknownA.unknown_set)
    assert {sampleA.id} == set(query_resultB_with_unknownA.unknown_set)
    assert 9 == len(query_resultB_with_unknownA.universe_set)

    query_result_onlyB_with_unknownA = query_resultB_with_unknownA.reset_universe(include_unknown=True)
    assert 1 == len(query_result_onlyB_with_unknownA)
    assert {sampleB.id} == set(query_result_onlyB_with_unknownA.sample_set)
    assert 1 == len(query_result_onlyB_with_unknownA.unknown_set)
    assert {sampleA.id} == set(query_result_onlyB_with_unknownA.unknown_set)
    assert 2 == len(query_result_onlyB_with_unknownA.universe_set)
    assert {sampleA.id, sampleB.id} == set(query_result_onlyB_with_unknownA.universe_set)

    query_resultC = query(loaded_database_connection).isa('SampleC')
    assert 1 == len(query_resultC)
    assert {sampleC.id} == set(query_resultC.sample_set)
    assert 9 == len(query_resultC.universe_set)

    query_result_onlyC = query_resultC.reset_universe()
    assert 1 == len(query_result_onlyC)
    assert {sampleC.id} == set(query_result_onlyC.sample_set)
    assert 1 == len(query_result_onlyC.universe_set)
    assert {sampleC.id} == set(query_result_onlyC.universe_set)

    query_resultAB = query(loaded_database_connection).isin(['SampleA', 'SampleB'])
    assert 2 == len(query_resultAB)
    assert {sampleA.id, sampleB.id} == set(query_resultAB.sample_set)
    assert 0 == len(query_resultC.unknown_set)
    assert 9 == len(query_resultAB.universe_set)

    # Test AND
    query_resultA_and_AB = query_resultA.and_(query_resultAB)
    assert 1 == len(query_resultA_and_AB)
    assert {sampleA.id} == set(query_resultA_and_AB.sample_set)
    assert 0 == len(query_resultA_and_AB.unknown_set)
    assert 9 == len(query_resultA_and_AB.universe_set)

    query_resultAB_and_A = query_resultAB.and_(query_resultA)
    assert 1 == len(query_resultAB_and_A)
    assert {sampleA.id} == set(query_resultAB_and_A.sample_set)
    assert 0 == len(query_resultAB_and_A.unknown_set)
    assert 9 == len(query_resultAB_and_A.universe_set)

    query_resultA_and_B = query_resultA.and_(query_resultB)
    assert 0 == len(query_resultA_and_B)
    assert 0 == len(query_resultA_and_B.unknown_set)
    assert 9 == len(query_resultA_and_B.universe_set)

    # Test AND with unknowns
    query_resultA_and_B_unknown_A = query_resultA.and_(query_resultB_with_unknownA)
    assert 0 == len(query_resultA_and_B_unknown_A)
    assert 1 == len(query_resultA_and_B_unknown_A.unknown_set)
    assert {sampleA.id} == set(query_resultA_and_B_unknown_A.unknown_set)
    assert 9 == len(query_resultA_and_B_unknown_A.universe_set)

    query_resultAB_and_B_unknown_A = query_resultAB.and_(query_resultB_with_unknownA)
    assert 1 == len(query_resultAB_and_B_unknown_A)
    assert {sampleB.id} == set(query_resultAB_and_B_unknown_A.sample_set)
    assert 1 == len(query_resultAB_and_B_unknown_A.unknown_set)
    assert {sampleA.id} == set(query_resultAB_and_B_unknown_A.unknown_set)
    assert 9 == len(query_resultAB_and_B_unknown_A.universe_set)

    query_resultB_unknown_A_and_AB = query_resultB_with_unknownA.and_(query_resultAB)
    assert 1 == len(query_resultB_unknown_A_and_AB)
    assert {sampleB.id} == set(query_resultB_unknown_A_and_AB.sample_set)
    assert 1 == len(query_resultB_unknown_A_and_AB.unknown_set)
    assert {sampleA.id} == set(query_resultB_unknown_A_and_AB.unknown_set)
    assert 9 == len(query_resultB_unknown_A_and_AB.universe_set)

    query_resultB_and_B_unknown_A = query_resultB.and_(query_resultB_with_unknownA)
    assert 1 == len(query_resultB_and_B_unknown_A)
    assert {sampleB.id} == set(query_resultB_and_B_unknown_A.sample_set)
    assert 0 == len(query_resultB_and_B_unknown_A.unknown_set)
    assert 9 == len(query_resultB_and_B_unknown_A.universe_set)

    query_resultB_unknown_A_and_B = query_resultB_with_unknownA.and_(query_resultB)
    assert 1 == len(query_resultB_unknown_A_and_B)
    assert {sampleB.id} == set(query_resultB_unknown_A_and_B.sample_set)
    assert 0 == len(query_resultB_unknown_A_and_B.unknown_set)
    assert 9 == len(query_resultB_unknown_A_and_B.universe_set)

    query_result_onlyB_and_onlyA = query_result_onlyB.and_(query_result_onlyA)
    assert 0 == len(query_result_onlyB_and_onlyA)
    assert 0 == len(query_result_onlyB_and_onlyA.unknown_set)
    assert {sampleA.id, sampleB.id} == set(query_result_onlyB_and_onlyA.absent_set)
    assert 2 == len(query_result_onlyB_and_onlyA.universe_set)
    assert {sampleA.id, sampleB.id} == set(query_result_onlyB_and_onlyA.universe_set)

    query_result_onlyB_unknown_A_and_onlyA = query_result_onlyB_with_unknownA.and_(query_result_onlyA)
    assert 0 == len(query_result_onlyB_unknown_A_and_onlyA)
    assert 1 == len(query_result_onlyB_unknown_A_and_onlyA.unknown_set)
    assert {sampleA.id} == set(query_result_onlyB_unknown_A_and_onlyA.unknown_set)
    assert 1 == len(query_result_onlyB_unknown_A_and_onlyA.absent_set)
    assert 2 == len(query_result_onlyB_unknown_A_and_onlyA.universe_set)
    assert {sampleA.id, sampleB.id} == set(query_result_onlyB_unknown_A_and_onlyA.universe_set)

    query_result_onlyB_unknown_A_and_onlyC = query_result_onlyB_with_unknownA.and_(query_result_onlyC)
    assert 0 == len(query_result_onlyB_unknown_A_and_onlyC)
    assert 0 == len(query_result_onlyB_unknown_A_and_onlyC.unknown_set)
    assert 3 == len(query_result_onlyB_unknown_A_and_onlyC.absent_set)
    assert 3 == len(query_result_onlyB_unknown_A_and_onlyC.universe_set)
    assert {sampleA.id, sampleB.id, sampleC.id} == set(query_result_onlyB_unknown_A_and_onlyC.universe_set)

    # Test Python & (bitwise and) operator
    query_result_onlyB_unknown_A_and_onlyA = query_result_onlyB_with_unknownA & query_result_onlyA
    assert 0 == len(query_result_onlyB_unknown_A_and_onlyA)
    assert 1 == len(query_result_onlyB_unknown_A_and_onlyA.unknown_set)
    assert {sampleA.id} == set(query_result_onlyB_unknown_A_and_onlyA.unknown_set)
    assert 1 == len(query_result_onlyB_unknown_A_and_onlyA.absent_set)
    assert 2 == len(query_result_onlyB_unknown_A_and_onlyA.universe_set)
    assert {sampleA.id, sampleB.id} == set(query_result_onlyB_unknown_A_and_onlyA.universe_set)

    query_result_onlyB_unknown_A_and_onlyC = query_result_onlyB_with_unknownA & query_result_onlyC
    assert 0 == len(query_result_onlyB_unknown_A_and_onlyC)
    assert 0 == len(query_result_onlyB_unknown_A_and_onlyC.unknown_set)
    assert 3 == len(query_result_onlyB_unknown_A_and_onlyC.absent_set)
    assert 3 == len(query_result_onlyB_unknown_A_and_onlyC.universe_set)
    assert {sampleA.id, sampleB.id, sampleC.id} == set(query_result_onlyB_unknown_A_and_onlyC.universe_set)

    # Test OR
    query_resultA_or_AB = query_resultA.or_(query_resultAB)
    assert 2 == len(query_resultA_or_AB)
    assert {sampleA.id, sampleB.id} == set(query_resultA_or_AB.sample_set)
    assert 7 == len(query_resultA_or_AB.absent_set)
    assert 0 == len(query_resultA_or_AB.unknown_set)
    assert 9 == len(query_resultA_or_AB.universe_set)

    query_resultAB_or_A = query_resultAB.or_(query_resultA)
    assert 2 == len(query_resultAB_or_A)
    assert {sampleA.id, sampleB.id} == set(query_resultAB_or_A.sample_set)
    assert 0 == len(query_resultAB_or_A.unknown_set)
    assert 7 == len(query_resultAB_or_A.absent_set)
    assert 9 == len(query_resultAB_or_A.universe_set)

    query_resultA_or_B = query_resultA.or_(query_resultB)
    assert 2 == len(query_resultA_or_B)
    assert {sampleA.id, sampleB.id} == set(query_resultA_or_B.sample_set)
    assert 0 == len(query_resultA_or_B.unknown_set)
    assert 7 == len(query_resultA_or_B.absent_set)
    assert 9 == len(query_resultA_or_B.universe_set)

    query_resultA_or_B_or_C = query_resultA.or_(query_resultB).or_(query_resultC)
    assert 3 == len(query_resultA_or_B_or_C)
    assert {sampleA.id, sampleB.id, sampleC.id} == set(query_resultA_or_B_or_C.sample_set)
    assert 0 == len(query_resultA_or_B_or_C.unknown_set)
    assert 6 == len(query_resultA_or_B_or_C.absent_set)
    assert 9 == len(query_resultA_or_B_or_C.universe_set)

    # Test OR with unknowns
    query_resultA_or_B_unknown_A = query_resultA.or_(query_resultB_with_unknownA)
    assert 2 == len(query_resultA_or_B_unknown_A)
    assert {sampleA.id, sampleB.id} == set(query_resultA_or_B_unknown_A.sample_set)
    assert 0 == len(query_resultA_or_B_unknown_A.unknown_set)
    assert 9 == len(query_resultA_or_B_unknown_A.universe_set)

    query_resultAB_or_B_unknown_A = query_resultAB.or_(query_resultB_with_unknownA)
    assert 2 == len(query_resultAB_or_B_unknown_A)
    assert {sampleA.id, sampleB.id} == set(query_resultAB_or_B_unknown_A.sample_set)
    assert 0 == len(query_resultAB_or_B_unknown_A.unknown_set)
    assert 9 == len(query_resultAB_or_B_unknown_A.universe_set)

    query_resultB_unknown_A_or_AB = query_resultB_with_unknownA.or_(query_resultAB)
    assert 2 == len(query_resultB_unknown_A_or_AB)
    assert {sampleA.id, sampleB.id} == set(query_resultB_unknown_A_or_AB.sample_set)
    assert 0 == len(query_resultB_unknown_A_or_AB.unknown_set)
    assert 9 == len(query_resultB_unknown_A_or_AB.universe_set)

    query_resultB_or_B_unknown_A = query_resultB.or_(query_resultB_with_unknownA)
    assert 1 == len(query_resultB_or_B_unknown_A)
    assert {sampleB.id} == set(query_resultB_or_B_unknown_A.sample_set)
    assert 1 == len(query_resultB_or_B_unknown_A.unknown_set)
    assert {sampleA.id} == set(query_resultB_or_B_unknown_A.unknown_set)
    assert 9 == len(query_resultB_or_B_unknown_A.universe_set)

    query_resultB_unknown_A_or_B = query_resultB_with_unknownA.or_(query_resultB)
    assert 1 == len(query_resultB_unknown_A_or_B)
    assert {sampleB.id} == set(query_resultB_unknown_A_or_B.sample_set)
    assert 1 == len(query_resultB_unknown_A_or_B.unknown_set)
    assert {sampleA.id} == set(query_resultB_unknown_A_or_B.unknown_set)
    assert 9 == len(query_resultB_unknown_A_or_B.universe_set)

    query_result_onlyB_or_onlyA = query_result_onlyB.or_(query_result_onlyA)
    assert 2 == len(query_result_onlyB_or_onlyA)
    assert 0 == len(query_result_onlyB_or_onlyA.unknown_set)
    assert {sampleA.id, sampleB.id} == set(query_result_onlyB_or_onlyA.sample_set)
    assert 2 == len(query_result_onlyB_or_onlyA.universe_set)
    assert {sampleA.id, sampleB.id} == set(query_result_onlyB_or_onlyA.universe_set)

    query_result_onlyB_unknown_A_or_onlyA = query_result_onlyB_with_unknownA.or_(query_result_onlyA)
    assert 2 == len(query_result_onlyB_unknown_A_or_onlyA)
    assert 0 == len(query_result_onlyB_unknown_A_or_onlyA.unknown_set)
    assert {sampleA.id, sampleB.id} == set(query_result_onlyB_unknown_A_or_onlyA.sample_set)
    assert 0 == len(query_result_onlyB_unknown_A_or_onlyA.absent_set)
    assert 2 == len(query_result_onlyB_unknown_A_or_onlyA.universe_set)
    assert {sampleA.id, sampleB.id} == set(query_result_onlyB_unknown_A_or_onlyA.universe_set)

    query_result_onlyB_unknown_A_or_onlyC = query_result_onlyB_with_unknownA.or_(query_result_onlyC)
    assert 2 == len(query_result_onlyB_unknown_A_or_onlyC)
    assert {sampleB.id, sampleC.id} == set(query_result_onlyB_unknown_A_or_onlyC.sample_set)
    assert 1 == len(query_result_onlyB_unknown_A_or_onlyC.unknown_set)
    assert {sampleA.id} == set(query_result_onlyB_unknown_A_or_onlyC.unknown_set)
    assert 0 == len(query_result_onlyB_unknown_A_or_onlyC.absent_set)
    assert 3 == len(query_result_onlyB_unknown_A_or_onlyC.universe_set)
    assert {sampleA.id, sampleB.id, sampleC.id} == set(query_result_onlyB_unknown_A_or_onlyC.universe_set)

    # Test Python | (bitwise or) operator
    query_result_onlyB_unknown_A_or_onlyA = query_result_onlyB_with_unknownA | query_result_onlyA
    assert 2 == len(query_result_onlyB_unknown_A_or_onlyA)
    assert 0 == len(query_result_onlyB_unknown_A_or_onlyA.unknown_set)
    assert {sampleA.id, sampleB.id} == set(query_result_onlyB_unknown_A_or_onlyA.sample_set)
    assert 0 == len(query_result_onlyB_unknown_A_or_onlyA.absent_set)
    assert 2 == len(query_result_onlyB_unknown_A_or_onlyA.universe_set)
    assert {sampleA.id, sampleB.id} == set(query_result_onlyB_unknown_A_or_onlyA.universe_set)

    query_result_onlyB_unknown_A_or_onlyC = query_result_onlyB_with_unknownA | query_result_onlyC
    assert 2 == len(query_result_onlyB_unknown_A_or_onlyC)
    assert {sampleB.id, sampleC.id} == set(query_result_onlyB_unknown_A_or_onlyC.sample_set)
    assert 1 == len(query_result_onlyB_unknown_A_or_onlyC.unknown_set)
    assert {sampleA.id} == set(query_result_onlyB_unknown_A_or_onlyC.unknown_set)
    assert 0 == len(query_result_onlyB_unknown_A_or_onlyC.absent_set)
    assert 3 == len(query_result_onlyB_unknown_A_or_onlyC.universe_set)
    assert {sampleA.id, sampleB.id, sampleC.id} == set(query_result_onlyB_unknown_A_or_onlyC.universe_set)


def test_query_reset_universe(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    # No unknowns
    query_result = query(loaded_database_connection).isin('SampleB')
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 0 == len(query_result.unknown_set)
    assert 9 == len(query_result.universe_set)

    query_result = query_result.reset_universe()
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 0 == len(query_result.unknown_set)
    assert 1 == len(query_result.universe_set)

    # With unknowns
    query_result = query(loaded_database_connection).hasa('reference:5061:G:A')
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 1 == len(query_result.unknown_set)
    assert {sampleA.id} == set(query_result.unknown_set)
    assert 9 == len(query_result.universe_set)

    query_result = query_result.reset_universe()
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 1 == len(query_result.unknown_set)
    assert {sampleA.id} == set(query_result.unknown_set)
    assert 2 == len(query_result.universe_set)
    assert {sampleA.id, sampleB.id} == set(query_result.universe_set)

    # With unknowns but excluding them
    query_result = query(loaded_database_connection).hasa('reference:5061:G:A')
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 1 == len(query_result.unknown_set)
    assert {sampleA.id} == set(query_result.unknown_set)
    assert 9 == len(query_result.universe_set)

    query_result = query_result.reset_universe(include_unknown=False)
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 0 == len(query_result.unknown_set)
    assert 1 == len(query_result.universe_set)
    assert {sampleB.id} == set(query_result.universe_set)


def test_query_set_universe(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    # No unknowns
    query_result = query(loaded_database_connection).isin('SampleB')
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 0 == len(query_result.unknown_set)
    assert 8 == len(query_result.absent_set)
    assert 9 == len(query_result.universe_set)

    query_result_su = query_result.set_universe(SampleSet(sample_ids=[sampleA.id, sampleB.id]))
    assert 1 == len(query_result_su)
    assert {sampleB.id} == set(query_result_su.sample_set)
    assert 0 == len(query_result_su.unknown_set)
    assert 1 == len(query_result_su.absent_set)
    assert {sampleA.id} == set(query_result_su.absent_set)
    assert 2 == len(query_result_su.universe_set)
    assert {sampleA.id, sampleB.id} == set(query_result_su.universe_set)

    query_result_su = query_result.set_universe(SampleSet(sample_ids=[sampleB.id]))
    assert 1 == len(query_result_su)
    assert {sampleB.id} == set(query_result_su.sample_set)
    assert 0 == len(query_result_su.unknown_set)
    assert 0 == len(query_result_su.absent_set)
    assert 1 == len(query_result_su.universe_set)
    assert {sampleB.id} == set(query_result_su.universe_set)

    # With unknowns
    query_result = query(loaded_database_connection).hasa('reference:5061:G:A')
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 1 == len(query_result.unknown_set)
    assert {sampleA.id} == set(query_result.unknown_set)
    assert 7 == len(query_result.absent_set)
    assert 9 == len(query_result.universe_set)

    query_result_su = query_result.set_universe(SampleSet(sample_ids=[sampleA.id, sampleB.id]))
    assert 1 == len(query_result_su)
    assert {sampleB.id} == set(query_result_su.sample_set)
    assert 1 == len(query_result_su.unknown_set)
    assert 0 == len(query_result_su.absent_set)
    assert 2 == len(query_result_su.universe_set)
    assert {sampleA.id, sampleB.id} == set(query_result_su.universe_set)

    query_result_su = query_result.set_universe(SampleSet(sample_ids=[sampleB.id]))
    assert 1 == len(query_result_su)
    assert {sampleB.id} == set(query_result_su.sample_set)
    assert 0 == len(query_result_su.unknown_set)
    assert 0 == len(query_result_su.absent_set)
    assert 1 == len(query_result_su.universe_set)
    assert {sampleB.id} == set(query_result_su.universe_set)

    query_result_su = query_result.set_universe(SampleSet(sample_ids=[sampleA.id, sampleB.id, sampleC.id]))
    assert 1 == len(query_result_su)
    assert {sampleB.id} == set(query_result_su.sample_set)
    assert 1 == len(query_result_su.unknown_set)
    assert 1 == len(query_result_su.absent_set)
    assert 3 == len(query_result_su.universe_set)
    assert {sampleA.id, sampleB.id, sampleC.id} == set(query_result_su.universe_set)


def test_query_isin_samples_multilple(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    query_result = query(loaded_database_connection).isin(['SampleA', 'SampleB'])
    assert 2 == len(query_result)
    assert {sampleA.id, sampleB.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


def test_query_isin_samples_multilple_samples_query(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    query_result_AB = query(loaded_database_connection).isin(['SampleA', 'SampleB'])

    query_result = query(loaded_database_connection).isin(query_result_AB)
    assert 2 == len(query_result)
    assert {sampleA.id, sampleB.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


def test_query_isa_sample_name(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    query_result = query(loaded_database_connection).isa('SampleB')
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)
    assert not query_result.has_tree()

    # Test case where there is an unknown in the result
    query_result = query(loaded_database_connection).hasa('reference:1:1:T').isa('SampleB')
    assert 0 == len(query_result)
    assert {sampleB.id} == set(query_result.unknown_set)
    assert 8 == len(query_result.absent_set)
    assert 9 == len(query_result.universe_set)

    # Test case where there is a match and unknown in the result
    query_result = query(loaded_database_connection).hasa('reference:5061:G:A').isa('SampleB')
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 0 == len(query_result.unknown_set)
    assert 8 == len(query_result.absent_set)
    assert 9 == len(query_result.universe_set)

    # Test isa no exist
    query_result = query(loaded_database_connection).isa('no_exist')
    assert 0 == len(query_result)
    assert 9 == len(query_result.universe_set)
    assert 0 == len(query_result.unknown_set)

    # Test isa no exist with unknown
    query_result = query(loaded_database_connection).hasa('reference:5061:G:A').isa('no_exist')
    assert 0 == len(query_result)
    assert 0 == len(query_result.unknown_set)
    assert 9 == len(query_result.absent_set)
    assert 9 == len(query_result.universe_set)


def test_query_isa_samples_query(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    query_result_B = query(loaded_database_connection).isa('SampleB')

    query_result = query(loaded_database_connection).isa(query_result_B)
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 0 == len(query_result.unknown_set)
    assert 8 == len(query_result.absent_set)
    assert 9 == len(query_result.universe_set)

    # with unknown in result
    query_result = query(loaded_database_connection).hasa('reference:1:1:T').isa(query_result_B)
    assert 0 == len(query_result)
    assert 1 == len(query_result.unknown_set)
    assert {sampleB.id} == set(query_result.unknown_set)
    assert 8 == len(query_result.absent_set)
    assert 9 == len(query_result.universe_set)

    # with unknown that gets excluded
    query_result = query(loaded_database_connection).hasa('reference:5061:G:A').isa(query_result_B)
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 0 == len(query_result.unknown_set)
    assert 8 == len(query_result.absent_set)
    assert 9 == len(query_result.universe_set)


def test_tolist_and_toset(loaded_database_connection: DataIndexConnection):
    conn = loaded_database_connection

    assert {'SampleA', 'SampleB'} == set(query(conn).isin(['SampleA', 'SampleB']).tolist())
    assert {'SampleA', 'SampleB'} == query(conn).isin(['SampleA', 'SampleB']).toset()

    assert {'SampleB'} == set(query(conn).hasa('reference:5061:G:A').tolist())
    assert {'SampleB'} == query(conn).hasa('reference:5061:G:A').toset()

    assert {'SampleA', 'SampleB'} == set(query(conn).hasa('reference:5061:G:A').tolist(include_unknown=True))
    assert {'SampleA', 'SampleB'} == query(conn).hasa('reference:5061:G:A').toset(include_unknown=True)

    assert {'SampleA', 'SampleB', 'SampleC', '2014C-3598', '2014C-3599', '2014D-0067', '2014D-0068',
            'CFSAN002349', 'CFSAN023463'
            } == set(query(conn).hasa('reference:5061:G:A').tolist(include_unknown=True, include_absent=True))
    assert {'SampleA', 'SampleB', 'SampleC', '2014C-3598', '2014C-3599', '2014D-0067', '2014D-0068',
            'CFSAN002349', 'CFSAN023463'
            } == query(conn).hasa('reference:5061:G:A').toset(include_unknown=True, include_absent=True)

    assert {'SampleA', 'SampleC', '2014C-3598', '2014C-3599', '2014D-0067', '2014D-0068',
            'CFSAN002349', 'CFSAN023463'
            } == set(query(conn).hasa('reference:5061:G:A').tolist(include_present=False,
                                                                   include_unknown=True, include_absent=True))
    assert {'SampleA', 'SampleC', '2014C-3598', '2014C-3599', '2014D-0067', '2014D-0068',
            'CFSAN002349', 'CFSAN023463'
            } == query(conn).hasa('reference:5061:G:A').toset(include_present=False,
                                                              include_unknown=True, include_absent=True)

    assert {'SampleC', '2014C-3598', '2014C-3599', '2014D-0067', '2014D-0068',
            'CFSAN002349', 'CFSAN023463'
            } == set(query(conn).hasa('reference:5061:G:A').tolist(include_present=False,
                                                                   include_unknown=False, include_absent=True))
    assert {'SampleC', '2014C-3598', '2014C-3599', '2014D-0067', '2014D-0068',
            'CFSAN002349', 'CFSAN023463'
            } == query(conn).hasa('reference:5061:G:A').toset(include_present=False,
                                                              include_unknown=False, include_absent=True)

    assert {'SampleA'} == set(query(conn).hasa('reference:5061:G:A').tolist(include_present=False,
                                                                            include_unknown=True, include_absent=False))
    assert {'SampleA'} == query(conn).hasa('reference:5061:G:A').toset(include_present=False,
                                                                       include_unknown=True, include_absent=False)


def test_query_isin_kmer(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    query_result = query(loaded_database_connection).isin('SampleA', kind='distance', distance=1.0,
                                                          units='kmer_jaccard')
    assert 3 == len(query_result)
    assert {sampleA.id, sampleB.id, sampleC.id} == set(query_result.sample_set)
    assert 0 == len(query_result.unknown_set)
    assert 9 == len(query_result.universe_set)
    assert "isin_kmer_jaccard('SampleA', dist=1.0, k=31)" == query_result.query_expression()

    # test with unknowns in query
    query_result = query(loaded_database_connection).hasa('reference:5061:G:A').isin(
        'SampleA', kind='distance', distance=1.0, units='kmer_jaccard')
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert {sampleA.id} == set(query_result.unknown_set)
    assert 9 == len(query_result.universe_set)
    assert "reference:5061:G:A AND isin_kmer_jaccard('SampleA', dist=1.0, k=31)" == query_result.query_expression()

    # test with unknowns in query, isin should remove unknowns
    query_result = query(loaded_database_connection).hasa('reference:5061:G:A').isin(
        'SampleB', kind='distance', distance=0, units='kmer_jaccard')
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 0 == len(query_result.unknown_set)
    assert 9 == len(query_result.universe_set)
    assert "reference:5061:G:A AND isin_kmer_jaccard('SampleB', dist=0, k=31)" == query_result.query_expression()


def test_query_isin_kmer_samples_query(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()
    query_result_A = query(loaded_database_connection).isa('SampleA')

    query_result = query(loaded_database_connection).isin(query_result_A, kind='distance', distance=1.0,
                                                          units='kmer_jaccard')
    assert 3 == len(query_result)
    assert {sampleA.id, sampleB.id, sampleC.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)
    assert "isin_kmer_jaccard(<SamplesQueryIndex[selected=11% (1/9) samples, " \
           "unknown=0% (0/9) samples]>, dist=1.0, k=31)" == query_result.query_expression()


def test_query_isin_kmer_samples_set(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()
    sample_set_a = SampleSet([sampleA.id])

    query_result = query(loaded_database_connection).isin(sample_set_a, kind='distance', distance=1.0,
                                                          units='kmer_jaccard')
    assert 3 == len(query_result)
    assert {sampleA.id, sampleB.id, sampleC.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)
    assert "isin_kmer_jaccard(set(1 samples), dist=1.0, k=31)" == query_result.query_expression()


def test_query_isin_kmer_samples_query_no_matches(loaded_database_connection: DataIndexConnection):
    query_result_empty = query(loaded_database_connection).isa('no_exist')

    query_result = query(loaded_database_connection).isin(query_result_empty, kind='distance', distance=1.0,
                                                          units='kmer_jaccard')
    assert 0 == len(query_result)
    assert 9 == len(query_result.universe_set)
    assert "isin_kmer_jaccard(<SamplesQueryIndex[selected=0% (0/9) samples, " \
           "unknown=0% (0/9) samples]>, dist=1.0, k=31)" == query_result.query_expression()


def test_query_within_kmer_default(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    query_result = query(loaded_database_connection).within('SampleA', distance=1.0)
    assert 3 == len(query_result)
    assert {sampleA.id, sampleB.id, sampleC.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)
    assert "isin_kmer_jaccard('SampleA', dist=1.0, k=31)" == query_result.query_expression()


def test_query_within_invalid_unit_with_no_tree(loaded_database_connection: DataIndexConnection):
    with pytest.raises(Exception) as execinfo:
        query(loaded_database_connection).within('SampleA', distance=1.0,
                                                 units='substitutions')
    assert 'units=[substitutions] is not supported' in str(execinfo.value)


def test_to_distances_kmer(loaded_database_connection: DataIndexConnection):
    query_result = query(loaded_database_connection).isin(['SampleA', 'SampleB', 'SampleC'], kind='samples')
    results_d, labels = query_result.to_distances(kind='kmer')

    assert (3, 3) == results_d.shape
    assert {'SampleA', 'SampleB', 'SampleC'} == set(labels)

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


def test_query_isin_kmer_2_matches(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    query_result = query(loaded_database_connection).isin('SampleA', kind='distances', distance=0.5,
                                                          units='kmer_jaccard')
    assert 2 == len(query_result)
    assert {sampleA.id, sampleC.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


def test_query_isin_kmer_1_match(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()

    query_result = query(loaded_database_connection).isin('SampleA', kind='distance', distance=0.49,
                                                          units='kmer_jaccard')
    assert 1 == len(query_result)
    assert {sampleA.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)
    assert "isin_kmer_jaccard('SampleA', dist=0.49, k=31)" == query_result.query_expression()


def test_query_single_mutation_serial_processing(loaded_database_connection: DataIndexConnection):
    do_test_query_single_mutation(loaded_database_connection)


def test_query_single_mutation_parallel_processing(loaded_database_connection_parallel_variants: DataIndexConnection):
    do_test_query_single_mutation(loaded_database_connection_parallel_variants)


def do_test_query_single_mutation(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    all_sample_ids = {i for i, in db.get_session().query(Sample.id).all()}

    query_result = query(loaded_database_connection).hasa(QueryFeatureMutationSPDI('reference:5061:G:A'))
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert {sampleA.id} == set(query_result.unknown_set)
    assert all_sample_ids - {sampleA.id} - {sampleB.id} == set(query_result.absent_set)
    assert 9 == len(query_result.universe_set)
    assert not query_result.has_tree()


def test_query_single_mutation_then_add_new_genomes_and_query(loaded_database_connection: DataIndexConnection,
                                                              snippy_data_package_2: NucleotideSampleDataPackage):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    all_sample_ids = {i for i, in db.get_session().query(Sample.id).all()}
    assert 9 == len(all_sample_ids)

    query_result = query(loaded_database_connection).hasa(QueryFeatureMutationSPDI('reference:5061:G:A'))
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert {sampleA.id} == set(query_result.unknown_set)
    assert all_sample_ids - {sampleA.id} - {sampleB.id} == set(query_result.absent_set)
    assert 9 == len(query_result.universe_set)

    # Insert new data which should increase the number of matches I get
    loaded_database_connection.variation_service.insert(data_package=snippy_data_package_2,
                                                        feature_scope_name='genome')

    sampleA_2 = db.get_session().query(Sample).filter(Sample.name == 'SampleA_2').one()
    sampleB_2 = db.get_session().query(Sample).filter(Sample.name == 'SampleB_2').one()
    all_sample_ids = {i for i, in db.get_session().query(Sample.id).all()}
    assert 12 == len(all_sample_ids)

    query_result = query(loaded_database_connection).hasa(QueryFeatureMutationSPDI('reference:5061:G:A'))
    assert 2 == len(query_result)
    assert {sampleB.id, sampleB_2.id} == set(query_result.sample_set)
    assert {sampleA.id, sampleA_2.id} == set(query_result.unknown_set)
    assert all_sample_ids - {sampleA.id, sampleA_2.id, sampleB.id, sampleB_2.id} == set(query_result.absent_set)
    assert 12 == len(query_result.universe_set)


def test_select_unknown_absent_present(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    all_sample_ids = {i for i, in db.get_session().query(Sample.id).all()}

    query_result = query(loaded_database_connection).hasa('reference:5061:G:A')
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 1 == len(query_result.unknown_set)
    assert {sampleA.id} == set(query_result.unknown_set)
    assert 7 == len(query_result.absent_set)
    assert all_sample_ids - {sampleA.id, sampleB.id} == set(query_result.absent_set)
    assert 9 == len(query_result.universe_set)

    # Select unknown
    selected_result = query_result.select_unknown()
    assert 1 == len(selected_result)
    assert {sampleA.id} == set(selected_result.sample_set)
    assert 0 == len(selected_result.unknown_set)
    assert 8 == len(selected_result.absent_set)
    assert all_sample_ids - {sampleA.id} == set(selected_result.absent_set)
    assert 9 == len(selected_result.universe_set)

    # Select absent
    selected_result = query_result.select_absent()
    assert 7 == len(selected_result)
    assert all_sample_ids - {sampleA.id, sampleB.id} == set(selected_result.sample_set)
    assert 0 == len(selected_result.unknown_set)
    assert 2 == len(selected_result.absent_set)
    assert {sampleA.id, sampleB.id} == set(selected_result.absent_set)
    assert 9 == len(selected_result.universe_set)

    # Select present
    selected_result = query_result.select_present()
    assert 1 == len(selected_result)
    assert {sampleB.id} == set(selected_result.sample_set)
    assert 0 == len(selected_result.unknown_set)
    assert 8 == len(selected_result.absent_set)
    assert all_sample_ids - {sampleB.id} == set(selected_result.absent_set)
    assert 9 == len(selected_result.universe_set)


def test_query_multiple_mutation_unknowns(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()
    all_sample_ids = {i for i, in db.get_session().query(Sample.id).all()}

    # Initial query object
    query_result = query(loaded_database_connection)
    assert 9 == len(query_result)
    assert 0 == len(query_result.unknown_set)
    assert 9 == len(query_result.universe_set)
    assert 0 == len(query_result.absent_set)

    # Test order 'reference:5061:G:A' then 'reference:190:A:G'
    query_result = query(loaded_database_connection).hasa('reference:5061:G:A')
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert {sampleA.id} == set(query_result.unknown_set)
    assert all_sample_ids - {sampleA.id, sampleB.id} == set(query_result.absent_set)
    assert 9 == len(query_result.universe_set)

    query_result = query_result.hasa('reference:190:A:G')
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert {sampleA.id} == set(query_result.unknown_set)
    assert all_sample_ids - {sampleA.id, sampleB.id} == set(query_result.absent_set)
    assert 9 == len(query_result.universe_set)

    # Test order 'reference:190:A:G' then 'reference:5061:G:A'
    query_result = query(loaded_database_connection).hasa('reference:190:A:G')
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert {sampleA.id, sampleC.id} == set(query_result.unknown_set)
    assert all_sample_ids - {sampleA.id, sampleB.id, sampleC.id} == set(query_result.absent_set)
    assert 9 == len(query_result.universe_set)

    query_result = query_result.hasa('reference:5061:G:A')
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert {sampleA.id} == set(query_result.unknown_set)
    assert all_sample_ids - {sampleA.id, sampleB.id} == set(query_result.absent_set)
    assert 9 == len(query_result.universe_set)

    # Reset query
    query_result = query(loaded_database_connection).hasa('reference:190:A:G')
    assert {sampleB.id} == set(query_result.sample_set)
    assert {sampleA.id, sampleC.id} == set(query_result.unknown_set)
    assert all_sample_ids - {sampleA.id, sampleB.id, sampleC.id} == set(query_result.absent_set)

    # Empty query with no unknowns
    query_result_empty = query_result.hasa('reference:800:1:G')
    assert 0 == len(query_result_empty)
    assert set() == set(query_result_empty.sample_set)
    assert 0 == len(query_result_empty.unknown_set)
    assert 9 == len(query_result_empty.universe_set)
    assert 9 == len(query_result_empty.absent_set)

    # Empty query with all unknowns
    query_result_unknown = query_result.hasa('reference:5160:1:G')
    assert 0 == len(query_result_unknown)
    assert {sampleA.id, sampleB.id, sampleC.id} == set(query_result_unknown.unknown_set)
    assert 9 == len(query_result_unknown.universe_set)
    assert all_sample_ids - {sampleA.id, sampleB.id, sampleC.id} == set(query_result.absent_set)

    # Query for A,B then for mutation with B and verify A ends up in unknowns
    query_result = query(loaded_database_connection).isin(['SampleA', 'SampleB'])
    assert 2 == len(query_result)
    assert {sampleA.id, sampleB.id} == set(query_result.sample_set)
    assert 0 == len(query_result.unknown_set)
    assert 9 == len(query_result.universe_set)
    assert all_sample_ids - {sampleA.id, sampleB.id} == set(query_result.absent_set)

    query_result = query_result.hasa('reference:5061:G:A')
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert {sampleA.id} == set(query_result.unknown_set)
    assert 9 == len(query_result.universe_set)
    assert all_sample_ids - {sampleA.id, sampleB.id} == set(query_result.absent_set)


def test_query_intersect_sample_set(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    query_result = query(loaded_database_connection).hasa('reference:5061:G:A')
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert {sampleA.id} == set(query_result.unknown_set)
    assert 9 == len(query_result.universe_set)

    # Test intersect exclude unknown
    query_intersect = query_result.intersect(SampleSet(sample_ids=[sampleB.id]))
    assert 1 == len(query_intersect)
    assert {sampleB.id} == set(query_intersect.sample_set)
    assert 0 == len(query_intersect.unknown_set)
    assert 9 == len(query_intersect.universe_set)

    # Test intersect include unknown
    query_intersect = query_result.intersect(SampleSet(sample_ids=[sampleA.id, sampleB.id]))
    assert 1 == len(query_intersect)
    assert {sampleB.id} == set(query_intersect.sample_set)
    assert 1 == len(query_intersect.unknown_set)
    assert {sampleA.id} == set(query_intersect.unknown_set)
    assert 9 == len(query_intersect.universe_set)

    # Test intersect include unknown and sample not in query
    query_intersect = query_result.intersect(SampleSet(sample_ids=[sampleA.id, sampleB.id, sampleC.id]))
    assert 1 == len(query_intersect)
    assert {sampleB.id} == set(query_intersect.sample_set)
    assert 1 == len(query_intersect.unknown_set)
    assert {sampleA.id} == set(query_intersect.unknown_set)
    assert 9 == len(query_intersect.universe_set)

    # Test intersect only unknown and sample not in query
    query_intersect = query_result.intersect(SampleSet(sample_ids=[sampleA.id, sampleC.id]))
    assert 0 == len(query_intersect)
    assert 1 == len(query_intersect.unknown_set)
    assert {sampleA.id} == set(query_intersect.unknown_set)
    assert 9 == len(query_intersect.universe_set)

    # Test intersect only unknown
    query_intersect = query_result.intersect(SampleSet(sample_ids=[sampleA.id]))
    assert 0 == len(query_intersect)
    assert 1 == len(query_intersect.unknown_set)
    assert {sampleA.id} == set(query_intersect.unknown_set)
    assert 9 == len(query_intersect.universe_set)


def test_query_mutation_hgvs(loaded_database_connection_annotations: DataIndexConnection):
    db = loaded_database_connection_annotations.database
    sample_sh14_001 = db.get_session().query(Sample).filter(Sample.name == 'SH14-001').one()
    sample_sh14_014 = db.get_session().query(Sample).filter(Sample.name == 'SH14-014').one()
    sample_sh10_014 = db.get_session().query(Sample).filter(Sample.name == 'SH10-014').one()

    # hgvs c (nucleotide)
    ## Test using QueryFeature object
    query_result = query(loaded_database_connection_annotations).hasa(
        QueryFeatureHGVS.create_from_id('hgvs:NC_011083:SEHA_RS04550:p.Ile224fs'))
    assert 2 == len(query_result)
    assert {sample_sh14_001.id, sample_sh14_014.id} == set(query_result.sample_set)
    assert 3 == len(query_result.universe_set)
    assert not query_result.has_tree()

    ## Test using string
    query_result = query(loaded_database_connection_annotations).hasa('hgvs:NC_011083:SEHA_RS04550:p.Ile224fs')
    assert 2 == len(query_result)
    assert {sample_sh14_001.id, sample_sh14_014.id} == set(query_result.sample_set)
    assert 3 == len(query_result.universe_set)
    assert not query_result.has_tree()

    # hgvs p (protein)
    ## Test using QueryFeature object
    query_result = query(loaded_database_connection_annotations).hasa(
        QueryFeatureHGVS.create_from_id('hgvs:NC_011083:SEHA_RS04550:c.670dupA'))
    assert 2 == len(query_result)
    assert {sample_sh14_001.id, sample_sh14_014.id} == set(query_result.sample_set)
    assert 3 == len(query_result.universe_set)
    assert not query_result.has_tree()

    ## Test using string
    query_result = query(loaded_database_connection_annotations).hasa('hgvs:NC_011083:SEHA_RS04550:c.670dupA')
    assert 2 == len(query_result)
    assert {sample_sh14_001.id, sample_sh14_014.id} == set(query_result.sample_set)
    assert 3 == len(query_result.universe_set)
    assert not query_result.has_tree()

    # hgvs n (nucleotide)
    ## Test using QueryFeature object
    query_result = query(loaded_database_connection_annotations).hasa(
        QueryFeatureHGVS.create_from_id('hgvs:NC_011083:n.882634G>A'))
    assert 1 == len(query_result)
    assert {sample_sh10_014.id} == set(query_result.sample_set)
    assert 3 == len(query_result.universe_set)
    assert not query_result.has_tree()

    ## Test using string
    query_result = query(loaded_database_connection_annotations).hasa('hgvs:NC_011083:n.882634G>A')
    assert 1 == len(query_result)
    assert {sample_sh10_014.id} == set(query_result.sample_set)
    assert 3 == len(query_result.universe_set)
    assert not query_result.has_tree()

    # Test find no results
    query_result = query(loaded_database_connection_annotations).hasa(
        QueryFeatureHGVS.create_from_id('hgvs:NC_011083:SEHA_RS04550:p.none'))
    assert 0 == len(query_result)


def test_query_single_mutation_complement(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    all_sample_ids = {i for i, in db.get_session().query(Sample.id).all()}

    query_result = query(loaded_database_connection).hasa(QueryFeatureMutationSPDI('reference:5061:G:A'))
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert {sampleA.id} == set(query_result.unknown_set)
    assert all_sample_ids - {sampleA.id, sampleB.id} == set(query_result.absent_set)
    assert 9 == len(query_result.universe_set)

    query_result_c = query_result.complement()

    assert 7 == len(query_result_c)
    assert all_sample_ids - {sampleA.id, sampleB.id} == set(query_result_c.sample_set)
    assert {sampleA.id} == set(query_result_c.unknown_set)
    assert {sampleB.id} == set(query_result_c.absent_set)
    assert 9 == len(query_result_c.universe_set)

    query_result_c = ~query_result

    assert 7 == len(query_result_c)
    assert all_sample_ids - {sampleA.id, sampleB.id} == set(query_result_c.sample_set)
    assert {sampleA.id} == set(query_result_c.unknown_set)
    assert {sampleB.id} == set(query_result_c.absent_set)
    assert 9 == len(query_result_c.universe_set)


def test_query_single_mutation_two_samples_complement(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()
    all_sample_ids = {i for i, in db.get_session().query(Sample.id).all()}

    query_result = query(loaded_database_connection).hasa(QueryFeatureMutationSPDI('reference:839:C:G'))
    assert 2 == len(query_result)
    assert {sampleB.id, sampleC.id} == set(query_result.sample_set)
    assert set() == set(query_result.unknown_set)
    assert all_sample_ids - {sampleB.id, sampleC.id} == set(query_result.absent_set)
    assert 9 == len(query_result.universe_set)

    query_result_c = query_result.complement()

    assert 7 == len(query_result_c)
    assert all_sample_ids - {sampleB.id, sampleC.id} == set(query_result_c.sample_set)
    assert set() == set(query_result_c.unknown_set)
    assert {sampleB.id, sampleC.id} == set(query_result_c.absent_set)
    assert 9 == len(query_result_c.universe_set)

    assert sampleB.id not in query_result_c.sample_set
    assert sampleC.id not in query_result_c.sample_set
    assert sampleA.id in query_result_c.sample_set

    query_result_c = ~query_result

    assert 7 == len(query_result_c)
    assert all_sample_ids - {sampleB.id, sampleC.id} == set(query_result_c.sample_set)
    assert set() == set(query_result_c.unknown_set)
    assert {sampleB.id, sampleC.id} == set(query_result_c.absent_set)
    assert 9 == len(query_result_c.universe_set)

    assert sampleB.id not in query_result_c.sample_set
    assert sampleC.id not in query_result_c.sample_set
    assert sampleA.id in query_result_c.sample_set


def test_query_single_mutation_summary(loaded_database_connection: DataIndexConnection):
    df = query(loaded_database_connection).hasa('reference:5061:G:A', kind='mutation').summary()
    assert 1 == len(df)
    assert ['Query', 'Present', 'Absent', 'Unknown', 'Total',
            '% Present', '% Absent', '% Unknown'] == df.columns.tolist()
    assert 'reference:5061:G:A' == df.iloc[0]['Query']
    assert 1 == df.iloc[0]['Present']
    assert 7 == df.iloc[0]['Absent']
    assert 1 == df.iloc[0]['Unknown']
    assert 9 == df.iloc[0]['Total']
    assert math.isclose((1 / 9) * 100, df.iloc[0]['% Present'])
    assert math.isclose((7 / 9) * 100, df.iloc[0]['% Absent'])
    assert math.isclose((1 / 9) * 100, df.iloc[0]['% Unknown'])


def test_query_single_mutation_two_samples(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    query_result = query(loaded_database_connection).hasa(QueryFeatureMutationSPDI('reference:839:C:G'))
    assert 2 == len(query_result)
    assert {sampleB.id, sampleC.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


def test_query_single_mutation_two_samples_kmer_one_sample(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    query_result = query(loaded_database_connection).hasa(QueryFeatureMutationSPDI('reference:839:C:G'))
    assert 2 == len(query_result)
    assert {sampleB.id, sampleC.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)

    query_result = query_result.isin('SampleA', kind='distance', distance=0.5,
                                     units='kmer_jaccard')
    assert 1 == len(query_result)
    assert {sampleC.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


def test_query_hasa_string_features(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    # Test default hasa SPDI
    query_result = query(loaded_database_connection).hasa('reference:839:C:G')
    assert 2 == len(query_result)
    assert {sampleB.id, sampleC.id} == set(query_result.sample_set)
    assert 0 == len(query_result.unknown_set)
    assert 9 == len(query_result.universe_set)

    # Test HGVS (should return no results since snpeff annotations don't exist)
    query_result = query(loaded_database_connection).hasa('hgvs:reference:n.839C>G')
    assert 0 == len(query_result)
    assert 0 == len(query_result.unknown_set)
    assert 9 == len(query_result.universe_set)


def test_query_hasa_string_features_snpeff(loaded_database_connection_annotations: DataIndexConnection):
    db = loaded_database_connection_annotations.database
    sample_sh14_001 = db.get_session().query(Sample).filter(Sample.name == 'SH14-001').one()
    sample_sh14_014 = db.get_session().query(Sample).filter(Sample.name == 'SH14-014').one()
    sample_sh10_014 = db.get_session().query(Sample).filter(Sample.name == 'SH10-014').one()

    query_result = query(loaded_database_connection_annotations)

    # Test HGVS with amino acid notation
    query_result_test = query_result.hasa('hgvs:NC_011083:SEHA_RS04550:p.Ile224fs')
    assert 2 == len(query_result_test)
    assert {sample_sh14_001.id, sample_sh14_014.id} == set(query_result_test.sample_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 3 == len(query_result_test.universe_set)

    # Test HGVS nucleotide coding notation
    query_result_test = query_result.hasa('hgvs:NC_011083:SEHA_RS04550:c.670dupA')
    assert 2 == len(query_result_test)
    assert {sample_sh14_001.id, sample_sh14_014.id} == set(query_result_test.sample_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 3 == len(query_result_test.universe_set)

    # Test SPDI deletion sequence
    query_result_test = query_result.hasa('NC_011083:835147:C:CA')
    assert 2 == len(query_result_test)
    assert {sample_sh14_001.id, sample_sh14_014.id} == set(query_result_test.sample_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 3 == len(query_result_test.universe_set)

    # Test SPDI deletion integer
    query_result_test = query_result.hasa('NC_011083:835147:1:CA')
    assert 2 == len(query_result_test)
    assert {sample_sh14_001.id, sample_sh14_014.id} == set(query_result_test.sample_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 3 == len(query_result_test.universe_set)

    # Test HGVS, intergenic region
    query_result_test = query_result.hasa('hgvs:NC_011083:n.298943A>T')
    assert 3 == len(query_result_test)
    assert {sample_sh14_001.id, sample_sh14_014.id, sample_sh10_014.id} == set(query_result_test.sample_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 3 == len(query_result_test.universe_set)

    # Test SPDI, sequence and intergenic region
    query_result_test = query_result.hasa('NC_011083:298943:A:T')
    assert 3 == len(query_result_test)
    assert {sample_sh14_001.id, sample_sh14_014.id, sample_sh10_014.id} == set(query_result_test.sample_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 3 == len(query_result_test.universe_set)

    # Test SPDI, deletion integer and intergenic region
    query_result_test = query_result.hasa('NC_011083:298943:1:T')
    assert 3 == len(query_result_test)
    assert {sample_sh14_001.id, sample_sh14_014.id, sample_sh10_014.id} == set(query_result_test.sample_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 3 == len(query_result_test.universe_set)

    # Test HGVS large deletion with amino acid notation
    query_result_test = query_result.hasa('hgvs:NC_011083:SEHA_RS15905:p.Asp140_His155del')
    assert 1 == len(query_result_test)
    assert {sample_sh10_014.id} == set(query_result_test.sample_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 3 == len(query_result_test.universe_set)

    # Test HGVS large deletion with nucleotide coding notation
    query_result_test = query_result.hasa(
        'hgvs:NC_011083:SEHA_RS15905:c.417_464delCGACCACGACCACGACCACGACCACGACCACGACCACGACCACGACCA')
    assert 1 == len(query_result_test)
    assert {sample_sh10_014.id} == set(query_result_test.sample_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 3 == len(query_result_test.universe_set)

    # Test SPDI, large deletion sequence
    query_result_test = query_result.hasa('NC_011083:3167187:AACCACGACCACGACCACGACCACGACCACGACCACGACCACGACCACG:A')
    assert 1 == len(query_result_test)
    assert {sample_sh10_014.id} == set(query_result_test.sample_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 3 == len(query_result_test.universe_set)

    # Test SPDI, large deletion integer
    query_result_test = query_result.hasa(
        f'NC_011083:3167187:{len("AACCACGACCACGACCACGACCACGACCACGACCACGACCACGACCACG")}:A')
    assert 1 == len(query_result_test)
    assert {sample_sh10_014.id} == set(query_result_test.sample_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 3 == len(query_result_test.universe_set)

    # Test HGVS smaller deletion in same region with amino acid notation
    query_result_test = query_result.hasa('hgvs:NC_011083:SEHA_RS15905:p.Asp144_His155del')
    assert 2 == len(query_result_test)
    assert {sample_sh14_001.id, sample_sh14_014.id} == set(query_result_test.sample_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 3 == len(query_result_test.universe_set)

    # Test HGVS smaller deletion in same region with nucleotide coding notation
    query_result_test = query_result.hasa(
        'hgvs:NC_011083:SEHA_RS15905:c.429_464delCGACCACGACCACGACCACGACCACGACCACGACCA')
    assert 2 == len(query_result_test)
    assert {sample_sh14_001.id, sample_sh14_014.id} == set(query_result_test.sample_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 3 == len(query_result_test.universe_set)

    # Test SPDI, smaller deletion in same region sequence
    query_result_test = query_result.hasa('NC_011083:3167187:AACCACGACCACGACCACGACCACGACCACGACCACG:A')
    assert 2 == len(query_result_test)
    assert {sample_sh14_001.id, sample_sh14_014.id} == set(query_result_test.sample_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 3 == len(query_result_test.universe_set)

    # Test SPDI, smaller deletion in same region integer
    query_result_test = query_result.hasa(f'NC_011083:3167187:{len("AACCACGACCACGACCACGACCACGACCACGACCACG")}:A')
    assert 2 == len(query_result_test)
    assert {sample_sh14_001.id, sample_sh14_014.id} == set(query_result_test.sample_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 3 == len(query_result_test.universe_set)

    # Test HGVSGN.c single mutation
    query_result_test = query_result.hasa('hgvs_gn:NC_011083:murF:c.497C>A')
    assert 3 == len(query_result_test)
    assert {sample_sh10_014.id, sample_sh14_001.id, sample_sh14_014.id} == set(query_result_test.sample_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 3 == len(query_result_test.universe_set)

    # Test HGVSGN.p single mutation
    query_result_test = query_result.hasa('hgvs_gn:NC_011083:murF:p.Ala166Glu')
    assert 3 == len(query_result_test)
    assert {sample_sh10_014.id, sample_sh14_001.id, sample_sh14_014.id} == set(query_result_test.sample_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 3 == len(query_result_test.universe_set)

    # Test HGVSGN.c single mutation 2 results
    query_result_test = query_result.hasa('hgvs_gn:NC_011083:oadA:c.609T>C')
    assert 2 == len(query_result_test)
    assert {sample_sh10_014.id, sample_sh14_001.id} == set(query_result_test.sample_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 3 == len(query_result_test.universe_set)

    # Test HGVSGN.p single mutation 2 results
    query_result_test = query_result.hasa('hgvs_gn:NC_011083:oadA:p.Cys203Cys')
    assert 2 == len(query_result_test)
    assert {sample_sh10_014.id, sample_sh14_001.id} == set(query_result_test.sample_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 3 == len(query_result_test.universe_set)

    # Test HGVS.c single mutation 2 results
    query_result_test = query_result.hasa('hgvs:NC_011083:SEHA_RS17780:c.609T>C')
    assert 2 == len(query_result_test)
    assert {sample_sh10_014.id, sample_sh14_001.id} == set(query_result_test.sample_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 3 == len(query_result_test.universe_set)

    # Test HGVS.c single mutation 2 results
    query_result_test = query_result.hasa('hgvs:NC_011083:SEHA_RS17780:p.Cys203Cys')
    assert 2 == len(query_result_test)
    assert {sample_sh10_014.id, sample_sh14_001.id} == set(query_result_test.sample_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 3 == len(query_result_test.universe_set)

    # Test equivalent SPDI identifier for above
    query_result_test = query_result.hasa('NC_011083:3535635:A:G')
    assert 2 == len(query_result_test)
    assert {sample_sh10_014.id, sample_sh14_001.id} == set(query_result_test.sample_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 3 == len(query_result_test.universe_set)


def test_query_hasa_string_features_snpeff_duplicate_genes(
        loaded_database_connection_annotations_duplicate_genes: DataIndexConnection):
    db = loaded_database_connection_annotations_duplicate_genes.database
    sample1 = db.get_session().query(Sample).filter(Sample.name == 'SH10-014-dup-gene-variant').one()
    sample2 = db.get_session().query(Sample).filter(Sample.name == 'SH10-014-dup-gene-variant-2').one()

    query_result = query(loaded_database_connection_annotations_duplicate_genes)

    # Test HGVSGN.c single mutation but there are two different copies of murF gene and so query should investigate both
    query_result_test = query_result.hasa('hgvs_gn:NC_011083:murF:c.497C>A')
    assert 2 == len(query_result_test)
    assert {sample1.id, sample2.id} == set(query_result_test.sample_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 2 == len(query_result_test.universe_set)

    # Test HGVSGN.p single mutation but there are two different copies of murF gene and so query should investigate both
    query_result_test = query_result.hasa('hgvs_gn:NC_011083:murF:p.Ala166Glu')
    assert 2 == len(query_result_test)
    assert {sample1.id, sample2.id} == set(query_result_test.sample_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 2 == len(query_result_test.universe_set)

    # Test HGVS.c single mutation which, because it's selecting by locus identifier which is unique, should give 1 result
    query_result_test = query_result.hasa('hgvs:NC_011083:SEHA_RS01180:c.497C>A')
    assert 1 == len(query_result_test)
    assert {sample1.id} == set(query_result_test.sample_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 2 == len(query_result_test.universe_set)

    # Test HGVS.p single mutation which, because it's selecting by locus identifier which is unique, should give 1 result
    query_result_test = query_result.hasa('hgvs:NC_011083:SEHA_RS01180:p.Ala166Glu')
    assert 1 == len(query_result_test)
    assert {sample1.id} == set(query_result_test.sample_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 2 == len(query_result_test.universe_set)


def test_query_single_mutation_no_results_is_empty(loaded_database_connection: DataIndexConnection):
    # Test is_empty for something with unknown positions
    query_result = query(loaded_database_connection).hasa(QueryFeatureMutationSPDI('reference:1:1:A'))
    assert 0 == len(query_result)
    assert query_result.is_empty()
    assert 9 == len(query_result.universe_set)
    assert not query_result.is_empty(include_unknown=True)

    # Test is_empty for something without unknown positions
    query_result = query(loaded_database_connection).hasa(QueryFeatureMutationSPDI('reference:3000:1:A'))
    assert 0 == len(query_result)
    assert query_result.is_empty()
    assert 9 == len(query_result.universe_set)
    assert query_result.is_empty(include_unknown=True)


def test_query_chained_mutation(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    query_result = query(loaded_database_connection).hasa(
        QueryFeatureMutationSPDI('reference:839:C:G')).hasa(
        QueryFeatureMutationSPDI('reference:5061:G:A'))
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


def test_query_chained_mutation_has_mutation(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    query_result = query(loaded_database_connection).hasa(
        'reference:839:C:G', kind='mutation').hasa(
        'reference:5061:G:A', kind='mutation')
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


def test_query_mlst_allele(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sample_CFSAN002349 = db.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()
    sample_CFSAN023463 = db.get_session().query(Sample).filter(Sample.name == 'CFSAN023463').one()
    sample_2014D_0067 = db.get_session().query(Sample).filter(Sample.name == '2014D-0067').one()
    sample_2014D_0068 = db.get_session().query(Sample).filter(Sample.name == '2014D-0068').one()
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    # No unknowns
    query_result = query(loaded_database_connection).hasa(QueryFeatureMLST('mlst:lmonocytogenes:abcZ:1'))
    assert 5 == len(query_result)
    assert {sample_CFSAN002349.id, sample_CFSAN023463.id, sampleA.id, sampleB.id, sampleC.id} == set(
        query_result.sample_set)
    assert 0 == len(query_result.unknown_set)
    assert 4 == len(query_result.absent_set)
    assert 9 == len(query_result.universe_set)

    assert {'CFSAN002349', 'CFSAN023463', 'SampleA', 'SampleB', 'SampleC'} == set(query_result.tolist())
    assert {sample_CFSAN002349.id, sample_CFSAN023463.id, sampleA.id, sampleB.id, sampleC.id} == set(
        query_result.tolist(names=False))

    # With unknown and present
    query_result = query(loaded_database_connection).hasa(QueryFeatureMLST('mlst:campylobacter:uncA:6'))
    assert 1 == len(query_result)
    assert {sample_2014D_0068.id} == set(query_result.sample_set)
    assert 1 == len(query_result.unknown_set)
    assert {sample_2014D_0067.id} == set(query_result.unknown_set)
    assert 7 == len(query_result.absent_set)
    assert 9 == len(query_result.universe_set)

    # With unknown and absent
    query_result = query(loaded_database_connection).hasa(QueryFeatureMLST('mlst:campylobacter:uncA:5'))
    assert 0 == len(query_result)
    assert 1 == len(query_result.unknown_set)
    assert {sample_2014D_0067.id} == set(query_result.unknown_set)
    assert 8 == len(query_result.absent_set)
    assert 9 == len(query_result.universe_set)

    # Direct from string
    query_result = query(loaded_database_connection).hasa('mlst:campylobacter:uncA:6')
    assert 1 == len(query_result)
    assert {sample_2014D_0068.id} == set(query_result.sample_set)
    assert 1 == len(query_result.unknown_set)
    assert {sample_2014D_0067.id} == set(query_result.unknown_set)
    assert 7 == len(query_result.absent_set)
    assert 9 == len(query_result.universe_set)


def test_query_chained_mlst_alleles(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sample1 = db.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()

    query_result = query(loaded_database_connection).hasa(
        QueryFeatureMLST('mlst:lmonocytogenes:abcZ:1')).hasa(
        QueryFeatureMLST('mlst:lmonocytogenes:lhkA:4'))
    assert 1 == len(query_result)
    assert {sample1.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


def test_query_chained_mlst_alleles_has_allele(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sample1 = db.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()

    query_result = query(loaded_database_connection) \
        .hasa('mlst:lmonocytogenes:abcZ:1', kind='mlst') \
        .hasa('mlst:lmonocytogenes:lhkA:4', kind='mlst')
    assert 1 == len(query_result)
    assert {sample1.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


def test_query_chained_mlst_nucleotide(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()
    all_sample_ids = {i for i, in db.get_session().query(Sample.id).all()}

    # Test query mutation then MLST
    query_result = query(loaded_database_connection) \
        .hasa('reference:839:C:G', kind='mutation') \
        .hasa('mlst:lmonocytogenes:cat:12', kind='mlst')
    assert 1 == len(query_result)
    assert {sampleC.id} == set(query_result.sample_set)
    assert 0 == len(query_result.unknown_set)
    assert 8 == len(query_result.absent_set)
    assert all_sample_ids - {sampleC.id} == set(query_result.absent_set)
    assert 9 == len(query_result.universe_set)

    assert ['SampleC'] == query_result.tolist()
    assert [sampleC.id] == query_result.tolist(names=False)

    # Test query MLST then mutation with a deletion that will be switched to unknown
    query_result = query(loaded_database_connection) \
        .hasa('mlst:lmonocytogenes:cat:11', kind='mlst') \
        .hasa('reference:3897:GCGCA:G', kind='mutation')
    assert 0 == len(query_result)
    assert 1 == len(query_result.unknown_set)
    assert {sampleB.id} == set(query_result.unknown_set)
    assert 8 == len(query_result.absent_set)
    assert all_sample_ids - {sampleB.id} == set(query_result.absent_set)
    assert 9 == len(query_result.universe_set)

    # Test query MLST (with unknown allele) then mutation with a deletion that will be switched to unknown
    query_result = query(loaded_database_connection) \
        .hasa('mlst:lmonocytogenes:ldh:5', kind='mlst') \
        .hasa('reference:3897:GCGCA:G', kind='mutation')
    assert 0 == len(query_result)
    assert 2 == len(query_result.unknown_set)
    assert {sampleB.id, sampleC.id} == set(query_result.unknown_set)
    assert 7 == len(query_result.absent_set)
    assert all_sample_ids - {sampleB.id, sampleC.id} == set(query_result.absent_set)
    assert 9 == len(query_result.universe_set)

    # Test the unknown allele of MLST with a deletion that will be switched to unknown
    query_result = query(loaded_database_connection) \
        .hasa('mlst:lmonocytogenes:ldh:?', kind='mlst') \
        .hasa('reference:3897:GCGCA:G', kind='mutation')
    assert 0 == len(query_result)
    assert 1 == len(query_result.unknown_set)
    assert {sampleB.id} == set(query_result.unknown_set)
    assert 8 == len(query_result.absent_set)
    assert all_sample_ids - {sampleB.id} == set(query_result.absent_set)
    assert 9 == len(query_result.universe_set)

    # Test query MLST (with unknown allele) then mutation (no issues with unknown/found overlap)
    query_result = query(loaded_database_connection) \
        .hasa('mlst:lmonocytogenes:ldh:5', kind='mlst') \
        .hasa('reference:839:C:G', kind='mutation')
    print(query_result.toframe(include_unknown=True)[['Sample Name', 'Status']])
    assert 1 == len(query_result)
    assert {sampleC.id} == set(query_result.sample_set)
    assert 1 == len(query_result.unknown_set)
    assert {sampleB.id} == set(query_result.unknown_set)
    assert 7 == len(query_result.absent_set)
    assert all_sample_ids - {sampleB.id, sampleC.id} == set(query_result.absent_set)
    assert 9 == len(query_result.universe_set)

    # Test query mutation (no issues with unknown/found overlap) then MLST (with unknown allele)
    query_result = query(loaded_database_connection) \
        .hasa('reference:839:C:G', kind='mutation') \
        .hasa('mlst:lmonocytogenes:ldh:5', kind='mlst')
    assert 1 == len(query_result)
    assert {sampleC.id} == set(query_result.sample_set)
    assert 1 == len(query_result.unknown_set)
    assert {sampleB.id} == set(query_result.unknown_set)
    assert 7 == len(query_result.absent_set)
    assert all_sample_ids - {sampleB.id, sampleC.id} == set(query_result.absent_set)
    assert 9 == len(query_result.universe_set)


def test_query_single_mutation_dataframe(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    # Case with no unknowns
    df = query(loaded_database_connection).hasa('reference:839:C:G', kind='mutation').toframe()
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    df = df.sort_values(['Sample Name'])
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert [sampleB.id, sampleC.id] == df['Sample ID'].tolist()
    assert {'reference:839:C:G'} == set(df['Query'].tolist())
    assert ['Present', 'Present'] == df['Status'].tolist()

    # Case with no unknowns, exclude present
    df = query(loaded_database_connection).hasa('reference:839:C:G', kind='mutation').toframe(include_present=False)
    assert 0 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()

    # Case with some unknowns
    df = query(loaded_database_connection).hasa('reference:5061:G:A', kind='mutation').toframe()
    assert 1 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    df = df.sort_values(['Sample Name'])
    assert ['SampleB'] == df['Sample Name'].tolist()
    assert [sampleB.id] == df['Sample ID'].tolist()
    assert {'reference:5061:G:A'} == set(df['Query'].tolist())
    assert ['Present'] == df['Status'].tolist()

    # Case with some unknowns, exclude present
    df = query(loaded_database_connection).hasa('reference:5061:G:A', kind='mutation').toframe(include_present=False)
    assert 0 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()


def test_query_single_mutation_dataframe_include_all(loaded_database_connection: DataIndexConnection):
    # Case of no unknowns
    df = query(loaded_database_connection).hasa(
        'reference:839:C:G', kind='mutation').toframe(include_absent=True, include_unknown=True)
    assert 9 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    df = df.sort_values(['Sample Name'])
    assert ['2014C-3598', '2014C-3599', '2014D-0067', '2014D-0068',
            'CFSAN002349', 'CFSAN023463',
            'SampleA', 'SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert ['Absent', 'Absent', 'Absent', 'Absent',
            'Absent', 'Absent',
            'Absent', 'Present', 'Present'] == df['Status'].tolist()
    assert {'reference:839:C:G'} == set(df['Query'].tolist())

    # Case of no unknowns but 'include_unknown' is False
    df = query(loaded_database_connection).hasa(
        'reference:839:C:G', kind='mutation').toframe(include_absent=True, include_unknown=False)
    assert 9 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    df = df.sort_values(['Sample Name'])
    assert ['2014C-3598', '2014C-3599', '2014D-0067', '2014D-0068',
            'CFSAN002349', 'CFSAN023463',
            'SampleA', 'SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert ['Absent', 'Absent', 'Absent', 'Absent',
            'Absent', 'Absent',
            'Absent', 'Present', 'Present'] == df['Status'].tolist()
    assert {'reference:839:C:G'} == set(df['Query'].tolist())

    # Case of some unknowns
    df = query(loaded_database_connection).hasa(
        'reference:5061:G:A', kind='mutation').toframe(include_absent=True, include_unknown=True)
    assert 9 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    df = df.sort_values(['Sample Name'])
    assert ['2014C-3598', '2014C-3599', '2014D-0067', '2014D-0068',
            'CFSAN002349', 'CFSAN023463',
            'SampleA', 'SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert ['Absent', 'Absent', 'Absent', 'Absent',
            'Absent', 'Absent',
            'Unknown', 'Present', 'Absent'] == df['Status'].tolist()
    assert {'reference:5061:G:A'} == set(df['Query'].tolist())

    # Case of some unknowns excluding unknowns
    df = query(loaded_database_connection).hasa(
        'reference:5061:G:A', kind='mutation').toframe(include_absent=True, include_unknown=False)
    assert 8 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    df = df.sort_values(['Sample Name'])
    assert ['2014C-3598', '2014C-3599', '2014D-0067', '2014D-0068',
            'CFSAN002349', 'CFSAN023463',
            'SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert ['Absent', 'Absent', 'Absent', 'Absent',
            'Absent', 'Absent',
            'Present', 'Absent'] == df['Status'].tolist()
    assert {'reference:5061:G:A'} == set(df['Query'].tolist())

    # Case of some unknowns, excluding absent
    df = query(loaded_database_connection).hasa(
        'reference:5061:G:A', kind='mutation').toframe(include_absent=False, include_unknown=True)
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    df = df.sort_values(['Sample Name'])
    assert ['SampleA', 'SampleB'] == df['Sample Name'].tolist()
    assert ['Unknown', 'Present'] == df['Status'].tolist()
    assert {'reference:5061:G:A'} == set(df['Query'].tolist())

    # Case of some unknowns, only unknowns
    df = query(loaded_database_connection).hasa(
        'reference:5061:G:A', kind='mutation').toframe(include_present=False,
                                                       include_absent=False, include_unknown=True)
    assert 1 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    df = df.sort_values(['Sample Name'])
    assert ['SampleA'] == df['Sample Name'].tolist()
    assert ['Unknown'] == df['Status'].tolist()
    assert {'reference:5061:G:A'} == set(df['Query'].tolist())

    # Case of some unknowns, only absent
    df = query(loaded_database_connection).hasa(
        'reference:5061:G:A', kind='mutation').toframe(include_present=False, include_absent=True,
                                                       include_unknown=False)
    assert 7 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    df = df.sort_values(['Sample Name'])
    assert ['2014C-3598', '2014C-3599', '2014D-0067', '2014D-0068',
            'CFSAN002349', 'CFSAN023463',
            'SampleC'] == df['Sample Name'].tolist()
    assert ['Absent', 'Absent', 'Absent', 'Absent',
            'Absent', 'Absent',
            'Absent'] == df['Status'].tolist()
    assert {'reference:5061:G:A'} == set(df['Query'].tolist())


def test_query_chained_allele_dataframe(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sample1 = db.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()

    df = query(loaded_database_connection) \
        .hasa('mlst:lmonocytogenes:abcZ:1', kind='mlst') \
        .hasa('mlst:lmonocytogenes:lhkA:4', kind='mlst').toframe()

    assert 1 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()

    df = df.sort_values(['Sample Name'])
    assert ['CFSAN002349'] == df['Sample Name'].tolist()
    assert [sample1.id] == df['Sample ID'].tolist()
    assert {'mlst:lmonocytogenes:abcZ:1 AND mlst:lmonocytogenes:lhkA:4'} == set(df['Query'].tolist())


def test_query_single_mutation_no_results_toframe(loaded_database_connection: DataIndexConnection):
    df = query(loaded_database_connection).hasa('reference:1:1:A', kind='mutation').toframe()
    assert 0 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()


def test_query_single_mutation_all_unknown_summary(loaded_database_connection: DataIndexConnection):
    df = query(loaded_database_connection).hasa('reference:1:1:A', kind='mutation').summary()
    assert 1 == len(df)
    assert ['Query', 'Present', 'Absent', 'Unknown', 'Total',
            '% Present', '% Absent', '% Unknown'] == df.columns.tolist()
    assert 'reference:1:1:A' == df.iloc[0]['Query']
    assert 0 == df.iloc[0]['Present']
    assert 6 == df.iloc[0]['Absent']
    assert 3 == df.iloc[0]['Unknown']
    assert 9 == df.iloc[0]['Total']
    assert math.isclose((0 / 9) * 100, df.iloc[0]['% Present'])
    assert math.isclose((6 / 9) * 100, df.iloc[0]['% Absent'])
    assert math.isclose((3 / 9) * 100, df.iloc[0]['% Unknown'])


def test_query_single_mutation_all_absent_summary(loaded_database_connection: DataIndexConnection):
    df = query(loaded_database_connection).hasa('reference:3000:1:A', kind='mutation').summary()
    assert 1 == len(df)
    assert ['Query', 'Present', 'Absent', 'Unknown', 'Total',
            '% Present', '% Absent', '% Unknown'] == df.columns.tolist()
    assert 'reference:3000:1:A' == df.iloc[0]['Query']
    assert 0 == df.iloc[0]['Present']
    assert 9 == df.iloc[0]['Absent']
    assert 0 == df.iloc[0]['Unknown']
    assert 9 == df.iloc[0]['Total']
    assert math.isclose((0 / 9) * 100, df.iloc[0]['% Present'])
    assert math.isclose((9 / 9) * 100, df.iloc[0]['% Absent'])
    assert math.isclose((0 / 9) * 100, df.iloc[0]['% Unknown'])


def test_all_samples_summary(loaded_database_connection: DataIndexConnection):
    df = query(loaded_database_connection).summary()
    assert 1 == len(df)
    assert ['Query', 'Present', 'Absent', 'Unknown', 'Total',
            '% Present', '% Absent', '% Unknown'] == df.columns.tolist()
    assert '' == df.iloc[0]['Query']
    assert 9 == df.iloc[0]['Present']
    assert 0 == df.iloc[0]['Absent']
    assert 0 == df.iloc[0]['Unknown']
    assert 9 == df.iloc[0]['Total']
    assert math.isclose((9 / 9) * 100, df.iloc[0]['% Present'])
    assert math.isclose((0 / 9) * 100, df.iloc[0]['% Absent'])
    assert math.isclose((0 / 9) * 100, df.iloc[0]['% Unknown'])


def test_join_custom_dataframe_no_query(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    df = pd.DataFrame([
        [sampleA.id, 'red'],
        [sampleB.id, 'green'],
        [sampleC.id, 'blue']
    ], columns=['Sample ID', 'Color'])

    query_result = query(loaded_database_connection, universe='dataframe',
                         data_frame=df, sample_ids_column='Sample ID')
    assert 3 == len(query_result)
    assert 3 == len(query_result.universe_set)
    assert {sampleA.id, sampleB.id, sampleC.id} == set(query_result.sample_set)
    assert {'dataframe(ids_col=[Sample ID])'} == set(query_result.toframe()['Query'].tolist())


def test_query_custom_dataframe_isin_samples(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    df = pd.DataFrame([
        [sampleA.id, 'red'],
        [sampleB.id, 'green'],
        [sampleC.id, 'blue']
    ], columns=['Sample ID', 'Color'])

    query_result = query(loaded_database_connection, universe='dataframe',
                         data_frame=df, sample_ids_column='Sample ID')
    query_result = query_result.isin(['SampleA', 'SampleC'])
    assert 2 == len(query_result)
    assert 3 == len(query_result.universe_set)
    assert {sampleA.id, sampleC.id} == set(query_result.sample_set)
    assert {"dataframe(ids_col=[Sample ID]) AND isin_samples(['SampleA', 'SampleC'])"} == set(
        query_result.toframe()['Query'].tolist())


def test_query_custom_dataframe_isin_kmer_distance(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    df = pd.DataFrame([
        [sampleA.id, 'red'],
        [sampleB.id, 'blue'],
        [sampleC.id, 'blue']
    ], columns=['Sample ID', 'Color'])

    main_query = query(loaded_database_connection, universe='dataframe',
                       data_frame=df, sample_ids_column='Sample ID')

    # Test isin with sample name
    query_result = main_query.isin('SampleA', kind='distance', distance=0.5,
                                   units='kmer_jaccard')
    assert 2 == len(query_result)
    assert 3 == len(query_result.universe_set)
    assert {sampleA.id, sampleC.id} == set(query_result.sample_set)

    # Test isin with sample set
    query_result_A = main_query.isa('SampleA', kind='sample')
    query_result = main_query.isin(query_result_A, kind='distance', distance=0.5,
                                   units='kmer_jaccard')
    assert 2 == len(query_result)
    assert 3 == len(query_result.universe_set)
    assert {sampleA.id, sampleC.id} == set(query_result.sample_set)

    # Test query with series expression
    query_result = query_result.isin(df['Color'] == 'blue', kind='dataframe')
    assert 1 == len(query_result)
    assert 3 == len(query_result.universe_set)
    assert {sampleC.id} == set(query_result.sample_set)


def test_query_custom_dataframe_kmer_to_distances(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    df = pd.DataFrame([
        [sampleA.id, 'red'],
        [sampleB.id, 'blue'],
        [sampleC.id, 'blue']
    ], columns=['Sample ID', 'Color'])

    query_result = query(loaded_database_connection, universe='dataframe',
                         data_frame=df, sample_ids_column='Sample ID')
    query_result = query_result.isin('SampleA', kind='distance', distance=0.5,
                                     units='kmer_jaccard')
    assert {sampleA.id, sampleC.id} == set(query_result.sample_set)

    results_d, labels = query_result.to_distances(kind='kmer')

    assert (2, 2) == results_d.shape
    assert {'SampleA', 'SampleC'} == set(labels)

    l = {element: idx for idx, element in enumerate(labels)}

    assert math.isclose(results_d[l['SampleA']][l['SampleA']], 0, rel_tol=1e-3)
    assert math.isclose(results_d[l['SampleA']][l['SampleC']], 0.5, rel_tol=1e-3)

    assert math.isclose(results_d[l['SampleC']][l['SampleA']], 0.5, rel_tol=1e-3)
    assert math.isclose(results_d[l['SampleC']][l['SampleC']], 0, rel_tol=1e-3)


def test_join_custom_dataframe_single_query(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    df = pd.DataFrame([
        [sampleA.id, 'red'],
        [sampleB.id, 'green'],
        [sampleC.id, 'blue']
    ], columns=['Sample ID', 'Color'])

    query_result = query(loaded_database_connection,
                         universe='dataframe',
                         data_frame=df,
                         sample_ids_column='Sample ID')

    assert 3 == len(query_result)
    assert 3 == len(query_result.universe_set)

    query_result = query_result.hasa('reference:839:C:G', kind='mutation')
    assert 2 == len(query_result)
    assert 3 == len(query_result.universe_set)
    assert {sampleB.id, sampleC.id} == set(query_result.sample_set)

    df = query_result.toframe()

    assert 2 == len(df)
    assert {'Query', 'Sample Name', 'Sample ID', 'Status', 'Color'} == set(df.columns.tolist())

    df = df.sort_values(['Sample Name'])
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert [sampleB.id, sampleC.id] == df['Sample ID'].tolist()
    assert ['green', 'blue'] == df['Color'].tolist()
    assert {'dataframe(ids_col=[Sample ID]) AND reference:839:C:G'} == set(df['Query'].tolist())


def test_join_custom_dataframe_query_reset_universe(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    df = pd.DataFrame([
        [sampleA.id, 'red'],
        [sampleB.id, 'green'],
        [sampleC.id, 'blue']
    ], columns=['Sample ID', 'Color'])

    query_result = query(loaded_database_connection,
                         universe='dataframe',
                         data_frame=df,
                         sample_ids_column='Sample ID')

    query_result = query_result.hasa('reference:839:C:G', kind='mutation')
    assert 2 == len(query_result)
    assert 3 == len(query_result.universe_set)
    assert {sampleB.id, sampleC.id} == set(query_result.sample_set)

    query_result = query_result.reset_universe()
    assert 2 == len(query_result)
    assert 2 == len(query_result.universe_set)


def test_join_custom_dataframe_single_query_sample_names_with_unknowns(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    df = pd.DataFrame([
        ['SampleA', 'red'],
        ['SampleB', 'green'],
        ['SampleC', 'blue']
    ], columns=['Samples', 'Color'])

    query_result = query(loaded_database_connection,
                         universe='dataframe',
                         data_frame=df,
                         sample_names_column='Samples')
    assert 3 == len(query_result)
    assert 3 == len(query_result.universe_set)

    # No unknowns
    query_result_test = query_result.hasa('reference:839:C:G', kind='mutation')
    assert 2 == len(query_result_test)
    assert 0 == len(query_result_test.unknown_set)
    assert 3 == len(query_result_test.universe_set)
    assert {sampleB.id, sampleC.id} == set(query_result_test.sample_set)
    df = query_result_test.toframe()
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status', 'Samples', 'Color'] == df.columns.tolist()
    df = df.sort_values(['Sample Name'])
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert [sampleB.id, sampleC.id] == df['Sample ID'].tolist()
    assert ['green', 'blue'] == df['Color'].tolist()
    assert {'dataframe(names_col=[Samples]) AND reference:839:C:G'} == set(df['Query'].tolist())

    # With one unknown and one present
    query_result_test = query_result.hasa('reference:5061:G:A', kind='mutation')
    assert 1 == len(query_result_test)
    assert 1 == len(query_result_test.unknown_set)
    assert 3 == len(query_result_test.universe_set)
    assert {sampleA.id} == set(query_result_test.unknown_set)
    assert {sampleB.id} == set(query_result_test.sample_set)
    df = query_result_test.toframe(include_unknown=True)
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status', 'Samples', 'Color'] == df.columns.tolist()
    df = df.sort_values(['Sample Name'])
    assert ['SampleA', 'SampleB'] == df['Sample Name'].tolist()
    assert ['Unknown', 'Present'] == df['Status'].tolist()
    assert [sampleA.id, sampleB.id] == df['Sample ID'].tolist()
    assert ['red', 'green'] == df['Color'].tolist()
    assert {'dataframe(names_col=[Samples]) AND reference:5061:G:A'} == set(df['Query'].tolist())

    # With all unknown
    query_result_test = query_result.hasa('reference:1:1:T', kind='mutation')
    assert 0 == len(query_result_test)
    assert 3 == len(query_result_test.unknown_set)
    assert 3 == len(query_result_test.universe_set)
    assert {sampleA.id, sampleB.id, sampleC.id} == set(query_result_test.unknown_set)
    df = query_result_test.toframe(include_unknown=True)
    assert 3 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status', 'Samples', 'Color'] == df.columns.tolist()
    df = df.sort_values(['Sample Name'])
    assert ['SampleA', 'SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert ['Unknown', 'Unknown', 'Unknown'] == df['Status'].tolist()
    assert [sampleA.id, sampleB.id, sampleC.id] == df['Sample ID'].tolist()
    assert ['red', 'green', 'blue'] == df['Color'].tolist()
    assert {'dataframe(names_col=[Samples]) AND reference:1:1:T'} == set(df['Query'].tolist())


def test_join_custom_dataframe_extra_sample_names(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    df = pd.DataFrame([
        ['SampleA', 'red'],
        ['SampleB', 'green'],
        ['SampleC', 'blue'],
        ['Extra', 'purple']
    ], columns=['Samples', 'Color'])

    query_result = query(loaded_database_connection,
                         universe='dataframe',
                         data_frame=df,
                         sample_names_column='Samples')
    assert 3 == len(query_result)
    assert 3 == len(query_result.universe_set)

    # No unknowns
    query_result_test = query_result.hasa('reference:839:C:G', kind='mutation')
    assert 2 == len(query_result_test)
    assert 0 == len(query_result_test.unknown_set)
    assert 1 == len(query_result_test.absent_set)
    assert 3 == len(query_result_test.universe_set)
    df = query_result_test.toframe(include_unknown=True)
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status', 'Samples', 'Color'] == df.columns.tolist()
    df = df.sort_values(['Sample Name'])
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert ['Present', 'Present'] == df['Status'].tolist()
    assert [sampleB.id, sampleC.id] == df['Sample ID'].tolist()
    assert ['green', 'blue'] == df['Color'].tolist()
    assert {'dataframe(names_col=[Samples]) AND reference:839:C:G'} == set(df['Query'].tolist())

    # With unknowns
    query_result_test = query_result.hasa('reference:5061:G:A', kind='mutation')
    assert 1 == len(query_result_test)
    assert 1 == len(query_result_test.unknown_set)
    assert 1 == len(query_result_test.absent_set)
    assert 3 == len(query_result_test.universe_set)
    df = query_result_test.toframe(include_unknown=True)
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status', 'Samples', 'Color'] == df.columns.tolist()
    df = df.sort_values(['Sample Name'])
    assert ['SampleA', 'SampleB'] == df['Sample Name'].tolist()
    assert ['Unknown', 'Present'] == df['Status'].tolist()
    assert [sampleA.id, sampleB.id] == df['Sample ID'].tolist()
    assert ['red', 'green'] == df['Color'].tolist()
    assert {'dataframe(names_col=[Samples]) AND reference:5061:G:A'} == set(df['Query'].tolist())


def test_join_custom_dataframe_missing_sample_names(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    df = pd.DataFrame([
        ['SampleA', 'red'],
        ['SampleC', 'blue'],
    ], columns=['Samples', 'Color'])

    query_result = query(loaded_database_connection,
                         universe='dataframe',
                         data_frame=df,
                         sample_names_column='Samples')

    assert 2 == len(df)
    assert 2 == len(query_result.universe_set)

    # No unknowns
    query_result_test = query_result.hasa('reference:839:C:G', kind='mutation')
    assert 2 == len(query_result_test.universe_set)
    assert 0 == len(query_result_test.unknown_set)
    assert 1 == len(query_result_test.absent_set)
    assert 1 == len(query_result_test)
    df = query_result_test.toframe(include_unknown=True).sort_values(['Sample Name'])
    assert ['Query', 'Sample Name', 'Sample ID', 'Status', 'Samples', 'Color'] == df.columns.tolist()
    assert 1 == len(df)
    assert ['SampleC'] == df['Sample Name'].tolist()
    assert ['Present'] == df['Status'].tolist()
    assert [sampleC.id] == df['Sample ID'].tolist()
    assert ['blue'] == df['Color'].tolist()
    assert {'dataframe(names_col=[Samples]) AND reference:839:C:G'} == set(df['Query'].tolist())

    # One unknown, one present but it's missing from dataframe
    query_result_test = query_result.hasa('reference:5061:G:A', kind='mutation')
    assert 2 == len(query_result_test.universe_set)
    assert 1 == len(query_result_test.unknown_set)
    assert 1 == len(query_result_test.absent_set)
    assert 0 == len(query_result_test)
    df = query_result_test.toframe(include_unknown=True).sort_values(['Sample Name'])
    assert ['Query', 'Sample Name', 'Sample ID', 'Status', 'Samples', 'Color'] == df.columns.tolist()
    assert 1 == len(df)
    assert ['SampleA'] == df['Sample Name'].tolist()
    assert ['Unknown'] == df['Status'].tolist()
    assert [sampleA.id] == df['Sample ID'].tolist()
    assert ['red'] == df['Color'].tolist()
    assert {'dataframe(names_col=[Samples]) AND reference:5061:G:A'} == set(df['Query'].tolist())


def test_query_then_join_dataframe_single_query(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    metadata_df = pd.DataFrame([
        [sampleA.id, 'red'],
        [sampleB.id, 'green'],
        [sampleC.id, 'blue']
    ], columns=['Sample ID', 'Color'])

    query_result = query(loaded_database_connection)

    assert 9 == len(query_result)
    assert 9 == len(query_result.universe_set)

    # Case: no unknowns
    query_result_test = query_result.hasa('reference:839:C:G', kind='mutation')
    assert 2 == len(query_result_test)
    assert 0 == len(query_result_test.unknown_set)
    assert 7 == len(query_result_test.absent_set)
    assert 9 == len(query_result_test.universe_set)
    assert {sampleB.id, sampleC.id} == set(query_result_test.sample_set)

    df = query_result_test.toframe()

    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert {'reference:839:C:G'} == set(df['Query'].tolist())

    # Now join data frame
    query_result_test = query_result_test.join(data_frame=metadata_df, sample_ids_column='Sample ID')
    assert 2 == len(query_result_test)
    assert 3 == len(query_result_test.universe_set)

    df = query_result_test.toframe(include_unknown=True)

    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status', 'Color'] == df.columns.tolist()
    assert {'reference:839:C:G AND dataframe(ids_col=[Sample ID])'} == set(df['Query'].tolist())

    df = df.sort_values(['Sample Name'])
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert ['Present', 'Present'] == df['Status'].tolist()
    assert [sampleB.id, sampleC.id] == df['Sample ID'].tolist()
    assert ['green', 'blue'] == df['Color'].tolist()

    # Case: some unknowns
    query_result_test = query_result.hasa('reference:5061:G:A', kind='mutation')
    assert 1 == len(query_result_test)
    assert 1 == len(query_result_test.unknown_set)
    assert 7 == len(query_result_test.absent_set)
    assert 9 == len(query_result_test.universe_set)

    df = query_result_test.toframe(include_unknown=True)

    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert {'reference:5061:G:A'} == set(df['Query'].tolist())

    # Now join data frame
    query_result_test = query_result_test.join(data_frame=metadata_df, sample_ids_column='Sample ID')
    assert 1 == len(query_result_test)
    assert 1 == len(query_result_test.unknown_set)
    assert 1 == len(query_result_test.absent_set)
    assert 3 == len(query_result_test.universe_set)

    df = query_result_test.toframe(include_unknown=True)

    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status', 'Color'] == df.columns.tolist()
    assert {'reference:5061:G:A AND dataframe(ids_col=[Sample ID])'} == set(df['Query'].tolist())

    df = df.sort_values(['Sample Name'])
    assert ['SampleA', 'SampleB'] == df['Sample Name'].tolist()
    assert ['Unknown', 'Present'] == df['Status'].tolist()
    assert [sampleA.id, sampleB.id] == df['Sample ID'].tolist()
    assert ['red', 'green'] == df['Color'].tolist()


def test_query_join_dataframe_isa_dataframe_column(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    metadata_df = pd.DataFrame([
        [sampleA.id, 'red'],
        [sampleB.id, 'green'],
        [sampleC.id, 'blue']
    ], columns=['Sample ID', 'Color'])

    query_result = query(loaded_database_connection).join(data_frame=metadata_df, sample_ids_column='Sample ID')
    assert 3 == len(query_result)
    assert 3 == len(query_result.universe_set)

    # By default, isa should select by 'sample'
    sub_result = query_result.isa('SampleB')
    assert 1 == len(sub_result)
    assert 3 == len(sub_result.universe_set)
    assert {sampleB.id} == set(sub_result.sample_set)

    df = sub_result.toframe()
    assert ['SampleB'] == df['Sample Name'].tolist()
    assert {"dataframe(ids_col=[Sample ID]) AND isa_sample('SampleB')"} == set(df['Query'].tolist())

    # Make sure isa also works when we pass the kind
    sub_result = query_result.isa('SampleB', kind='sample')
    assert 1 == len(sub_result)
    assert 3 == len(sub_result.universe_set)
    assert {sampleB.id} == set(sub_result.sample_set)

    df = sub_result.toframe()
    assert ['SampleB'] == df['Sample Name'].tolist()
    assert {"dataframe(ids_col=[Sample ID]) AND isa_sample('SampleB')"} == set(df['Query'].tolist())

    # If we explicitly pass kind='dataframe' should select by column in dataframe
    sub_result = query_result.isa('red', isa_column='Color', kind='dataframe')
    assert 1 == len(sub_result)
    assert 3 == len(sub_result.universe_set)
    assert {sampleA.id} == set(sub_result.sample_set)

    df = sub_result.toframe()
    assert ['SampleA'] == df['Sample Name'].tolist()
    assert {"dataframe(ids_col=[Sample ID]) AND isa('Color' is 'red')"} == set(df['Query'].tolist())

    # If we pass default_isa_kind when joining, should be able to override defaults
    query_result = query(loaded_database_connection).join(
        data_frame=metadata_df, sample_ids_column='Sample ID', default_isa_kind='dataframe',
        default_isa_column='Color')

    sub_result = query_result.isa('red')
    assert 1 == len(sub_result)
    assert 3 == len(sub_result.universe_set)
    assert {sampleA.id} == set(sub_result.sample_set)

    df = sub_result.toframe()
    assert ['SampleA'] == df['Sample Name'].tolist()
    assert {"dataframe(ids_col=[Sample ID]) AND isa('Color' is 'red')"} == set(df['Query'].tolist())

    # Setting regex=True should let me pass regexes
    sub_result = query_result.isa(r'^re', regex=True)
    assert 1 == len(sub_result)
    assert 3 == len(sub_result.universe_set)
    assert {sampleA.id} == set(sub_result.sample_set)

    df = sub_result.toframe()
    assert ['SampleA'] == df['Sample Name'].tolist()
    assert {"dataframe(ids_col=[Sample ID]) AND isa('Color' contains '^re')"} == set(df['Query'].tolist())

    # Nothing should match this below regex
    sub_result = query_result.isa(r'^ed', regex=True)
    assert 0 == len(sub_result)
    assert 3 == len(sub_result.universe_set)


def test_query_with_unknown_join_dataframe_isa_dataframe_column(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    metadata_df = pd.DataFrame([
        [sampleA.id, 'red', 'small'],
        [sampleB.id, 'red', 'big'],
        [sampleC.id, 'blue', 'big']
    ], columns=['Sample ID', 'Color', 'Size'])

    query_result = query(loaded_database_connection).hasa('reference:5061:G:A').join(data_frame=metadata_df,
                                                                                     sample_ids_column='Sample ID')
    assert 1 == len(query_result)
    assert 1 == len(query_result.unknown_set)
    assert 1 == len(query_result.absent_set)
    assert 3 == len(query_result.universe_set)

    # By default, isa should select by 'sample'
    sub_result = query_result.isa('SampleB')
    assert 1 == len(sub_result)
    assert 0 == len(sub_result.unknown_set)
    assert 2 == len(sub_result.absent_set)
    assert 3 == len(sub_result.universe_set)
    assert {sampleB.id} == set(sub_result.sample_set)

    # isa by sample where sample is unknown
    sub_result = query_result.isa('SampleA')
    assert 0 == len(sub_result)
    assert 1 == len(sub_result.unknown_set)
    assert 2 == len(sub_result.absent_set)
    assert 3 == len(sub_result.universe_set)
    assert {sampleA.id} == set(sub_result.unknown_set)

    # If we explicitly pass kind='dataframe' should select by column in dataframe
    sub_result = query_result.isa('red', isa_column='Color', kind='dataframe')
    assert 1 == len(sub_result)
    assert 1 == len(sub_result.unknown_set)
    assert 1 == len(sub_result.absent_set)
    assert 3 == len(sub_result.universe_set)
    assert {sampleA.id} == set(sub_result.unknown_set)
    assert {sampleB.id} == set(sub_result.sample_set)

    # isa where we exclude unknown
    sub_result = query_result.isa('big', isa_column='Size', kind='dataframe')
    assert 1 == len(sub_result)
    assert 0 == len(sub_result.unknown_set)
    assert 2 == len(sub_result.absent_set)
    assert 3 == len(sub_result.universe_set)
    assert {sampleB.id} == set(sub_result.sample_set)


def test_query_join_dataframe_isin_dataframe_column(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    metadata_df = pd.DataFrame([
        [sampleA.id, 'red'],
        [sampleB.id, 'green'],
        [sampleC.id, 'blue']
    ], columns=['Sample ID', 'Color'])

    query_result = query(loaded_database_connection).join(data_frame=metadata_df, sample_ids_column='Sample ID')
    assert not query_result.has_tree()

    assert 3 == len(query_result)
    assert 3 == len(query_result.universe_set)

    query_result = query_result.isin(metadata_df['Color'] == 'red', kind='dataframe')

    assert 1 == len(query_result)
    assert 3 == len(query_result.universe_set)

    df = query_result.toframe()

    assert ['SampleA'] == df['Sample Name'].tolist()
    assert {'dataframe(ids_col=[Sample ID]) AND isin(subset from series)'} == set(df['Query'].tolist())


def test_query_join_dataframe_isin_dataframe_column_select_two_samples(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    metadata_df = pd.DataFrame([
        [sampleA.id, 'red'],
        [sampleB.id, 'red'],
        [sampleC.id, 'blue']
    ], columns=['Sample ID', 'Color'])

    query_result = query(loaded_database_connection).join(data_frame=metadata_df, sample_ids_column='Sample ID')

    assert 3 == len(query_result)
    assert 3 == len(query_result.universe_set)

    query_result = query_result.isin(metadata_df['Color'] == 'red', kind='dataframe')

    assert 2 == len(query_result)
    assert 3 == len(query_result.universe_set)

    df = query_result.toframe()

    df = df.sort_values(['Sample Name'])
    assert ['SampleA', 'SampleB'] == df['Sample Name'].tolist()
    assert {'dataframe(ids_col=[Sample ID]) AND isin(subset from series)'} == set(df['Query'].tolist())


def test_query_join_dataframe_isin_dataframe_column_invalid_series_index(
        loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    metadata_df = pd.DataFrame([
        [sampleA.id, 'red'],
        [sampleB.id, 'green'],
        [sampleC.id, 'blue']
    ], columns=['Sample ID', 'Color'])

    invalid_series_select = pd.Series([True, True], dtype=bool)

    query_result = query(loaded_database_connection).join(data_frame=metadata_df, sample_ids_column='Sample ID')

    with pytest.raises(Exception) as execinfo:
        query_result.isin(data=invalid_series_select, kind='dataframe')
    assert 'does not have same index as internal data frame' in str(execinfo.value)


def test_query_and_build_mutation_tree(loaded_database_connection: DataIndexConnection):
    query_result = query(loaded_database_connection).hasa('reference:839:C:G', kind='mutation')
    assert 2 == len(query_result)
    assert 9 == len(query_result.universe_set)

    query_result = query_result.build_tree(kind='mutation', scope='genome', include_reference=True)
    assert 2 == len(query_result)
    assert 9 == len(query_result.universe_set)

    assert isinstance(query_result, MutationTreeSamplesQuery)
    query_result = cast(MutationTreeSamplesQuery, query_result)
    assert query_result.reference_included
    assert 'genome' == query_result.reference_name

    assert query_result.tree is not None
    assert {'SampleB', 'SampleC', 'genome'} == set(query_result.tree.get_leaf_names())


def test_build_mutation_tree_include_reference(loaded_database_connection: DataIndexConnection):
    query_result = query(loaded_database_connection)
    assert 9 == len(query_result)
    assert 9 == len(query_result.universe_set)

    query_result = query_result.build_tree(kind='mutation', scope='genome', include_reference=True)
    assert 3 == len(query_result)
    assert 9 == len(query_result.universe_set)

    assert isinstance(query_result, MutationTreeSamplesQuery)
    query_result = cast(MutationTreeSamplesQuery, query_result)
    assert query_result.reference_included
    assert 'genome' == query_result.reference_name

    assert query_result.tree is not None
    assert {'SampleA', 'SampleB', 'SampleC', 'genome'} == set(query_result.tree.get_leaf_names())


def test_build_mutation_tree_no_include_reference(loaded_database_connection: DataIndexConnection):
    query_result = query(loaded_database_connection)
    assert 9 == len(query_result)
    assert 9 == len(query_result.universe_set)

    query_result = query_result.build_tree(kind='mutation', scope='genome', include_reference=False)
    assert 3 == len(query_result)
    assert 9 == len(query_result.universe_set)

    assert isinstance(query_result, MutationTreeSamplesQuery)
    query_result = cast(MutationTreeSamplesQuery, query_result)
    assert not query_result.reference_included
    assert 'genome' == query_result.reference_name

    assert query_result.tree is not None
    assert {'SampleA', 'SampleB', 'SampleC'} == set(query_result.tree.get_leaf_names())


def test_query_build_tree_and_query(loaded_database_connection: DataIndexConnection):
    query_result = query(loaded_database_connection).hasa('reference:839:C:G', kind='mutation')
    assert 2 == len(query_result)
    assert 9 == len(query_result.universe_set)

    query_result = query_result.build_tree(kind='mutation', scope='genome', include_reference=True)
    assert 2 == len(query_result)
    assert 9 == len(query_result.universe_set)

    query_result = query_result.hasa('reference:5061:G:A', kind='mutation')
    assert 1 == len(query_result)
    assert 9 == len(query_result.universe_set)

    # Tree should still be complete
    assert {'SampleB', 'SampleC', 'genome'} == set(query_result.tree.get_leaf_names())


def test_query_build_tree_and_within_kmer(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    query_result = query(loaded_database_connection).hasa('reference:839:C:G', kind='mutation') \
        .build_tree(kind='mutation', scope='genome', include_reference=True)
    assert 2 == len(query_result)
    assert 9 == len(query_result.universe_set)
    assert {sampleB.id, sampleC.id} == set(query_result.sample_set)
    assert query_result.has_tree()

    query_result = query_result.within('SampleA', distance=0.5,
                                       units='kmer_jaccard')
    assert 1 == len(query_result)
    assert 9 == len(query_result.universe_set)
    assert {sampleC.id} == set(query_result.sample_set)


def test_query_build_tree_dataframe(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    df = query(loaded_database_connection).hasa(
        'reference:839:C:G', kind='mutation').build_tree(
        kind='mutation', scope='genome', include_reference=True).hasa(
        'reference:5061:G:A', kind='mutation').toframe()

    assert 1 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()

    assert ['SampleB'] == df['Sample Name'].tolist()
    assert ['Present'] == df['Status'].tolist()
    assert [sampleB.id] == df['Sample ID'].tolist()
    assert ['reference:839:C:G AND mutation_tree(genome) AND reference:5061:G:A'] == df['Query'].tolist()


def test_query_then_build_tree_then_join_dataframe(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    metadata_df = pd.DataFrame([
        [sampleA.id, 'red'],
        [sampleB.id, 'green'],
        [sampleC.id, 'blue']
    ], columns=['Sample ID', 'Color'])

    query_result = query(loaded_database_connection).hasa('reference:839:C:G', kind='mutation')
    assert 2 == len(query_result)
    assert 9 == len(query_result.universe_set)

    query_result = query_result.build_tree(kind='mutation', scope='genome', include_reference=True)
    assert 2 == len(query_result)
    assert 9 == len(query_result.universe_set)
    assert isinstance(query_result, TreeSamplesQuery)

    df = query_result.toframe()

    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert {'reference:839:C:G AND mutation_tree(genome)'} == set(df['Query'].tolist())

    # Now join data frame
    query_result_join = query_result.join(data_frame=metadata_df, sample_ids_column='Sample ID')
    assert 2 == len(query_result_join)
    assert 3 == len(query_result_join.universe_set)
    assert isinstance(query_result_join, TreeSamplesQuery)

    df = query_result_join.toframe()

    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status', 'Color'] == df.columns.tolist()
    assert {'reference:839:C:G AND mutation_tree(genome) AND dataframe(ids_col=[Sample ID])'} == set(
        df['Query'].tolist())

    df = df.sort_values(['Sample Name'])
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert [sampleB.id, sampleC.id] == df['Sample ID'].tolist()
    assert ['green', 'blue'] == df['Color'].tolist()

    # I should still be able to perform within queries since I have a tree attached
    query_result = query_result_join.isin(['SampleB', 'SampleC'], kind='mrca')
    assert 2 == len(query_result)
    assert 3 == len(query_result.universe_set)
    assert {'SampleB', 'SampleC'} == set(query_result.tolist())

    # mrca from samples query
    query_result_BC = query_result_join.isin(['SampleB', 'SampleC'], kind='samples')
    query_result = query_result_join.isin(query_result_BC, kind='mrca')
    assert 2 == len(query_result)
    assert 3 == len(query_result.universe_set)
    assert {'SampleB', 'SampleC'} == set(query_result.tolist())

    # mrca from samples set
    sample_set_BC = SampleSet([sampleB.id, sampleC.id])
    query_result = query_result_join.isin(sample_set_BC, kind='mrca')
    assert 2 == len(query_result)
    assert 3 == len(query_result.universe_set)
    assert {'SampleB', 'SampleC'} == set(query_result.tolist())

    # mrca from samples query, empty result
    query_result_empty = query_result_join.isin([], kind='samples')
    query_result = query_result_join.isin(query_result_empty, kind='mrca')
    assert 0 == len(query_result)
    assert 3 == len(query_result.universe_set)

    # Resetting universe should work properly
    query_result_BC = query_result_join.isin(['SampleB', 'SampleC'], kind='samples')
    query_result = query_result_join.isin(query_result_BC, kind='mrca')
    assert 2 == len(query_result)
    assert 3 == len(query_result.universe_set)
    query_result = query_result.reset_universe()
    assert 2 == len(query_result)
    assert 2 == len(query_result.universe_set)


def test_empty_universe_tree_query(prebuilt_tree: Tree, loaded_database_connection: DataIndexConnection):
    query_result = query(loaded_database_connection).join_tree(
        tree=prebuilt_tree, kind='mutation',
        reference_name='genome',
        alignment_length=5180)

    query_result = query_result.isa('invalid_name')

    assert 0 == len(query_result)
    assert 0 == len(query_result.unknown_set)
    assert 9 == len(query_result.absent_set)
    assert 9 == len(query_result.universe_set)

    query_result = query_result.reset_universe()
    assert 0 == len(query_result)
    assert 0 == len(query_result.unknown_set)
    assert 0 == len(query_result.absent_set)
    assert 0 == len(query_result.universe_set)

    assert '<MutationTreeSamplesQuery[' in str(query_result)
    assert 'selected=NA%' in str(query_result)
    assert 'unknown=NA%' in str(query_result)


def test_query_tree_join_dataframe_isa_dataframe_column(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    metadata_df = pd.DataFrame([
        [sampleA.id, 'red'],
        [sampleB.id, 'green'],
        [sampleC.id, pd.NA]
    ], columns=['Sample ID', 'Color'])

    query_result = query(loaded_database_connection).hasa('reference:839:C:G', kind='mutation')
    assert 2 == len(query_result)
    assert {sampleB.id, sampleC.id} == set(query_result.sample_set)
    assert not query_result.has_tree()

    query_result = query_result.build_tree(kind='mutation', scope='genome', include_reference=True)
    assert 2 == len(query_result)
    assert isinstance(query_result, TreeSamplesQuery)
    assert query_result.has_tree()

    query_result = query_result.join(data_frame=metadata_df, sample_ids_column='Sample ID')
    assert 2 == len(query_result)
    assert query_result.has_tree()

    # Now try to do isa on tree + dataframe query
    query_result = query_result.isa('green', kind='dataframe', isa_column='Color')
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)

    # Do isa with regex and an NA in the dataframe
    query_result = query_result.isa(r'^gr', kind='dataframe', isa_column='Color', regex=True)
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)


def test_query_join_tree_join_dataframe(prebuilt_tree: Tree, loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    metadata_df = pd.DataFrame([
        [sampleA.id, 'red'],
        [sampleB.id, 'green'],
        [sampleC.id, 'blue']
    ], columns=['Sample ID', 'Color'])

    query_result = query(loaded_database_connection).join_tree(tree=prebuilt_tree, kind='mutation',
                                                               reference_name='genome',
                                                               alignment_length=5180)
    assert 3 == len(query_result)
    assert 9 == len(query_result.universe_set)
    assert {sampleA.id, sampleB.id, sampleC.id} == set(query_result.sample_set)
    assert query_result.has_tree()

    query_result = query_result.join(data_frame=metadata_df, sample_ids_column='Sample ID')
    assert 3 == len(query_result)
    assert query_result.has_tree()

    # Now try to do isa on tree + dataframe query
    query_result_color = query_result.isa('green', kind='dataframe', isa_column='Color')
    assert 1 == len(query_result_color)
    assert {sampleB.id} == set(query_result_color.sample_set)

    # mrca A and B
    df = query_result.isin(['SampleA', 'SampleC'], kind='mrca').toframe().sort_values('Sample Name')
    assert 3 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status', 'Color'] == df.columns.tolist()
    assert ['SampleA', 'SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert ['red', 'green', 'blue'] == df['Color'].tolist()


def test_query_tree_join_dataframe_hasa_snpeff(loaded_database_connection_annotations: DataIndexConnection):
    db = loaded_database_connection_annotations.database
    sample_sh14_001 = db.get_session().query(Sample).filter(Sample.name == 'SH14-001').one()
    sample_sh14_014 = db.get_session().query(Sample).filter(Sample.name == 'SH14-014').one()
    sample_sh10_014 = db.get_session().query(Sample).filter(Sample.name == 'SH10-014').one()

    metadata_df = pd.DataFrame([
        [sample_sh14_001.id, 'red'],
        [sample_sh14_014.id, 'red'],
        [sample_sh10_014.id, 'blue']
    ], columns=['Sample ID', 'Color'])

    tree = Tree(str(snpeff_tree_file))

    # Query hasa with tree query
    query_result = query(loaded_database_connection_annotations).join_tree(tree,
                                                                           kind='mutation',
                                                                           alignment_length=4888768,
                                                                           reference_name='NC_011083')
    assert 3 == len(query_result)
    assert isinstance(query_result, TreeSamplesQuery)
    assert query_result.has_tree()
    query_result = query_result.hasa('hgvs:NC_011083:SEHA_RS04550:p.Ile224fs')
    assert 2 == len(query_result)
    assert {sample_sh14_001.id, sample_sh14_014.id} == set(query_result.sample_set)

    # Query hasa with database query
    query_result = query(loaded_database_connection_annotations).join(data_frame=metadata_df,
                                                                      sample_ids_column='Sample ID')
    assert 3 == len(query_result)
    assert isinstance(query_result, DataFrameSamplesQuery)
    assert not query_result.has_tree()
    query_result = query_result.hasa('hgvs:NC_011083:SEHA_RS04550:p.Ile224fs')
    assert 2 == len(query_result)
    assert {sample_sh14_001.id, sample_sh14_014.id} == set(query_result.sample_set)

    # Query hasa with tree joined to dataframe
    query_result = query(loaded_database_connection_annotations).join_tree(tree,
                                                                           kind='mutation',
                                                                           alignment_length=4888768,
                                                                           reference_name='NC_011083')
    query_result = query_result.join(data_frame=metadata_df, sample_ids_column='Sample ID')
    assert 3 == len(query_result)
    assert isinstance(query_result, TreeSamplesQuery)
    assert query_result.has_tree()
    query_result = query_result.hasa('hgvs:NC_011083:SEHA_RS04550:p.Ile224fs')
    assert 2 == len(query_result)
    assert {sample_sh14_001.id, sample_sh14_014.id} == set(query_result.sample_set)

    # Query hasa with experimental tree query
    query_result = query(loaded_database_connection_annotations).join_tree(tree,
                                                                           kind='mutation',
                                                                           alignment_length=4888768,
                                                                           reference_name='NC_011083')
    query_result = cast(MutationTreeSamplesQuery, query_result)
    query_result = ExperimentalTreeSamplesQuery.from_tree_query(query_result)
    assert 3 == len(query_result)
    assert isinstance(query_result, ExperimentalTreeSamplesQuery)
    assert query_result.has_tree()
    query_result = query_result.hasa('hgvs:NC_011083:SEHA_RS04550:p.Ile224fs')
    assert 2 == len(query_result)
    assert {sample_sh14_001.id, sample_sh14_014.id} == set(query_result.sample_set)


def test_within_constructed_tree(loaded_database_connection: DataIndexConnection):
    query_result = query(loaded_database_connection).hasa(
        'reference:839:C:G', kind='mutation').build_tree(
        kind='mutation', scope='genome', include_reference=True, extra_params='--seed 42 -m GTR')
    assert 2 == len(query_result)
    assert 9 == len(query_result.universe_set)

    # subs/site
    df = query_result.isin('SampleC', kind='distance', distance=0.005,
                           units='substitutions/site').toframe().sort_values('Sample Name')
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {"reference:839:C:G AND mutation_tree(genome) AND within(0.005 substitutions/site of 'SampleC')"
            } == set(df['Query'].tolist())

    # subs/site using within
    df = query_result.within('SampleC', distance=0.005, units='substitutions/site').toframe().sort_values('Sample Name')
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {"reference:839:C:G AND mutation_tree(genome) AND within(0.005 substitutions/site of 'SampleC')"
            } == set(df['Query'].tolist())

    # subs
    df = query_result.isin('SampleC', kind='distance', distance=26, units='substitutions').toframe().sort_values(
        'Sample Name')
    assert 2 == len(df)
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {"reference:839:C:G AND mutation_tree(genome) AND within(26 substitutions of 'SampleC')"
            } == set(df['Query'].tolist())

    # subs using samples query
    query_result_C = query_result.isin(['SampleC'], kind='samples')
    df = query_result.isin(query_result_C, kind='distance', distance=26, units='substitutions').toframe().sort_values(
        'Sample Name')
    assert 2 == len(df)
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {'reference:839:C:G AND mutation_tree(genome) AND within(26 substitutions of '
            '<MutationTreeSamplesQuery[selected=11% (1/9) samples, unknown=0% (0/9) samples]>)'
            } == set(df['Query'].tolist())

    # subs using samples query, empty result
    query_result_empty = query_result.isin([], kind='samples')
    df = query_result.isin(query_result_empty, kind='distance', distance=26,
                           units='substitutions').toframe().sort_values(
        'Sample Name')
    assert 0 == len(df)

    # should not include reference genome
    df = query_result.isin('SampleC', kind='distance', distance=100, units='substitutions').toframe().sort_values(
        'Sample Name')
    assert 2 == len(df)
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {"reference:839:C:G AND mutation_tree(genome) AND within(100 substitutions of 'SampleC')"
            } == set(df['Query'].tolist())

    # should have only query sample
    df = query_result.isin('SampleC', kind='distance', distance=1, units='substitutions').toframe().sort_values(
        'Sample Name')
    assert 1 == len(df)
    assert ['SampleC'] == df['Sample Name'].tolist()
    assert {"reference:839:C:G AND mutation_tree(genome) AND within(1 substitutions of 'SampleC')"
            } == set(df['Query'].tolist())

    # should have only query sample, using samples query as input
    query_result_C = query_result.isin(['SampleC'], kind='samples')
    df = query_result.isin(query_result_C, kind='distance', distance=1, units='substitutions').toframe().sort_values(
        'Sample Name')
    assert 1 == len(df)
    assert ['SampleC'] == df['Sample Name'].tolist()
    assert {"reference:839:C:G AND mutation_tree(genome) AND within(1 substitutions of "
            "<MutationTreeSamplesQuery[selected=11% (1/9) samples, unknown=0% (0/9) samples]>)"
            } == set(df['Query'].tolist())

    # mrca
    df = query_result.isin(['SampleB', 'SampleC'], kind='mrca').toframe().sort_values('Sample Name')
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {"reference:839:C:G AND mutation_tree(genome) AND within(mrca of ['SampleB', 'SampleC'])"
            } == set(df['Query'].tolist())

    # mrca of samples query
    query_result_BC = query_result.isin(['SampleB', 'SampleC'], kind='samples')
    df = query_result.isin(query_result_BC, kind='mrca').toframe().sort_values('Sample Name')
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {"reference:839:C:G AND mutation_tree(genome) AND within(mrca of "
            "<MutationTreeSamplesQuery[selected=22% (2/9) samples, unknown=0% (0/9) samples]>)"
            } == set(df['Query'].tolist())

    # mrca of samples query, single sample as result
    query_result_B = query_result.isin(['SampleB'], kind='samples')
    df = query_result.isin(query_result_B, kind='mrca').toframe().sort_values('Sample Name')
    assert 1 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleB'] == df['Sample Name'].tolist()
    assert {"reference:839:C:G AND mutation_tree(genome) AND within(mrca of "
            "<MutationTreeSamplesQuery[selected=11% (1/9) samples, unknown=0% (0/9) samples]>)"
            } == set(df['Query'].tolist())

    # mrca of samples query, empty_results
    query_result_empty = query_result.isin([], kind='samples')
    df = query_result.isin(query_result_empty, kind='mrca').toframe().sort_values('Sample Name')
    assert 0 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()

    # Samples isin()
    df = query_result.isin(['SampleA', 'SampleC']).toframe().sort_values('Sample Name')
    assert 1 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleC'] == df['Sample Name'].tolist()
    assert {"reference:839:C:G AND mutation_tree(genome) AND isin_samples(['SampleA', 'SampleC'])"
            } == set(df['Query'].tolist())

    # Samples isin from samples query()
    query_result_AC = query_result.isin(['SampleA', 'SampleC'], kind='samples')
    df = query_result.isin(query_result_AC).toframe().sort_values('Sample Name')
    assert 1 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleC'] == df['Sample Name'].tolist()
    assert {
               "reference:839:C:G AND mutation_tree(genome) AND "
               "isin_samples(<MutationTreeSamplesQuery[selected=11% (1/9) samples, unknown=0% (0/9) samples]>)"
           } == set(df['Query'].tolist())

    # Sample isa()
    df = query_result.isa('SampleC').toframe().sort_values('Sample Name')
    assert 1 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleC'] == df['Sample Name'].tolist()
    assert {"reference:839:C:G AND mutation_tree(genome) AND isa_sample('SampleC')"
            } == set(df['Query'].tolist())

    # Sample isa() empty
    df = query_result.isa('SampleA').toframe().sort_values('Sample Name')
    assert 0 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()

    # kmer jaccard still works
    df = query_result.within('SampleA', distance=0.5, units='kmer_jaccard').toframe().sort_values('Sample Name')
    assert 1 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleC'] == df['Sample Name'].tolist()
    assert {"reference:839:C:G AND mutation_tree(genome) AND isin_kmer_jaccard('SampleA', dist=0.5, k=31)"
            } == set(df['Query'].tolist())


def test_build_tree_experimental(loaded_database_connection: DataIndexConnection):
    query_result = query(loaded_database_connection).hasa(
        'reference:839:C:G', kind='mutation').build_tree(
        kind='mutation_experimental', scope='genome', include_reference=True, extra_params='--seed 42 -m GTR')
    assert 2 == len(query_result)
    assert isinstance(query_result, ExperimentalTreeSamplesQuery)

    # isin should still work with ExperimentalTreeSamplesQuery
    df = query_result.isin('SampleC', kind='distance', distance=0.005,
                           units='substitutions/site').toframe().sort_values('Sample Name')
    assert 2 == len(df)


def test_within_constructed_tree_larger_tree(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    # Construct new tree with all the samples
    query_result = query(loaded_database_connection).build_tree(
        kind='mutation', scope='genome', include_reference=True, extra_params='--seed 42 -m GTR')
    assert 3 == len(query_result)

    # mrca B and C
    df = query_result.isin(['SampleB', 'SampleC'], kind='mrca').toframe().sort_values('Sample Name')
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {"mutation_tree(genome) AND within(mrca of ['SampleB', 'SampleC'])"
            } == set(df['Query'].tolist())

    # mrca of B and C, samples query
    query_result_BC = query_result.isin(['SampleB', 'SampleC'], kind='samples')
    df = query_result.isin(query_result_BC, kind='mrca').toframe().sort_values('Sample Name')
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {"mutation_tree(genome) AND within(mrca of "
            "<MutationTreeSamplesQuery[selected=22% (2/9) samples, unknown=0% (0/9) samples]>)"
            } == set(df['Query'].tolist())

    # mrca A and B
    df = query_result.isin(['SampleA', 'SampleC'], kind='mrca').toframe().sort_values('Sample Name')
    assert 3 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleA', 'SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {"mutation_tree(genome) AND within(mrca of ['SampleA', 'SampleC'])"
            } == set(df['Query'].tolist())

    # mrca of A and B, samples query
    query_result_AB = query_result.isin(['SampleA', 'SampleB'], kind='samples')
    df = query_result.isin(query_result_AB, kind='mrca').toframe().sort_values('Sample Name')
    assert 3 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleA', 'SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {"mutation_tree(genome) AND within(mrca of "
            "<MutationTreeSamplesQuery[selected=22% (2/9) samples, unknown=0% (0/9) samples]>)"
            } == set(df['Query'].tolist())

    # mrca of A and B, samples set
    sample_set_AB = SampleSet([sampleA.id, sampleB.id])
    df = query_result.isin(sample_set_AB, kind='mrca').toframe().sort_values('Sample Name')
    assert 3 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleA', 'SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {"mutation_tree(genome) AND within(mrca of "
            "set(2 samples))"
            } == set(df['Query'].tolist())

    # hasa putting all in unknown and then mrca of A and B, samples query
    query_result_AB = query_result.isin(['SampleA', 'SampleB'], kind='samples')
    df = query_result.hasa('reference:1:1:T').isin(
        query_result_AB, kind='mrca').toframe(include_unknown=True).sort_values('Sample Name')
    assert 3 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleA', 'SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert ['Unknown', 'Unknown', 'Unknown'] == df['Status'].tolist()
    assert {"mutation_tree(genome) AND reference:1:1:T AND within(mrca of "
            "<MutationTreeSamplesQuery[selected=22% (2/9) samples, unknown=0% (0/9) samples]>)"
            } == set(df['Query'].tolist())

    # hasa putting all in unknown and then mrca of A and B, by name
    df = query_result.hasa('reference:1:1:T').isin(
        ['SampleA', 'SampleB'], kind='mrca').toframe(include_unknown=True).sort_values('Sample Name')
    assert 3 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleA', 'SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert ['Unknown', 'Unknown', 'Unknown'] == df['Status'].tolist()
    assert {"mutation_tree(genome) AND reference:1:1:T AND within(mrca of ['SampleA', 'SampleB'])"
            } == set(df['Query'].tolist())

    # hasa mrca of A and B (by name) and then putting all in unknown
    df = query_result.isin(
        ['SampleA', 'SampleB'], kind='mrca').hasa(
        'reference:1:1:T').toframe(include_unknown=True).sort_values('Sample Name')
    assert 3 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleA', 'SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert ['Unknown', 'Unknown', 'Unknown'] == df['Status'].tolist()
    assert {"mutation_tree(genome) AND within(mrca of ['SampleA', 'SampleB']) AND reference:1:1:T"
            } == set(df['Query'].tolist())

    # hasa only A in unknown and then mrca of A and B, samples query
    query_result_AB = query_result.isin(['SampleA', 'SampleB'], kind='samples')
    df = query_result.hasa('reference:5061:G:A').isin(
        query_result_AB, kind='mrca').toframe(include_unknown=True).sort_values('Sample Name')
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleA', 'SampleB'] == df['Sample Name'].tolist()
    assert ['Unknown', 'Present'] == df['Status'].tolist()
    assert {"mutation_tree(genome) AND reference:5061:G:A AND within(mrca of "
            "<MutationTreeSamplesQuery[selected=22% (2/9) samples, unknown=0% (0/9) samples]>)"
            } == set(df['Query'].tolist())

    # hasa only A in unknown and then mrca of A and B, by name
    df = query_result.hasa('reference:5061:G:A').isin(
        ['SampleA', 'SampleB'], kind='mrca').toframe(include_unknown=True).sort_values('Sample Name')
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleA', 'SampleB'] == df['Sample Name'].tolist()
    assert ['Unknown', 'Present'] == df['Status'].tolist()
    assert {"mutation_tree(genome) AND reference:5061:G:A AND within(mrca of ['SampleA', 'SampleB'])"
            } == set(df['Query'].tolist())

    # mrca of A and B (by name) and hasa putting only A in unknown (and matching B)
    df = query_result.isin(
        ['SampleA', 'SampleB'], kind='mrca').hasa(
        'reference:5061:G:A').toframe(include_unknown=True).sort_values('Sample Name')
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleA', 'SampleB'] == df['Sample Name'].tolist()
    assert ['Unknown', 'Present'] == df['Status'].tolist()
    assert {"mutation_tree(genome) AND within(mrca of ['SampleA', 'SampleB']) AND reference:5061:G:A"
            } == set(df['Query'].tolist())

    # hasa only A in unknown and then mrca of A, by name
    df = query_result.hasa('reference:5061:G:A').isin(
        'SampleA', kind='mrca').toframe(include_unknown=True).sort_values('Sample Name')
    assert 1 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleA'] == df['Sample Name'].tolist()
    assert ['Unknown'] == df['Status'].tolist()
    assert {"mutation_tree(genome) AND reference:5061:G:A AND within(mrca of ['SampleA'])"
            } == set(df['Query'].tolist())


def test_within_joined_mutations_tree(prebuilt_tree: Tree, loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    query_result = query(loaded_database_connection).join_tree(tree=prebuilt_tree, kind='mutation',
                                                               reference_name='genome',
                                                               alignment_length=5180)
    assert 3 == len(query_result)

    # mrca B and C
    df = query_result.isin(['SampleB', 'SampleC'], kind='mrca').toframe().sort_values('Sample Name')
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {"join_tree(4 leaves) AND within(mrca of ['SampleB', 'SampleC'])"
            } == set(df['Query'].tolist())

    # mrca of B and C, samples query
    query_result_BC = query_result.isin(['SampleB', 'SampleC'], kind='samples')
    df = query_result.isin(query_result_BC, kind='mrca').toframe().sort_values('Sample Name')
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {"join_tree(4 leaves) AND within(mrca of "
            "<MutationTreeSamplesQuery[selected=22% (2/9) samples, unknown=0% (0/9) samples]>)"
            } == set(df['Query'].tolist())

    # mrca A and B
    df = query_result.isin(['SampleA', 'SampleC'], kind='mrca').toframe().sort_values('Sample Name')
    assert 3 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleA', 'SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {"join_tree(4 leaves) AND within(mrca of ['SampleA', 'SampleC'])"
            } == set(df['Query'].tolist())

    # mrca of A and B, samples query
    query_result_AB = query_result.isin(['SampleA', 'SampleB'], kind='samples')
    df = query_result.isin(query_result_AB, kind='mrca').toframe().sort_values('Sample Name')
    assert 3 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleA', 'SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {"join_tree(4 leaves) AND within(mrca of "
            "<MutationTreeSamplesQuery[selected=22% (2/9) samples, unknown=0% (0/9) samples]>)"
            } == set(df['Query'].tolist())

    # mrca of A and B, samples set
    sample_set_AB = SampleSet([sampleA.id, sampleB.id])
    df = query_result.isin(sample_set_AB, kind='mrca').toframe().sort_values('Sample Name')
    assert 3 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleA', 'SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {"join_tree(4 leaves) AND within(mrca of "
            "set(2 samples))"
            } == set(df['Query'].tolist())

    # subs/site using samples query
    query_result_C = query_result.isin(['SampleC'], kind='samples')
    df = query_result.isin(query_result_C, kind='distance', distance=2,
                           units='substitutions/site').toframe().sort_values(
        'Sample Name')
    assert 2 == len(df)
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {'join_tree(4 leaves) AND within(2 substitutions/site of '
            '<MutationTreeSamplesQuery[selected=11% (1/9) samples, unknown=0% (0/9) samples]>)'
            } == set(df['Query'].tolist())

    # hasa to put A, B, and C in unknown and then subs/site using samples query to only select B and C
    query_result_C = query_result.isin(['SampleC'], kind='samples')
    df = query_result.hasa('reference:1:1:T').isin(query_result_C, kind='distance', distance=2,
                                                   units='substitutions/site').toframe(
        include_unknown=True).sort_values('Sample Name')
    assert 2 == len(df)
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert ['Unknown', 'Unknown'] == df['Status'].tolist()
    assert {'join_tree(4 leaves) AND reference:1:1:T AND within(2 substitutions/site of '
            '<MutationTreeSamplesQuery[selected=11% (1/9) samples, unknown=0% (0/9) samples]>)'
            } == set(df['Query'].tolist())

    # hasa to put A, B, and C in unknown and then subs/site using samples query to select A, B, and C
    query_result_C = query_result.isin(['SampleC'], kind='samples')
    df = query_result.hasa('reference:1:1:T').isin(query_result_C, kind='distance', distance=4,
                                                   units='substitutions/site').toframe(
        include_unknown=True).sort_values('Sample Name')
    assert 3 == len(df)
    assert ['SampleA', 'SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert ['Unknown', 'Unknown', 'Unknown'] == df['Status'].tolist()
    assert {'join_tree(4 leaves) AND reference:1:1:T AND within(4 substitutions/site of '
            '<MutationTreeSamplesQuery[selected=11% (1/9) samples, unknown=0% (0/9) samples]>)'
            } == set(df['Query'].tolist())

    # hasa to put A, B, and C in unknown and then subs/site using sample names to only select B and C
    df = query_result.hasa('reference:1:1:T').isin(['SampleC'], kind='distance', distance=2,
                                                   units='substitutions/site').toframe(
        include_unknown=True).sort_values('Sample Name')
    assert 2 == len(df)
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert ['Unknown', 'Unknown'] == df['Status'].tolist()
    assert {"join_tree(4 leaves) AND reference:1:1:T AND within(2 substitutions/site of ['SampleC'])"
            } == set(df['Query'].tolist())

    # hasa to put A, B, and C in unknown and then subs/site using sample names to select A, B, and C
    df = query_result.hasa('reference:1:1:T').isin(['SampleA', 'SampleC'], kind='distance', distance=2,
                                                   units='substitutions/site').toframe(
        include_unknown=True).sort_values('Sample Name')
    assert 3 == len(df)
    assert ['SampleA', 'SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert ['Unknown', 'Unknown', 'Unknown'] == df['Status'].tolist()
    assert {"join_tree(4 leaves) AND reference:1:1:T AND within(2 substitutions/site of ['SampleA', 'SampleC'])"
            } == set(df['Query'].tolist())

    # hasa to put A in unknown, select B and then subs/site using sample names to select A, B, and C
    df = query_result.hasa('reference:5061:G:A').isin(['SampleA', 'SampleC'], kind='distance', distance=2,
                                                      units='substitutions/site').toframe(
        include_unknown=True).sort_values('Sample Name')
    assert 2 == len(df)
    assert ['SampleA', 'SampleB'] == df['Sample Name'].tolist()
    assert ['Unknown', 'Present'] == df['Status'].tolist()
    assert {"join_tree(4 leaves) AND reference:5061:G:A AND within(2 substitutions/site of ['SampleA', 'SampleC'])"
            } == set(df['Query'].tolist())

    # subs/site using sample names to select A, B, and C, and then hasa to put A in unknown, select B
    df = query_result.isin(['SampleA', 'SampleC'], kind='distance', distance=2, units='substitutions/site').hasa(
        'reference:5061:G:A').toframe(include_unknown=True).sort_values('Sample Name')
    assert 2 == len(df)
    assert ['SampleA', 'SampleB'] == df['Sample Name'].tolist()
    assert ['Unknown', 'Present'] == df['Status'].tolist()
    assert {"join_tree(4 leaves) AND within(2 substitutions/site of ['SampleA', 'SampleC']) AND reference:5061:G:A"
            } == set(df['Query'].tolist())

    # subs samples query
    df = query_result.isin(['SampleC'], kind='distance', distance=2 * 5180,
                           units='substitutions').toframe().sort_values('Sample Name')
    assert 2 == len(df)
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {f'join_tree(4 leaves) AND within({2 * 5180} substitutions of '
            "['SampleC'])"
            } == set(df['Query'].tolist())

    # kmer jaccard still works
    df = query_result.within('SampleA', distance=0.5, units='kmer_jaccard').toframe().sort_values('Sample Name')
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleA', 'SampleC'] == df['Sample Name'].tolist()
    assert {"join_tree(4 leaves) AND isin_kmer_jaccard('SampleA', dist=0.5, k=31)"
            } == set(df['Query'].tolist())


def test_join_kmer_tree(prebuilt_tree: Tree, loaded_database_connection: DataIndexConnection):
    query_result = query(loaded_database_connection).join_tree(tree=prebuilt_tree, kind='kmer')
    assert 3 == len(query_result)
    assert 9 == len(query_result.universe_set)
    assert isinstance(query_result.tree, ClusterTree)
    assert {'SampleA', 'SampleB', 'SampleC'} == set(query_result.tree.get_leaf_names())


def test_summary_features_kindmutations(loaded_database_connection: DataIndexConnection):
    dfA = pd.read_csv(snippy_all_dataframes['SampleA'], sep='\t')
    dfB = pd.read_csv(snippy_all_dataframes['SampleB'], sep='\t')
    dfC = pd.read_csv(snippy_all_dataframes['SampleC'], sep='\t')
    expected_df = pd.concat([dfA, dfB, dfC])
    expected_df = expected_df.groupby('Mutation').agg({
        'Sequence': 'first',
        'Position': 'first',
        'Deletion': 'first',
        'Insertion': 'first',
        'Mutation': 'count',
    }).rename(columns={'Mutation': 'Count'}).sort_index()
    expected_df['Total'] = 9
    expected_df['Percent'] = 100 * (expected_df['Count'] / expected_df['Total'])

    f = ['reference:461:AAAT:G']
    warnings.warn(f'Removing {f} from expected until I can figure out how to properly handle these')
    expected_df = expected_df.drop(f)

    mutations_df = query(loaded_database_connection).features_summary(ignore_annotations=True)
    mutations_df = mutations_df.sort_index()

    assert len(expected_df) == len(mutations_df)
    assert list(expected_df.columns) == list(mutations_df.columns)
    assert list(expected_df.index) == list(mutations_df.index)
    assert list(expected_df['Count']) == list(mutations_df['Count'])
    assert list(expected_df['Total']) == list(mutations_df['Total'])
    assert math.isclose(100 * (2 / 9), mutations_df.loc['reference:619:G:C', 'Percent'])

    # Test including unknowns
    mutations_df = query(loaded_database_connection).features_summary(ignore_annotations=True,
                                                                      include_unknown_features=True)
    mutations_df['Percent'] = mutations_df['Percent'].astype(int)  # Convert to int for easier comparison
    mutations_df = mutations_df.sort_index()

    assert 632 == len(mutations_df)
    assert list(expected_df.columns) == list(mutations_df.columns)
    assert 2 == mutations_df.loc['reference:619:G:C', 'Count']
    assert 2 == mutations_df.loc['reference:3063:A:ATGCAGC', 'Count']
    assert 1 == mutations_df.loc['reference:1984:GTGATTG:TTGA', 'Count']
    assert 1 == mutations_df.loc['reference:866:GCCAGATCC:G', 'Count']
    assert 3 == mutations_df.loc['reference:90:T:?', 'Count']
    assert 2 == mutations_df.loc['reference:190:A:?', 'Count']
    assert 1 == mutations_df.loc['reference:887:T:?', 'Count']

    # Test only include unknowns
    mutations_df = query(loaded_database_connection).features_summary(ignore_annotations=True,
                                                                      include_unknown_features=True,
                                                                      include_present_features=False)
    mutations_df['Percent'] = mutations_df['Percent'].astype(int)  # Convert to int for easier comparison
    mutations_df = mutations_df.sort_index()

    assert 521 == len(mutations_df)
    assert list(expected_df.columns) == list(mutations_df.columns)
    assert 3 == mutations_df.loc['reference:90:T:?', 'Count']
    assert 2 == mutations_df.loc['reference:190:A:?', 'Count']
    assert 1 == mutations_df.loc['reference:887:T:?', 'Count']
    assert 'reference:619:G:C' not in mutations_df


def test_summary_features_kindmutations_unique(loaded_database_connection: DataIndexConnection):
    dfA = pd.read_csv(snippy_all_dataframes['SampleA'], sep='\t')
    dfB = pd.read_csv(snippy_all_dataframes['SampleB'], sep='\t')
    dfC = pd.read_csv(snippy_all_dataframes['SampleC'], sep='\t')

    # Unique to A
    expected_df = dfA
    expected_df = expected_df.groupby('Mutation').agg({
        'Sequence': 'first',
        'Position': 'first',
        'Deletion': 'first',
        'Insertion': 'first',
        'Mutation': 'count',
    }).rename(columns={'Mutation': 'Count'}).sort_index()
    expected_df['Total'] = 1
    expected_df['Percent'] = 100 * (expected_df['Count'] / expected_df['Total'])

    f = ['reference:461:AAAT:G']
    warnings.warn(f'Removing {f} from expected until I can figure out how to properly handle these')
    expected_df = expected_df.drop(f)

    q = query(loaded_database_connection)

    mutations_df = q.isa('SampleA').features_summary(selection='unique', ignore_annotations=True)
    mutations_df = mutations_df.sort_index()

    assert len(expected_df) == len(mutations_df)
    assert 45 == len(mutations_df)  # Check length against independently generated length
    assert list(expected_df.index) == list(mutations_df.index)
    assert list(expected_df['Count']) == list(mutations_df['Count'])
    assert list(expected_df['Total']) == list(mutations_df['Total'])
    assert math.isclose(100 * (1 / 1), mutations_df.loc['reference:3656:CATT:C', 'Percent'])

    # Unique to B
    dfAC = pd.concat([dfA, dfC])
    expected_df = dfB[~dfB['Mutation'].isin(list(dfAC['Mutation']))]
    expected_df = expected_df.groupby('Mutation').agg({
        'Sequence': 'first',
        'Position': 'first',
        'Deletion': 'first',
        'Insertion': 'first',
        'Mutation': 'count',
    }).rename(columns={'Mutation': 'Count'}).sort_index()
    expected_df['Total'] = 1
    expected_df['Percent'] = 100 * (expected_df['Count'] / expected_df['Total'])

    mutations_df = q.isa('SampleB').features_summary(selection='unique', ignore_annotations=True)
    mutations_df = mutations_df.sort_index()

    assert len(expected_df) == len(mutations_df)
    assert list(expected_df.index) == list(mutations_df.index)
    assert list(expected_df['Count']) == list(mutations_df['Count'])
    assert list(expected_df['Total']) == list(mutations_df['Total'])
    assert math.isclose(100 * (1 / 1), mutations_df.loc['reference:349:AAGT:A', 'Percent'])

    # Unique to C
    dfAB = pd.concat([dfA, dfB])
    expected_df = dfC[~dfC['Mutation'].isin(list(dfAB['Mutation']))]
    expected_df = expected_df.groupby('Mutation').agg({
        'Sequence': 'first',
        'Position': 'first',
        'Deletion': 'first',
        'Insertion': 'first',
        'Mutation': 'count',
    }).rename(columns={'Mutation': 'Count'}).sort_index()
    expected_df['Total'] = 1
    expected_df['Percent'] = 100 * (expected_df['Count'] / expected_df['Total'])

    mutations_df = q.isa('SampleC').features_summary(selection='unique', ignore_annotations=True)
    mutations_df = mutations_df.sort_index()

    assert len(expected_df) == len(mutations_df)
    assert list(expected_df.index) == list(mutations_df.index)
    assert list(expected_df['Count']) == list(mutations_df['Count'])
    assert list(expected_df['Total']) == list(mutations_df['Total'])
    assert math.isclose(100 * (1 / 1), mutations_df.loc['reference:866:GCCAGATCC:G', 'Percent'])

    # Unique to BC
    dfBC = pd.concat([dfB, dfC])
    expected_df = dfBC[~dfBC['Mutation'].isin(list(dfA['Mutation']))]
    expected_df = expected_df.groupby('Mutation').agg({
        'Sequence': 'first',
        'Position': 'first',
        'Deletion': 'first',
        'Insertion': 'first',
        'Mutation': 'count',
    }).rename(columns={'Mutation': 'Count'}).sort_index()
    expected_df['Total'] = 2
    expected_df['Percent'] = 100 * (expected_df['Count'] / expected_df['Total'])

    mutations_df = q.isin(['SampleB', 'SampleC']).features_summary(selection='unique', ignore_annotations=True)
    mutations_df = mutations_df.sort_index()

    assert len(expected_df) == len(mutations_df)
    assert 66 == len(mutations_df)  # Check length against independently generated length
    assert list(expected_df.index) == list(mutations_df.index)
    assert list(expected_df['Count']) == list(mutations_df['Count'])
    assert list(expected_df['Total']) == list(mutations_df['Total'])
    assert math.isclose(100 * (2 / 2), mutations_df.loc['reference:619:G:C', 'Percent'])
    assert math.isclose(100 * (1 / 2), mutations_df.loc['reference:866:GCCAGATCC:G', 'Percent'])
    assert math.isclose(100 * (1 / 2), mutations_df.loc['reference:349:AAGT:A', 'Percent'])

    # Unique to ABC (all with mutations)
    dfABC = pd.concat([dfA, dfB, dfC])
    expected_df = dfABC
    expected_df = expected_df.groupby('Mutation').agg({
        'Sequence': 'first',
        'Position': 'first',
        'Deletion': 'first',
        'Insertion': 'first',
        'Mutation': 'count',
    }).rename(columns={'Mutation': 'Count'}).sort_index()
    expected_df['Total'] = 3
    expected_df['Percent'] = 100 * (expected_df['Count'] / expected_df['Total'])

    f = ['reference:461:AAAT:G']
    warnings.warn(f'Removing {f} from expected until I can figure out how to properly handle these')
    expected_df = expected_df.drop(f)

    mutations_df = q.isin(['SampleA', 'SampleB', 'SampleC']).features_summary(selection='unique',
                                                                              ignore_annotations=True)
    mutations_df = mutations_df.sort_index()

    assert len(expected_df) == len(mutations_df)
    assert 111 == len(mutations_df)  # Check length against independently generated length
    assert list(expected_df.index) == list(mutations_df.index)
    assert list(expected_df['Count']) == list(mutations_df['Count'])
    assert list(expected_df['Total']) == list(mutations_df['Total'])
    assert math.isclose(100 * (2 / 3), mutations_df.loc['reference:619:G:C', 'Percent'])
    assert math.isclose(100 * (1 / 3), mutations_df.loc['reference:866:GCCAGATCC:G', 'Percent'])
    assert math.isclose(100 * (1 / 3), mutations_df.loc['reference:349:AAGT:A', 'Percent'])
    assert math.isclose(100 * (1 / 3), mutations_df.loc['reference:3656:CATT:C', 'Percent'])

    # Unique to None
    mutations_df = q.isin([]).features_summary(selection='unique', ignore_annotations=True)
    mutations_df = mutations_df.sort_index()

    assert 0 == len(mutations_df)
    assert ['Sequence', 'Position', 'Deletion', 'Insertion',
            'Count', 'Total', 'Percent'] == list(mutations_df.columns)


def test_summary_features_kindmutations_annotations(loaded_database_connection_annotations: DataIndexConnection):
    q = query(loaded_database_connection_annotations)

    # 1 sample
    mutations_df = q.isa('SH10-014').features_summary(ignore_annotations=False)

    assert ['Sequence', 'Position', 'Deletion', 'Insertion',
            'Count', 'Total', 'Percent', 'Annotation', 'Annotation_Impact',
            'Gene_Name', 'Gene_ID', 'Feature_Type', 'Transcript_BioType',
            'HGVS.c', 'HGVS.p', 'ID_HGVS.c', 'ID_HGVS.p', 'ID_HGVS_GN.c', 'ID_HGVS_GN.p'] == list(mutations_df.columns)
    assert 139 == len(mutations_df)

    ## Convert percent to int to make it easier to compare in assert statements
    mutations_df['Percent'] = mutations_df['Percent'].astype(int)

    ## missense variant
    assert ['NC_011083', 140658, 'C', 'A', 1, 1, 100,
            'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
            'c.497C>A', 'p.Ala166Glu',
            'hgvs:NC_011083:SEHA_RS01180:c.497C>A', 'hgvs:NC_011083:SEHA_RS01180:p.Ala166Glu',
            'hgvs_gn:NC_011083:murF:c.497C>A', 'hgvs_gn:NC_011083:murF:p.Ala166Glu'] == list(
        mutations_df.loc['NC_011083:140658:C:A'])

    ## inframe deletion
    assert ['NC_011083', 4465400, 'GGCCGAA', 'G', 1, 1, 100,
            'conservative_inframe_deletion', 'MODERATE', 'tyrB', 'SEHA_RS22180', 'transcript', 'protein_coding',
            'c.157_162delGAAGCC', 'p.Glu53_Ala54del',
            'hgvs:NC_011083:SEHA_RS22180:c.157_162delGAAGCC', 'hgvs:NC_011083:SEHA_RS22180:p.Glu53_Ala54del',
            'hgvs_gn:NC_011083:tyrB:c.157_162delGAAGCC', 'hgvs_gn:NC_011083:tyrB:p.Glu53_Ala54del'] == list(
        mutations_df.loc['NC_011083:4465400:GGCCGAA:G'])

    ## Intergenic variant (with some NA values in fields)
    assert ['NC_011083', 4555461, 'T', 'TC', 1, 1, 100,
            'intergenic_region', 'MODIFIER', 'SEHA_RS22510-SEHA_RS26685', 'SEHA_RS22510-SEHA_RS26685',
            'intergenic_region', 'NA',
            'n.4555461_4555462insC', 'NA',
            'hgvs:NC_011083:n.4555461_4555462insC', 'NA',
            'hgvs_gn:NC_011083:n.4555461_4555462insC', 'NA'] == list(
        mutations_df.loc['NC_011083:4555461:T:TC'].fillna('NA'))

    # 3 samples
    mutations_df = q.isin(['SH10-014', 'SH14-001', 'SH14-014']).features_summary(ignore_annotations=False)

    ## Convert percent to int to make it easier to compare in assert statements
    mutations_df['Percent'] = mutations_df['Percent'].astype(int)

    assert ['Sequence', 'Position', 'Deletion', 'Insertion',
            'Count', 'Total', 'Percent', 'Annotation', 'Annotation_Impact',
            'Gene_Name', 'Gene_ID', 'Feature_Type', 'Transcript_BioType',
            'HGVS.c', 'HGVS.p', 'ID_HGVS.c', 'ID_HGVS.p', 'ID_HGVS_GN.c', 'ID_HGVS_GN.p'] == list(mutations_df.columns)
    assert 177 == len(mutations_df)

    ## missense variant (3/3)
    assert ['NC_011083', 140658, 'C', 'A', 3, 3, 100,
            'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
            'c.497C>A', 'p.Ala166Glu',
            'hgvs:NC_011083:SEHA_RS01180:c.497C>A', 'hgvs:NC_011083:SEHA_RS01180:p.Ala166Glu',
            'hgvs_gn:NC_011083:murF:c.497C>A', 'hgvs_gn:NC_011083:murF:p.Ala166Glu'] == list(
        mutations_df.loc['NC_011083:140658:C:A'])

    ## Intergenic variant (1/3)
    assert ['NC_011083', 4555461, 'T', 'TC', 1, 3, 33,
            'intergenic_region', 'MODIFIER', 'SEHA_RS22510-SEHA_RS26685', 'SEHA_RS22510-SEHA_RS26685',
            'intergenic_region', 'NA',
            'n.4555461_4555462insC', 'NA',
            'hgvs:NC_011083:n.4555461_4555462insC', 'NA',
            'hgvs_gn:NC_011083:n.4555461_4555462insC', 'NA'] == list(
        mutations_df.loc['NC_011083:4555461:T:TC'].fillna('NA'))

    # Test ignore annotations
    mutations_df = q.isin(['SH10-014', 'SH14-001', 'SH14-014']).features_summary(ignore_annotations=True)

    assert ['Sequence', 'Position', 'Deletion', 'Insertion',
            'Count', 'Total', 'Percent'] == list(mutations_df.columns)
    assert 177 == len(mutations_df)

    ## Test unique
    mutations_df = q.isa('SH10-014').features_summary(selection='unique', ignore_annotations=False)

    ## Convert percent to int to make it easier to compare in assert statements
    mutations_df['Percent'] = mutations_df['Percent'].astype(int)

    assert ['Sequence', 'Position', 'Deletion', 'Insertion',
            'Count', 'Total', 'Percent', 'Annotation', 'Annotation_Impact',
            'Gene_Name', 'Gene_ID', 'Feature_Type', 'Transcript_BioType',
            'HGVS.c', 'HGVS.p', 'ID_HGVS.c', 'ID_HGVS.p', 'ID_HGVS_GN.c', 'ID_HGVS_GN.p'] == list(mutations_df.columns)
    assert 60 == len(mutations_df)

    ## missense variant
    assert ['NC_011083', 2049576, 'A', 'C', 1, 1, 100,
            'missense_variant', 'MODERATE', 'cutC', 'SEHA_RS10675', 'transcript', 'protein_coding',
            'c.536T>G', 'p.Val179Gly',
            'hgvs:NC_011083:SEHA_RS10675:c.536T>G', 'hgvs:NC_011083:SEHA_RS10675:p.Val179Gly',
            'hgvs_gn:NC_011083:cutC:c.536T>G', 'hgvs_gn:NC_011083:cutC:p.Val179Gly'] == list(
        mutations_df.loc['NC_011083:2049576:A:C'])


def test_summary_features_two(loaded_database_connection: DataIndexConnection):
    dfB = pd.read_csv(snippy_all_dataframes['SampleB'], sep='\t')
    dfC = pd.read_csv(snippy_all_dataframes['SampleC'], sep='\t')
    expected_df = pd.concat([dfB, dfC])
    expected_df = expected_df.groupby('Mutation').agg({
        'Sequence': 'first',
        'Position': 'first',
        'Deletion': 'first',
        'Insertion': 'first',
        'Mutation': 'count',
    }).rename(columns={'Mutation': 'Count'}).sort_index()
    expected_df['Total'] = 2
    expected_df['Percent'] = 100 * (expected_df['Count'] / expected_df['Total'])

    mutations_df = query(loaded_database_connection).hasa(
        'reference:839:C:G', kind='mutation').features_summary()
    mutations_df = mutations_df.sort_index()

    assert len(expected_df) == len(mutations_df)
    assert list(expected_df.index) == list(mutations_df.index)
    assert list(expected_df['Count']) == list(mutations_df['Count'])
    assert list(expected_df['Total']) == list(mutations_df['Total'])
    assert math.isclose(100 * (2 / 2), mutations_df.loc['reference:839:C:G', 'Percent'])


def test_summary_features_kindmlst(loaded_database_connection: DataIndexConnection):
    # Test case of summary of single sample
    summary_df = query(loaded_database_connection).isa('SampleA').features_summary(kind='mlst')
    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    assert 7 == len(summary_df)
    assert {'lmonocytogenes'} == set(summary_df['Scheme'].tolist())
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Total', 'Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'abcZ', '1', 1, 1, 100] == summary_df.loc['mlst:lmonocytogenes:abcZ:1'].tolist()
    assert ['lmonocytogenes', 'bglA', '51', 1, 1, 100] == summary_df.loc['mlst:lmonocytogenes:bglA:51'].tolist()

    # Test samples across multiple schemes
    summary_df = query(loaded_database_connection).isin(['SampleA', 'SampleB', 'CFSAN002349',
                                                         '2014D-0067']).features_summary(kind='mlst')
    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    assert 15 == len(summary_df)
    assert {'lmonocytogenes', 'campylobacter'} == set(summary_df['Scheme'].tolist())
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Total', 'Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'abcZ', '1', 3, 4, 75] == summary_df.loc['mlst:lmonocytogenes:abcZ:1'].tolist()
    assert ['lmonocytogenes', 'bglA', '51', 2, 4, 50] == summary_df.loc['mlst:lmonocytogenes:bglA:51'].tolist()
    assert ['lmonocytogenes', 'bglA', '52', 1, 4, 25] == summary_df.loc['mlst:lmonocytogenes:bglA:52'].tolist()
    assert ['lmonocytogenes', 'ldh', '5', 2, 4, 50] == summary_df.loc['mlst:lmonocytogenes:ldh:5'].tolist()
    assert ['lmonocytogenes', 'lhkA', '5', 2, 4, 50] == summary_df.loc['mlst:lmonocytogenes:lhkA:5'].tolist()
    assert ['lmonocytogenes', 'lhkA', '4', 1, 4, 25] == summary_df.loc['mlst:lmonocytogenes:lhkA:4'].tolist()
    assert ['campylobacter', 'aspA', '2', 1, 4, 25] == summary_df.loc['mlst:campylobacter:aspA:2'].tolist()
    assert ['campylobacter', 'glyA', '3', 1, 4, 25] == summary_df.loc['mlst:campylobacter:glyA:3'].tolist()
    assert 6 == len(summary_df[summary_df['Scheme'] == 'campylobacter'])  # Missing one feature since it's unknown

    # Test only unknown
    summary_df = query(loaded_database_connection).isin(
        ['SampleA', 'SampleB', 'CFSAN002349', '2014D-0067']).features_summary(
        kind='mlst', include_present_features=False, include_unknown_features=True)
    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    assert 2 == len(summary_df)
    assert {'lmonocytogenes', 'campylobacter'} == set(summary_df['Scheme'].tolist())
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Total', 'Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'ldh', '?', 1, 4, 25] == summary_df.loc['mlst:lmonocytogenes:ldh:?'].tolist()
    assert ['campylobacter', 'uncA', '?', 1, 4, 25] == summary_df.loc['mlst:campylobacter:uncA:?'].tolist()

    # Test only unknown, restrict scheme
    summary_df = query(loaded_database_connection).isin(
        ['SampleA', 'SampleB', 'CFSAN002349', '2014D-0067']).features_summary(
        kind='mlst', scheme='lmonocytogenes', include_present_features=False, include_unknown_features=True)
    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    assert 1 == len(summary_df)
    assert {'lmonocytogenes'} == set(summary_df['Scheme'].tolist())
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Total', 'Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'ldh', '?', 1, 4, 25] == summary_df.loc['mlst:lmonocytogenes:ldh:?'].tolist()


def test_tofeaturesset_all(loaded_database_only_snippy: DataIndexConnection):
    dfA = pd.read_csv(snippy_all_dataframes['SampleA'], sep='\t')
    dfB = pd.read_csv(snippy_all_dataframes['SampleB'], sep='\t')
    dfC = pd.read_csv(snippy_all_dataframes['SampleC'], sep='\t')
    expected_df = pd.concat([dfA, dfB, dfC])
    expected_df = expected_df.groupby('Mutation').agg({
        'Sequence': 'first',
        'Position': 'first',
        'Deletion': 'first',
        'Insertion': 'first',
        'Mutation': 'count',
    }).rename(columns={'Mutation': 'Count'}).sort_index()

    f = ['reference:461:AAAT:G']
    warnings.warn(f'Removing {f} from expected until I can figure out how to properly handle these')
    expected_df = expected_df.drop(f)

    expected_set = set(expected_df.index)

    # Test default (no unknown features)
    mutations = query(loaded_database_only_snippy).tofeaturesset()
    assert 111 == len(mutations)
    assert set(expected_set) == set(mutations)

    # Test include unknowns
    mutations = query(loaded_database_only_snippy).tofeaturesset(include_unknown_features=True)
    assert 632 == len(mutations)

    # Test only unknowns
    mutations = query(loaded_database_only_snippy).tofeaturesset(include_present_features=False,
                                                                 include_unknown_features=True)
    assert 521 == len(mutations)


def test_tofeaturesset_unique_all_selected(loaded_database_connection: DataIndexConnection):
    unique_mutations = query(loaded_database_connection).tofeaturesset(selection='unique')
    assert 111 == len(unique_mutations)

    # Include unknowns
    unique_mutations = query(loaded_database_connection).tofeaturesset(selection='unique',
                                                                       include_unknown_features=True)
    assert 632 == len(unique_mutations)

    # Include only unknowns
    unique_mutations = query(loaded_database_connection).tofeaturesset(selection='unique',
                                                                       include_present_features=False,
                                                                       include_unknown_features=True)
    assert 521 == len(unique_mutations)


def test_tofeaturesset_unique_none_selected(loaded_database_connection: DataIndexConnection):
    unique_mutations = query(loaded_database_connection).intersect(
        SampleSet.create_empty()).tofeaturesset(selection='unique')
    assert 0 == len(unique_mutations)


def test_tofeaturesset_unique_one_sample(loaded_database_only_snippy: DataIndexConnection):
    db = loaded_database_only_snippy.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()

    with open(data_dir / 'features_in_A_not_BC.txt', 'r') as fh:
        expected_set = {line.rstrip() for line in fh}

    f = {'reference:461:AAAT:G'}
    warnings.warn(f'Removing {f} from expected until I can figure out how to properly handle these')
    expected_set = expected_set - f

    sample_setA = SampleSet([sampleA.id])

    # No unknowns
    query_A = query(loaded_database_only_snippy).intersect(sample_setA)
    unique_mutations_A = query_A.tofeaturesset(selection='unique')
    assert 45 == len(unique_mutations_A)
    assert expected_set == unique_mutations_A

    # Include unknowns
    query_A = query(loaded_database_only_snippy).intersect(sample_setA)
    unique_mutations_A = query_A.tofeaturesset(selection='unique', include_unknown_features=True)
    assert 45 + 138 == len(unique_mutations_A)

    # Include only unknowns
    query_A = query(loaded_database_only_snippy).intersect(sample_setA)
    unique_mutations_A = query_A.tofeaturesset(selection='unique', include_unknown_features=True,
                                               include_present_features=False)
    assert 138 == len(unique_mutations_A)


def test_tofeaturesset_unique_two_samples(loaded_database_only_snippy: DataIndexConnection):
    db = loaded_database_only_snippy.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    with open(data_dir / 'features_in_BC_not_A.txt', 'r') as fh:
        expected_set = {line.rstrip() for line in fh}

    sample_setBC = SampleSet([sampleB.id, sampleC.id])

    # No include unknowns
    query_BC = query(loaded_database_only_snippy).intersect(sample_setBC)
    unique_mutations_BC = query_BC.tofeaturesset(selection='unique')
    assert 66 == len(unique_mutations_BC)
    assert expected_set == unique_mutations_BC

    # Include unknowns
    query_BC = query(loaded_database_only_snippy).intersect(sample_setBC)
    unique_mutations_BC = query_BC.tofeaturesset(selection='unique', include_unknown_features=True)
    assert 66 + 84 == len(unique_mutations_BC)

    # Include only unknowns
    query_BC = query(loaded_database_only_snippy).intersect(sample_setBC)
    unique_mutations_BC = query_BC.tofeaturesset(selection='unique', include_unknown_features=True,
                                                 include_present_features=False)
    assert 84 == len(unique_mutations_BC)
