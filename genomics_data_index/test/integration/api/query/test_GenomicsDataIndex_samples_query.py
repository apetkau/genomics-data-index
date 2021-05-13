import math
from typing import cast

import pandas as pd
import pytest
from ete3 import Tree, ClusterTree

from genomics_data_index.api.query.GenomicsDataIndex import GenomicsDataIndex
from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.impl.ExperimentalTreeSamplesQuery import ExperimentalTreeSamplesQuery
from genomics_data_index.api.query.impl.MutationTreeSamplesQuery import MutationTreeSamplesQuery
from genomics_data_index.api.query.impl.TreeSamplesQuery import TreeSamplesQuery
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.QueryFeatureMLST import QueryFeatureMLST
from genomics_data_index.storage.model.QueryFeatureMutation import QueryFeatureMutation
from genomics_data_index.storage.model.db import Sample
from genomics_data_index.test.integration import snippy_all_dataframes, data_dir


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


def test_query_isin_samples(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    query_result = query(loaded_database_connection).isin('SampleB')
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


def test_query_isin_samples_no_exist(loaded_database_connection: DataIndexConnection):
    query_result = query(loaded_database_connection).isin('no_exist')
    assert 0 == len(query_result)
    assert 9 == len(query_result.universe_set)


def test_query_isin_sample_set_single(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleBSet = SampleSet([sampleB.id])

    query_result = query(loaded_database_connection).isin(sampleBSet)
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


def test_query_isin_samples_query_single(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    query_result_B = query(loaded_database_connection).isin('SampleB')

    query_result = query(loaded_database_connection).isin(query_result_B)
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


def test_query_isin_samples_query_no_matches(loaded_database_connection: DataIndexConnection):
    query_result_empty = query(loaded_database_connection).isin('no_exist')

    query_result = query(loaded_database_connection).isin(query_result_empty)
    assert 0 == len(query_result)
    assert 9 == len(query_result.universe_set)


def test_query_and_or(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    query_resultA = query(loaded_database_connection).isa('SampleA')
    assert 1 == len(query_resultA)
    assert {sampleA.id} == set(query_resultA.sample_set)
    assert 9 == len(query_resultA.universe_set)

    query_resultB = query(loaded_database_connection).isa('SampleB')
    assert 1 == len(query_resultB)
    assert {sampleB.id} == set(query_resultB.sample_set)
    assert 9 == len(query_resultB.universe_set)

    query_resultC = query(loaded_database_connection).isa('SampleC')
    assert 1 == len(query_resultC)
    assert {sampleC.id} == set(query_resultC.sample_set)
    assert 9 == len(query_resultC.universe_set)

    query_resultAB = query(loaded_database_connection).isin(['SampleA', 'SampleB'])
    assert 2 == len(query_resultAB)
    assert {sampleA.id, sampleB.id} == set(query_resultAB.sample_set)
    assert 9 == len(query_resultAB.universe_set)

    # Test AND
    query_resultA_and_AB = query_resultA.and_(query_resultAB)
    assert 1 == len(query_resultA_and_AB)
    assert {sampleA.id} == set(query_resultA_and_AB.sample_set)
    assert 9 == len(query_resultA_and_AB.universe_set)

    query_resultAB_and_A = query_resultAB.and_(query_resultA)
    assert 1 == len(query_resultAB_and_A)
    assert {sampleA.id} == set(query_resultAB_and_A.sample_set)
    assert 9 == len(query_resultAB_and_A.universe_set)

    query_resultA_and_B = query_resultA.and_(query_resultB)
    assert 0 == len(query_resultA_and_B)
    assert 9 == len(query_resultA_and_B.universe_set)

    # Test OR
    query_resultA_or_AB = query_resultA.or_(query_resultAB)
    assert 2 == len(query_resultA_or_AB)
    assert {sampleA.id, sampleB.id} == set(query_resultA_or_AB.sample_set)
    assert 9 == len(query_resultA_or_AB.universe_set)

    query_resultAB_or_A = query_resultAB.or_(query_resultA)
    assert 2 == len(query_resultAB_or_A)
    assert {sampleA.id, sampleB.id} == set(query_resultAB_or_A.sample_set)
    assert 9 == len(query_resultAB_or_A.universe_set)

    query_resultA_or_B = query_resultA.or_(query_resultB)
    assert 2 == len(query_resultA_or_B)
    assert {sampleA.id, sampleB.id} == set(query_resultA_or_B.sample_set)
    assert 9 == len(query_resultA_or_B.universe_set)

    query_resultA_or_B_or_C = query_resultA.or_(query_resultB).or_(query_resultC)
    assert 3 == len(query_resultA_or_B_or_C)
    assert {sampleA.id, sampleB.id, sampleC.id} == set(query_resultA_or_B_or_C.sample_set)
    assert 9 == len(query_resultA_or_B_or_C.universe_set)


def test_query_reset_universe(loaded_database_connection: DataIndexConnection):
    query_result = query(loaded_database_connection).isin('SampleB')
    assert 1 == len(query_result)
    assert 9 == len(query_result.universe_set)

    query_result = query_result.reset_universe()
    assert 1 == len(query_result)
    assert 1 == len(query_result.universe_set)


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


def test_query_isa_sample_no_exist(loaded_database_connection: DataIndexConnection):
    query_result = query(loaded_database_connection).isa('no_exist')
    assert 0 == len(query_result)
    assert 9 == len(query_result.universe_set)


def test_query_isa_samples_query(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    query_result_B = query(loaded_database_connection).isa('SampleB')

    query_result = query(loaded_database_connection).isa(query_result_B)
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


def test_query_isin_kmer(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    query_result = query(loaded_database_connection).isin('SampleA', kind='distance', distance=1.0,
                                                          units='kmer_jaccard')
    assert 3 == len(query_result)
    assert {sampleA.id, sampleB.id, sampleC.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)
    assert "isin_kmer_jaccard('SampleA', dist=1.0, k=31)" == query_result.query_expression()


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
    assert "isin_kmer_jaccard(<SamplesQueryIndex[11% (1/9) samples]>, dist=1.0, k=31)" == query_result.query_expression()


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
    assert "isin_kmer_jaccard(<SamplesQueryIndex[0% (0/9) samples]>, dist=1.0, k=31)" == query_result.query_expression()


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


def test_query_single_mutation(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    query_result = query(loaded_database_connection).hasa(QueryFeatureMutation('reference:5061:G:A'))
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)
    assert not query_result.has_tree()


def test_query_single_mutation_complement(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    query_result = query(loaded_database_connection).hasa(QueryFeatureMutation('reference:5061:G:A'))
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)
    assert [sampleB.id] == query_result.tolist(names=False)
    assert ['SampleB'] == query_result.tolist(names=True)
    assert ['SampleB'] == query_result.tolist()

    query_result = query_result.complement()
    assert 8 == len(query_result)
    assert 9 == len(query_result.universe_set)
    assert not query_result.has_tree()

    assert sampleB.id not in query_result.sample_set

    assert {'SampleA', 'SampleC', 'CFSAN002349', 'CFSAN023463',
            '2014C-3598', '2014C-3599', '2014D-0067', '2014D-0068'} == set(query_result.tolist())

    df = query_result.toframe()
    assert {'reference:5061:G:A AND complement'} == set(df['Query'].tolist())
    assert {'SampleA', 'SampleC', 'CFSAN002349', 'CFSAN023463',
            '2014C-3598', '2014C-3599', '2014D-0067', '2014D-0068'} == set(df['Sample Name'].tolist())
    assert {'Present'} == set(df['Status'].tolist())


def test_query_single_mutation_summary(loaded_database_connection: DataIndexConnection):
    df = query(loaded_database_connection).hasa('reference:5061:G:A', kind='mutation').summary()
    assert 1 == len(df)
    assert ['Query', 'Present', 'Absent', 'Unknown', 'Total',
            '% Present', '% Absent', '% Unknown'] == df.columns.tolist()
    assert 'reference:5061:G:A' == df.iloc[0]['Query']
    assert 1 == df.iloc[0]['Present']
    assert 8 == df.iloc[0]['Absent']
    assert pd.isna(df.iloc[0]['Unknown'])
    assert 9 == df.iloc[0]['Total']
    assert math.isclose((1 / 9) * 100, df.iloc[0]['% Present'])
    assert math.isclose((8 / 9) * 100, df.iloc[0]['% Absent'])
    assert pd.isna(df.iloc[0]['% Unknown'])


def test_query_single_mutation_two_samples(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    query_result = query(loaded_database_connection).hasa(QueryFeatureMutation('reference:839:C:G'))
    assert 2 == len(query_result)
    assert {sampleB.id, sampleC.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


def test_query_single_mutation_two_samples_kmer_one_sample(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    query_result = query(loaded_database_connection).hasa(QueryFeatureMutation('reference:839:C:G'))
    assert 2 == len(query_result)
    assert {sampleB.id, sampleC.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)

    query_result = query_result.isin('SampleA', kind='distance', distance=0.5,
                                     units='kmer_jaccard')
    assert 1 == len(query_result)
    assert {sampleC.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


def test_query_single_mutation_default_kind(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    # Should default to kind='mutation'
    query_result = query(loaded_database_connection).hasa('reference:839:C:G')
    assert 2 == len(query_result)
    assert {sampleB.id, sampleC.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


def test_query_single_mutation_two_samples_complement(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    query_result = query(loaded_database_connection).hasa(QueryFeatureMutation('reference:839:C:G'))
    assert 2 == len(query_result)
    assert {sampleB.id, sampleC.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)

    query_result = query_result.complement()
    assert 7 == len(query_result)
    assert sampleB.id not in query_result.sample_set
    assert sampleC.id not in query_result.sample_set
    assert sampleA.id in query_result.sample_set


def test_query_single_mutation_no_results(loaded_database_connection: DataIndexConnection):
    query_result = query(loaded_database_connection).hasa(QueryFeatureMutation('reference:1:1:A'))
    assert 0 == len(query_result)
    assert query_result.is_empty()
    assert 9 == len(query_result.universe_set)


def test_query_chained_mutation(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    query_result = query(loaded_database_connection).hasa(
        QueryFeatureMutation('reference:839:C:G')).hasa(
        QueryFeatureMutation('reference:5061:G:A'))
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
    sample1 = db.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()
    sample2 = db.get_session().query(Sample).filter(Sample.name == 'CFSAN023463').one()

    query_result = query(loaded_database_connection).hasa(QueryFeatureMLST('lmonocytogenes:abcZ:1'))
    assert 2 == len(query_result)
    assert {sample1.id, sample2.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)

    assert {'CFSAN002349', 'CFSAN023463'} == set(query_result.tolist())
    assert {sample1.id, sample2.id} == set(query_result.tolist(names=False))


def test_query_chained_mlst_alleles(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sample1 = db.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()

    query_result = query(loaded_database_connection).hasa(
        QueryFeatureMLST('lmonocytogenes:abcZ:1')).hasa(
        QueryFeatureMLST('lmonocytogenes:lhkA:4'))
    assert 1 == len(query_result)
    assert {sample1.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


def test_query_chained_mlst_alleles_has_allele(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sample1 = db.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()

    query_result = query(loaded_database_connection) \
        .hasa('lmonocytogenes:abcZ:1', kind='mlst') \
        .hasa('lmonocytogenes:lhkA:4', kind='mlst')
    assert 1 == len(query_result)
    assert {sample1.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


@pytest.mark.skip
def test_query_chained_mlst_nucleotide(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sample1 = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    query_result = query(loaded_database_connection) \
        .hasa('reference:839:C:G', kind='mutation') \
        .hasa('lmonocytogenes:cat:12', kind='mlst')
    assert 1 == len(query_result)
    assert {sample1.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)

    assert ['SampleC'] == query_result.tolist()
    assert [sample1.id] == query_result.tolist(names=False)


def test_query_single_mutation_dataframe(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    df = query(loaded_database_connection).hasa('reference:839:C:G', kind='mutation').toframe()

    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()

    df = df.sort_values(['Sample Name'])
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert [sampleB.id, sampleC.id] == df['Sample ID'].tolist()
    assert {'reference:839:C:G'} == set(df['Query'].tolist())


def test_query_single_mutation_dataframe_include_all(loaded_database_connection: DataIndexConnection):
    df = query(loaded_database_connection).hasa(
        'reference:839:C:G', kind='mutation').toframe(exclude_absent=False)

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


def test_query_chained_allele_dataframe(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sample1 = db.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()

    df = query(loaded_database_connection) \
        .hasa('lmonocytogenes:abcZ:1', kind='mlst') \
        .hasa('lmonocytogenes:lhkA:4', kind='mlst').toframe()

    assert 1 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()

    df = df.sort_values(['Sample Name'])
    assert ['CFSAN002349'] == df['Sample Name'].tolist()
    assert [sample1.id] == df['Sample ID'].tolist()
    assert {'lmonocytogenes:abcZ:1 AND lmonocytogenes:lhkA:4'} == set(df['Query'].tolist())


def test_query_single_mutation_no_results_toframe(loaded_database_connection: DataIndexConnection):
    df = query(loaded_database_connection).hasa('reference:1:1:A', kind='mutation').toframe()
    assert 0 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()


def test_query_single_mutation_no_results_summary(loaded_database_connection: DataIndexConnection):
    df = query(loaded_database_connection).hasa('reference:1:1:A', kind='mutation').summary()
    assert 1 == len(df)
    assert ['Query', 'Present', 'Absent', 'Unknown', 'Total',
            '% Present', '% Absent', '% Unknown'] == df.columns.tolist()
    assert 'reference:1:1:A' == df.iloc[0]['Query']
    assert 0 == df.iloc[0]['Present']
    assert 9 == df.iloc[0]['Absent']
    assert pd.isna(df.iloc[0]['Unknown'])
    assert 9 == df.iloc[0]['Total']
    assert math.isclose((0 / 9) * 100, df.iloc[0]['% Present'])
    assert math.isclose((9 / 9) * 100, df.iloc[0]['% Absent'])
    assert pd.isna(df.iloc[0]['% Unknown'])


def test_all_samples_summary(loaded_database_connection: DataIndexConnection):
    df = query(loaded_database_connection).summary()
    assert 1 == len(df)
    assert ['Query', 'Present', 'Absent', 'Unknown', 'Total',
            '% Present', '% Absent', '% Unknown'] == df.columns.tolist()
    assert '' == df.iloc[0]['Query']
    assert 9 == df.iloc[0]['Present']
    assert 0 == df.iloc[0]['Absent']
    assert pd.isna(df.iloc[0]['Unknown'])
    assert 9 == df.iloc[0]['Total']
    assert math.isclose((9 / 9) * 100, df.iloc[0]['% Present'])
    assert math.isclose((0 / 9) * 100, df.iloc[0]['% Absent'])
    assert pd.isna(df.iloc[0]['% Unknown'])


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


def test_join_custom_dataframe_single_query_sample_names(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
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

    query_result = query_result.hasa('reference:839:C:G', kind='mutation')

    assert 2 == len(query_result)
    assert 3 == len(query_result.universe_set)
    assert {sampleB.id, sampleC.id} == set(query_result.sample_set)

    df = query_result.toframe()

    assert 2 == len(df)
    assert {'Query', 'Sample Name', 'Sample ID', 'Status', 'Color', 'Samples'} == set(df.columns.tolist())

    df = df.sort_values(['Sample Name'])
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert [sampleB.id, sampleC.id] == df['Sample ID'].tolist()
    assert ['green', 'blue'] == df['Color'].tolist()
    assert {'dataframe(names_col=[Samples]) AND reference:839:C:G'} == set(df['Query'].tolist())


def test_join_custom_dataframe_extra_sample_names(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
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

    query_result = query_result.hasa('reference:839:C:G', kind='mutation')
    df = query_result.toframe()

    assert 2 == len(df)
    assert 3 == len(query_result.universe_set)
    assert {'Query', 'Sample Name', 'Sample ID', 'Status', 'Color', 'Samples'} == set(df.columns.tolist())

    df = df.sort_values(['Sample Name'])
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert [sampleB.id, sampleC.id] == df['Sample ID'].tolist()
    assert ['green', 'blue'] == df['Color'].tolist()
    assert {'dataframe(names_col=[Samples]) AND reference:839:C:G'} == set(df['Query'].tolist())


def test_join_custom_dataframe_missing_sample_names(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
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

    df = query_result.hasa('reference:839:C:G', kind='mutation').toframe()

    assert 1 == len(df)
    assert 2 == len(query_result.universe_set)
    assert {'Query', 'Sample Name', 'Sample ID', 'Status', 'Color', 'Samples'} == set(df.columns.tolist())

    df = df.sort_values(['Sample Name'])
    assert ['SampleC'] == df['Sample Name'].tolist()
    assert [sampleC.id] == df['Sample ID'].tolist()
    assert ['blue'] == df['Color'].tolist()
    assert {'dataframe(names_col=[Samples]) AND reference:839:C:G'} == set(df['Query'].tolist())


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

    query_result = query_result.hasa('reference:839:C:G', kind='mutation')
    assert 2 == len(query_result)
    assert 9 == len(query_result.universe_set)
    assert {sampleB.id, sampleC.id} == set(query_result.sample_set)

    df = query_result.toframe()

    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert {'reference:839:C:G'} == set(df['Query'].tolist())

    # Now join data frame
    query_result = query_result.join(data_frame=metadata_df, sample_ids_column='Sample ID')
    assert 2 == len(query_result)
    assert 3 == len(query_result.universe_set)

    df = query_result.toframe()

    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status', 'Color'] == df.columns.tolist()
    assert {'reference:839:C:G AND dataframe(ids_col=[Sample ID])'} == set(df['Query'].tolist())

    df = df.sort_values(['Sample Name'])
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert [sampleB.id, sampleC.id] == df['Sample ID'].tolist()
    assert ['green', 'blue'] == df['Color'].tolist()


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
            '<MutationTreeSamplesQuery[11% (1/9) samples]>)'
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
            "<MutationTreeSamplesQuery[11% (1/9) samples]>)"
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
            "<MutationTreeSamplesQuery[22% (2/9) samples]>)"
            } == set(df['Query'].tolist())

    # mrca of samples query, single sample as result
    query_result_B = query_result.isin(['SampleB'], kind='samples')
    df = query_result.isin(query_result_B, kind='mrca').toframe().sort_values('Sample Name')
    assert 1 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleB'] == df['Sample Name'].tolist()
    assert {"reference:839:C:G AND mutation_tree(genome) AND within(mrca of "
            "<MutationTreeSamplesQuery[11% (1/9) samples]>)"
            } == set(df['Query'].tolist())

    # mrca of samples query, empty_results
    query_result_empty = query_result.isin([], kind='samples')
    df = query_result.isin(query_result_empty, kind='mrca').toframe().sort_values('Sample Name')
    assert 0 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()

    # Samples isin()
    df = query_result.isin(['SampleA', 'SampleC']).toframe().sort_values('Sample Name')
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleA', 'SampleC'] == df['Sample Name'].tolist()
    assert {"reference:839:C:G AND mutation_tree(genome) AND isin_samples(['SampleA', 'SampleC'])"
            } == set(df['Query'].tolist())

    # Samples isin from samples query()
    query_result_AC = query_result.isin(['SampleA', 'SampleC'], kind='samples')
    df = query_result.isin(query_result_AC).toframe().sort_values('Sample Name')
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleA', 'SampleC'] == df['Sample Name'].tolist()
    assert {
               "reference:839:C:G AND mutation_tree(genome) AND isin_samples(<MutationTreeSamplesQuery[22% (2/9) samples]>)"
           } == set(df['Query'].tolist())

    # Sample isa()
    df = query_result.isa('SampleA').toframe().sort_values('Sample Name')
    assert 1 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleA'] == df['Sample Name'].tolist()
    assert {"reference:839:C:G AND mutation_tree(genome) AND isa_sample('SampleA')"
            } == set(df['Query'].tolist())

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
            "<MutationTreeSamplesQuery[22% (2/9) samples]>)"
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
            "<MutationTreeSamplesQuery[22% (2/9) samples]>)"
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
            "<MutationTreeSamplesQuery[22% (2/9) samples]>)"
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
            "<MutationTreeSamplesQuery[22% (2/9) samples]>)"
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
            '<MutationTreeSamplesQuery[11% (1/9) samples]>)'
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


def test_summary_features_all(loaded_database_connection: DataIndexConnection):
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

    mutations_df = query(loaded_database_connection).summary_features()
    mutations_df = mutations_df.sort_index()

    assert len(expected_df) == len(mutations_df)
    assert list(expected_df.index) == list(mutations_df.index)
    assert list(expected_df['Count']) == list(mutations_df['Count'])


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

    mutations_df = query(loaded_database_connection).hasa(
        'reference:839:C:G', kind='mutation').summary_features()
    mutations_df = mutations_df.sort_index()

    assert len(expected_df) == len(mutations_df)
    assert list(expected_df.index) == list(mutations_df.index)
    assert list(expected_df['Count']) == list(mutations_df['Count'])


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
    expected_set = set(expected_df.index)

    mutations = query(loaded_database_only_snippy).tofeaturesset()

    assert 112 == len(mutations)
    assert set(expected_set) == set(mutations)


def test_tofeaturesset_unique_all_selected(loaded_database_connection: DataIndexConnection):
    unique_mutations = query(loaded_database_connection).tofeaturesset(selection='unique')
    assert 112 == len(unique_mutations)


def test_tofeaturesset_unique_none_selected(loaded_database_connection: DataIndexConnection):
    unique_mutations = query(loaded_database_connection).intersect(
        SampleSet.create_empty()).tofeaturesset(selection='unique')
    assert 0 == len(unique_mutations)


def test_tofeaturesset_unique_one_sample(loaded_database_only_snippy: DataIndexConnection):
    db = loaded_database_only_snippy.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()

    with open(data_dir / 'features_in_A_not_BC.txt', 'r') as fh:
        expected_set = {line.rstrip() for line in fh}

    sample_setA = SampleSet([sampleA.id])

    query_A = query(loaded_database_only_snippy).intersect(sample_setA)
    unique_mutations_A = query_A.tofeaturesset(selection='unique')

    assert 46 == len(unique_mutations_A)
    assert expected_set == unique_mutations_A


def test_tofeaturesset_unique_two_samples(loaded_database_only_snippy: DataIndexConnection):
    db = loaded_database_only_snippy.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    with open(data_dir / 'features_in_BC_not_A.txt', 'r') as fh:
        expected_set = {line.rstrip() for line in fh}

    sample_setBC = SampleSet([sampleB.id, sampleC.id])

    query_BC = query(loaded_database_only_snippy).intersect(sample_setBC)
    unique_mutations_BC = query_BC.tofeaturesset(selection='unique')

    assert 66 == len(unique_mutations_BC)
    assert expected_set == unique_mutations_BC
