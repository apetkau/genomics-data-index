from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd
import pytest

from storage.api.query import query, connect
from storage.connector.DataIndexConnection import DataIndexConnection
from storage.variant.SampleSet import AllSampleSet
from storage.variant.model.QueryFeatureMLST import QueryFeatureMLST
from storage.variant.model.QueryFeatureMutation import QueryFeatureMutation
from storage.variant.model.db import Sample
from storage.api.impl.TreeSamplesQuery import TreeSamplesQuery


def test_connect():
    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)

        connection = connect(database_connection='sqlite:///:memory:', database_dir=tmp_file)
        assert connection is not None
        assert connection.reference_service is not None
        assert connection.filesystem_storage.variation_dir.parent == tmp_file


def test_initialized_query(loaded_database_connection: DataIndexConnection):
    initial_query = query(loaded_database_connection)

    assert len(initial_query) == 9


def test_query_single_mutation(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    query_result = query(loaded_database_connection).has(QueryFeatureMutation('reference:5061:G:A'))
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)


def test_query_single_mutation_two_samples(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    query_result = query(loaded_database_connection).has(QueryFeatureMutation('reference:839:C:G'))
    assert 2 == len(query_result)
    assert {sampleB.id, sampleC.id} == set(query_result.sample_set)


def test_query_single_mutation_no_results(loaded_database_connection: DataIndexConnection):
    query_result = query(loaded_database_connection).has(QueryFeatureMutation('reference:1:1:A'))
    assert 0 == len(query_result)
    assert query_result.is_empty()


def test_query_chained_mutation(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    query_result = query(loaded_database_connection).has(
        QueryFeatureMutation('reference:839:C:G')).has(
        QueryFeatureMutation('reference:5061:G:A'))
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)


def test_query_chained_mutation_has_mutation(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    query_result = query(loaded_database_connection).has(
        'reference:839:C:G', kind='mutation').has(
        'reference:5061:G:A', kind='mutation')
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)


def test_query_mlst_allele(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sample1 = db.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()
    sample2 = db.get_session().query(Sample).filter(Sample.name == 'CFSAN023463').one()

    query_result = query(loaded_database_connection).has(QueryFeatureMLST('lmonocytogenes:abcZ:1'))
    assert 2 == len(query_result)
    assert {sample1.id, sample2.id} == set(query_result.sample_set)


def test_query_chained_mlst_alleles(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sample1 = db.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()

    query_result = query(loaded_database_connection).has(
        QueryFeatureMLST('lmonocytogenes:abcZ:1')).has(
        QueryFeatureMLST('lmonocytogenes:lhkA:4'))
    assert 1 == len(query_result)
    assert {sample1.id} == set(query_result.sample_set)


def test_query_chained_mlst_alleles_has_allele(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sample1 = db.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()

    query_result = query(loaded_database_connection) \
        .has('lmonocytogenes:abcZ:1', kind='mlst') \
        .has('lmonocytogenes:lhkA:4', kind='mlst')
    assert 1 == len(query_result)
    assert {sample1.id} == set(query_result.sample_set)


@pytest.mark.skip
def test_query_chained_mlst_nucleotide(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sample1 = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    query_result = query(loaded_database_connection) \
        .has('reference:839:C:G', kind='mutation') \
        .has('lmonocytogenes:cat:12', kind='mlst')
    assert 1 == len(query_result)
    assert {sample1.id} == set(query_result.sample_set)


def test_query_single_mutation_dataframe(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    df = query(loaded_database_connection).has('reference:839:C:G', kind='mutation').toframe()

    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()

    df = df.sort_values(['Sample Name'])
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert [sampleB.id, sampleC.id] == df['Sample ID'].tolist()
    assert {'reference:839:C:G'} == set(df['Query'].tolist())


def test_query_chained_allele_dataframe(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sample1 = db.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()

    df = query(loaded_database_connection) \
        .has('lmonocytogenes:abcZ:1', kind='mlst') \
        .has('lmonocytogenes:lhkA:4', kind='mlst').toframe()

    assert 1 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()

    df = df.sort_values(['Sample Name'])
    assert ['CFSAN002349'] == df['Sample Name'].tolist()
    assert [sample1.id] == df['Sample ID'].tolist()
    assert {'lmonocytogenes:abcZ:1 AND lmonocytogenes:lhkA:4'} == set(df['Query'].tolist())


def test_query_single_mutation_no_results_toframe(loaded_database_connection: DataIndexConnection):
    df = query(loaded_database_connection).has('reference:1:1:A', kind='mutation').toframe()
    assert 0 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()


def test_join_custom_dataframe_no_query(loaded_database_connection: DataIndexConnection):
    df = pd.DataFrame([
        [1, 'red'],
        [2, 'green'],
        [3, 'blue']
    ], columns=['Sample ID', 'Color'])

    query_result = query(loaded_database_connection, data_frame=df, sample_ids_column='Sample ID')
    assert 3 == len(query_result)
    assert {1, 2, 3} == set(query_result.sample_set)


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
                         data_frame=df,
                         sample_ids_column='Sample ID')
    query_result = query_result.has('reference:839:C:G', kind='mutation')
    assert 2 == len(query_result)
    assert {sampleB.id, sampleC.id} == set(query_result.sample_set)

    df = query_result.toframe()

    assert 2 == len(df)
    assert {'Query', 'Sample Name', 'Sample ID', 'Status', 'Color'} == set(df.columns.tolist())

    df = df.sort_values(['Sample Name'])
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert [sampleB.id, sampleC.id] == df['Sample ID'].tolist()
    assert ['green', 'blue'] == df['Color'].tolist()
    assert {'reference:839:C:G'} == set(df['Query'].tolist())


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
                         data_frame=df,
                         sample_names_column='Samples')
    query_result = query_result.has('reference:839:C:G', kind='mutation')
    assert 2 == len(query_result)
    assert {sampleB.id, sampleC.id} == set(query_result.sample_set)

    df = query_result.toframe()

    assert 2 == len(df)
    assert {'Query', 'Sample Name', 'Sample ID', 'Status', 'Color', 'Samples'} == set(df.columns.tolist())

    df = df.sort_values(['Sample Name'])
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert [sampleB.id, sampleC.id] == df['Sample ID'].tolist()
    assert ['green', 'blue'] == df['Color'].tolist()
    assert {'reference:839:C:G'} == set(df['Query'].tolist())


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
                         data_frame=df,
                         sample_names_column='Samples')
    df = query_result.has('reference:839:C:G', kind='mutation').toframe()

    assert 2 == len(df)
    assert {'Query', 'Sample Name', 'Sample ID', 'Status', 'Color', 'Samples'} == set(df.columns.tolist())

    df = df.sort_values(['Sample Name'])
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert [sampleB.id, sampleC.id] == df['Sample ID'].tolist()
    assert ['green', 'blue'] == df['Color'].tolist()
    assert {'reference:839:C:G'} == set(df['Query'].tolist())


def test_join_custom_dataframe_missing_sample_names(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    df = pd.DataFrame([
        ['SampleA', 'red'],
        ['SampleC', 'blue'],
    ], columns=['Samples', 'Color'])

    query_result = query(loaded_database_connection,
                         data_frame=df,
                         sample_names_column='Samples')
    df = query_result.has('reference:839:C:G', kind='mutation').toframe()

    assert 1 == len(df)
    assert {'Query', 'Sample Name', 'Sample ID', 'Status', 'Color', 'Samples'} == set(df.columns.tolist())

    df = df.sort_values(['Sample Name'])
    assert ['SampleC'] == df['Sample Name'].tolist()
    assert [sampleC.id] == df['Sample ID'].tolist()
    assert ['blue'] == df['Color'].tolist()
    assert {'reference:839:C:G'} == set(df['Query'].tolist())


def test_query_and_build_mutation_tree(loaded_database_connection: DataIndexConnection):
    query_result = query(loaded_database_connection).has('reference:839:C:G', kind='mutation')
    assert 2 == len(query_result)

    query_result = query_result.build_tree(kind='mutation', scope='genome', include_reference=True)
    assert 2 == len(query_result)

    assert isinstance(query_result, TreeSamplesQuery)
    assert query_result.tree is not None

    assert {'SampleB', 'SampleC', 'genome'} == set(query_result.tree.get_leaf_names())


def test_query_build_tree_and_query(loaded_database_connection: DataIndexConnection):
    query_result = query(loaded_database_connection).has('reference:839:C:G', kind='mutation')
    assert 2 == len(query_result)

    query_result = query_result.build_tree(kind='mutation', scope='genome', include_reference=True)
    assert 2 == len(query_result)

    query_result = query_result.has('reference:5061:G:A', kind='mutation')
    assert 1 == len(query_result)

    # Tree should still be complete
    assert {'SampleB', 'SampleC', 'genome'} == set(query_result.tree.get_leaf_names())


def test_query_build_tree_dataframe(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    df = query(loaded_database_connection).has(
        'reference:839:C:G', kind='mutation').build_tree(
        kind='mutation', scope='genome', include_reference=True).has(
        'reference:5061:G:A', kind='mutation').toframe()

    assert 1 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()

    assert ['SampleB'] == df['Sample Name'].tolist()
    assert ['Present'] == df['Status'].tolist()
    assert [sampleB.id] == df['Sample ID'].tolist()
    assert ['reference:839:C:G AND mutation_tree(genome) AND reference:5061:G:A'] == df['Query'].tolist()


def test_within_constructed_tree(loaded_database_connection: DataIndexConnection):
    query_result = query(loaded_database_connection).has(
        'reference:839:C:G', kind='mutation').build_tree(
        kind='mutation', scope='genome', include_reference=True, extra_params='--seed 42 -m GTR')
    assert 2 == len(query_result)

    # subs/site
    df = query_result.within(0.005, 'SampleC', units='substitutions/site').toframe().sort_values('Sample Name')
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {'reference:839:C:G AND mutation_tree(genome) AND within(0.005 substitutions/site of SampleC)'
            } == set(df['Query'].tolist())

    # subs
    df = query_result.within(26, 'SampleC', units='substitutions').toframe().sort_values('Sample Name')
    assert 2 == len(df)
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {'reference:839:C:G AND mutation_tree(genome) AND within(26 substitutions of SampleC)'
            } == set(df['Query'].tolist())

    # should not include reference genome
    df = query_result.within(100, 'SampleC', units='substitutions').toframe().sort_values('Sample Name')
    assert 2 == len(df)
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {'reference:839:C:G AND mutation_tree(genome) AND within(100 substitutions of SampleC)'
            } == set(df['Query'].tolist())

    # should have only query sample
    df = query_result.within(1, 'SampleC', units='substitutions').toframe().sort_values('Sample Name')
    assert 1 == len(df)
    assert ['SampleC'] == df['Sample Name'].tolist()
    assert {'reference:839:C:G AND mutation_tree(genome) AND within(1 substitutions of SampleC)'
            } == set(df['Query'].tolist())
