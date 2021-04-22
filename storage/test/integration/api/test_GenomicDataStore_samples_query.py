import math
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd
import pytest

from storage.test.integration import snippy_all_dataframes, data_dir
from storage.variant.SampleSet import SampleSet
from storage.api.impl.TreeSamplesQuery import TreeSamplesQuery
from storage.api.SamplesQuery import SamplesQuery
from storage.api.GenomicDataStore import GenomicDataStore
from storage.connector.DataIndexConnection import DataIndexConnection
from storage.variant.model.QueryFeatureMLST import QueryFeatureMLST
from storage.variant.model.QueryFeatureMutation import QueryFeatureMutation
from storage.variant.model.db import Sample

# wrapper methods to simplify writing tests
def query(connection: DataIndexConnection, **kwargs) -> SamplesQuery:
    return GenomicDataStore(connection=connection).samples_query(**kwargs)


def test_connect():
    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)

        ds = GenomicDataStore.connect(database_connection='sqlite:///:memory:', database_dir=tmp_file)
        assert ds is not None
        assert ds.connection.reference_service is not None
        assert ds.connection.filesystem_storage.variation_dir.parent == tmp_file


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


def test_query_single_mutation(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    query_result = query(loaded_database_connection).has(QueryFeatureMutation('reference:5061:G:A'))
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


def test_query_single_mutation_complement(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    query_result = query(loaded_database_connection).has(QueryFeatureMutation('reference:5061:G:A'))
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)
    assert [sampleB.id] == query_result.tolist(names=False)
    assert ['SampleB'] == query_result.tolist(names=True)
    assert ['SampleB'] == query_result.tolist()

    query_result = query_result.complement()
    assert 8 == len(query_result)
    assert 9 == len(query_result.universe_set)

    assert sampleB.id not in query_result.sample_set

    assert {'SampleA', 'SampleC', 'CFSAN002349', 'CFSAN023463',
            '2014C-3598', '2014C-3599', '2014D-0067', '2014D-0068'} == set(query_result.tolist())

    df = query_result.toframe()
    assert {'reference:5061:G:A AND complement'} == set(df['Query'].tolist())
    assert {'SampleA', 'SampleC', 'CFSAN002349', 'CFSAN023463',
            '2014C-3598', '2014C-3599', '2014D-0067', '2014D-0068'} == set(df['Sample Name'].tolist())
    assert {'Present'} == set(df['Status'].tolist())


def test_query_single_mutation_summary(loaded_database_connection: DataIndexConnection):
    df = query(loaded_database_connection).has('reference:5061:G:A', kind='mutation').summary()
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

    query_result = query(loaded_database_connection).has(QueryFeatureMutation('reference:839:C:G'))
    assert 2 == len(query_result)
    assert {sampleB.id, sampleC.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


def test_query_single_mutation_two_samples_complement(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    query_result = query(loaded_database_connection).has(QueryFeatureMutation('reference:839:C:G'))
    assert 2 == len(query_result)
    assert {sampleB.id, sampleC.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)

    query_result = query_result.complement()
    assert 7 == len(query_result)
    assert sampleB.id not in query_result.sample_set
    assert sampleC.id not in query_result.sample_set
    assert sampleA.id in query_result.sample_set


def test_query_single_mutation_no_results(loaded_database_connection: DataIndexConnection):
    query_result = query(loaded_database_connection).has(QueryFeatureMutation('reference:1:1:A'))
    assert 0 == len(query_result)
    assert query_result.is_empty()
    assert 9 == len(query_result.universe_set)


def test_query_chained_mutation(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    query_result = query(loaded_database_connection).has(
        QueryFeatureMutation('reference:839:C:G')).has(
        QueryFeatureMutation('reference:5061:G:A'))
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


def test_query_chained_mutation_has_mutation(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    query_result = query(loaded_database_connection).has(
        'reference:839:C:G', kind='mutation').has(
        'reference:5061:G:A', kind='mutation')
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


def test_query_mlst_allele(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sample1 = db.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()
    sample2 = db.get_session().query(Sample).filter(Sample.name == 'CFSAN023463').one()

    query_result = query(loaded_database_connection).has(QueryFeatureMLST('lmonocytogenes:abcZ:1'))
    assert 2 == len(query_result)
    assert {sample1.id, sample2.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)

    assert {'CFSAN002349', 'CFSAN023463'} == set(query_result.tolist())
    assert {sample1.id, sample2.id} == set(query_result.tolist(names=False))


def test_query_chained_mlst_alleles(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sample1 = db.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()

    query_result = query(loaded_database_connection).has(
        QueryFeatureMLST('lmonocytogenes:abcZ:1')).has(
        QueryFeatureMLST('lmonocytogenes:lhkA:4'))
    assert 1 == len(query_result)
    assert {sample1.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


def test_query_chained_mlst_alleles_has_allele(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sample1 = db.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()

    query_result = query(loaded_database_connection) \
        .has('lmonocytogenes:abcZ:1', kind='mlst') \
        .has('lmonocytogenes:lhkA:4', kind='mlst')
    assert 1 == len(query_result)
    assert {sample1.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


@pytest.mark.skip
def test_query_chained_mlst_nucleotide(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sample1 = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    query_result = query(loaded_database_connection) \
        .has('reference:839:C:G', kind='mutation') \
        .has('lmonocytogenes:cat:12', kind='mlst')
    assert 1 == len(query_result)
    assert {sample1.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)

    assert ['SampleC'] == query_result.tolist()
    assert [sample1.id] == query_result.tolist(names=False)


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


def test_query_single_mutation_dataframe_include_all(loaded_database_connection: DataIndexConnection):
    df = query(loaded_database_connection).has(
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


def test_query_single_mutation_no_results_summary(loaded_database_connection: DataIndexConnection):
    df = query(loaded_database_connection).has('reference:1:1:A', kind='mutation').summary()
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

    query_result = query_result.has('reference:839:C:G', kind='mutation')
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
                         universe='dataframe',
                         data_frame=df,
                         sample_names_column='Samples')
    assert 3 == len(query_result)
    assert 3 == len(query_result.universe_set)

    query_result = query_result.has('reference:839:C:G', kind='mutation')

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
                         universe='dataframe',
                         data_frame=df,
                         sample_names_column='Samples')
    assert 3 == len(query_result)
    assert 3 == len(query_result.universe_set)

    query_result = query_result.has('reference:839:C:G', kind='mutation')
    df = query_result.toframe()

    assert 2 == len(df)
    assert 3 == len(query_result.universe_set)
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
                         universe='dataframe',
                         data_frame=df,
                         sample_names_column='Samples')

    assert 2 == len(df)
    assert 2 == len(query_result.universe_set)

    df = query_result.has('reference:839:C:G', kind='mutation').toframe()

    assert 1 == len(df)
    assert 2 == len(query_result.universe_set)
    assert {'Query', 'Sample Name', 'Sample ID', 'Status', 'Color', 'Samples'} == set(df.columns.tolist())

    df = df.sort_values(['Sample Name'])
    assert ['SampleC'] == df['Sample Name'].tolist()
    assert [sampleC.id] == df['Sample ID'].tolist()
    assert ['blue'] == df['Color'].tolist()
    assert {'reference:839:C:G'} == set(df['Query'].tolist())


def test_query_and_build_mutation_tree(loaded_database_connection: DataIndexConnection):
    query_result = query(loaded_database_connection).has('reference:839:C:G', kind='mutation')
    assert 2 == len(query_result)
    assert 9 == len(query_result.universe_set)

    query_result = query_result.build_tree(kind='mutation', scope='genome', include_reference=True)
    assert 2 == len(query_result)
    assert 9 == len(query_result.universe_set)

    assert isinstance(query_result, TreeSamplesQuery)
    assert query_result.tree is not None

    assert {'SampleB', 'SampleC', 'genome'} == set(query_result.tree.get_leaf_names())


def test_query_build_tree_and_query(loaded_database_connection: DataIndexConnection):
    query_result = query(loaded_database_connection).has('reference:839:C:G', kind='mutation')
    assert 2 == len(query_result)
    assert 9 == len(query_result.universe_set)

    query_result = query_result.build_tree(kind='mutation', scope='genome', include_reference=True)
    assert 2 == len(query_result)
    assert 9 == len(query_result.universe_set)

    query_result = query_result.has('reference:5061:G:A', kind='mutation')
    assert 1 == len(query_result)
    assert 9 == len(query_result.universe_set)

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
    df = query_result.within('SampleC', distance=0.005, units='substitutions/site').toframe().sort_values('Sample Name')
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {'reference:839:C:G AND mutation_tree(genome) AND within(0.005 substitutions/site of SampleC)'
            } == set(df['Query'].tolist())

    # subs
    df = query_result.within('SampleC', distance=26, units='substitutions').toframe().sort_values('Sample Name')
    assert 2 == len(df)
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {'reference:839:C:G AND mutation_tree(genome) AND within(26 substitutions of SampleC)'
            } == set(df['Query'].tolist())

    # should not include reference genome
    df = query_result.within('SampleC', distance=100, units='substitutions').toframe().sort_values('Sample Name')
    assert 2 == len(df)
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {'reference:839:C:G AND mutation_tree(genome) AND within(100 substitutions of SampleC)'
            } == set(df['Query'].tolist())

    # should have only query sample
    df = query_result.within('SampleC', distance=1, units='substitutions').toframe().sort_values('Sample Name')
    assert 1 == len(df)
    assert ['SampleC'] == df['Sample Name'].tolist()
    assert {'reference:839:C:G AND mutation_tree(genome) AND within(1 substitutions of SampleC)'
            } == set(df['Query'].tolist())

    # mrca
    df = query_result.within(['SampleB', 'SampleC'], kind='mrca').toframe().sort_values('Sample Name')
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {"reference:839:C:G AND mutation_tree(genome) AND within(mrca of ['SampleB', 'SampleC'])"
            } == set(df['Query'].tolist())


def test_within_constructed_tree_larger_tree(loaded_database_connection: DataIndexConnection):
    # Construct new tree with all the samples
    query_result = query(loaded_database_connection).build_tree(
        kind='mutation', scope='genome', include_reference=True, extra_params='--seed 42 -m GTR')
    assert 3 == len(query_result)

    # mrca B and C
    df = query_result.within(['SampleB', 'SampleC'], kind='mrca').toframe().sort_values('Sample Name')
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {"mutation_tree(genome) AND within(mrca of ['SampleB', 'SampleC'])"
            } == set(df['Query'].tolist())

    # mrca A and B
    df = query_result.within(['SampleA', 'SampleC'], kind='mrca').toframe().sort_values('Sample Name')
    assert 3 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleA', 'SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {"mutation_tree(genome) AND within(mrca of ['SampleA', 'SampleC'])"
            } == set(df['Query'].tolist())


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

    mutations_df = query(loaded_database_connection).has(
        'reference:839:C:G', kind='mutation').summary_features()
    mutations_df = mutations_df.sort_index()

    assert len(expected_df) == len(mutations_df)
    assert list(expected_df.index) == list(mutations_df.index)
    assert list(expected_df['Count']) == list(mutations_df['Count'])


def test_tofeaturesset_all(loaded_database_connection: DataIndexConnection):
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

    mutations = query(loaded_database_connection).tofeaturesset()

    assert 112 == len(mutations)
    assert set(expected_set) == set(mutations)


def test_tofeaturesset_unique_all_selected(loaded_database_connection: DataIndexConnection):
    unique_mutations = query(loaded_database_connection).tofeaturesset(selection='unique')
    assert 0 == len(unique_mutations)


def test_tofeaturesset_unique_one_sample(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    with open(data_dir / 'features_in_A_not_BC.txt', 'r') as fh:
        expected_set = {line for line in fh}

    sample_setA = SampleSet([sampleA.id])

    query_A = query(loaded_database_connection).intersect(sample_setA)
    unique_mutations_A = query_A.tofeaturesset(selection='unique')

    assert 46 == len(unique_mutations_A)
    assert expected_set == unique_mutations_A
