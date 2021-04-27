import math

import pandas as pd
import pytest

from genomics_data_index.api.query.GenomicsDataIndex import GenomicsDataIndex
from genomics_data_index.api.query.SamplesQuery import SamplesQuery
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


def test_query_isin_sample_names(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    query_result = query(loaded_database_connection).isin('SampleB')
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


def test_query_isin_sample_names_multilple(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    query_result = query(loaded_database_connection).isin(['SampleA', 'SampleB'])
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


def test_query_single_mutation(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()

    query_result = query(loaded_database_connection).hasa(QueryFeatureMutation('reference:5061:G:A'))
    assert 1 == len(query_result)
    assert {sampleB.id} == set(query_result.sample_set)
    assert 9 == len(query_result.universe_set)


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


def test_query_custom_dataframe_isin_sample_names(loaded_database_connection: DataIndexConnection):
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
    assert {"dataframe(ids_col=[Sample ID]) AND isin(['SampleA', 'SampleC'])"} == set(query_result.toframe()['Query'].tolist())


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

    # By default, isa should select by 'name'
    sub_result = query_result.isa('SampleB')
    assert 1 == len(sub_result)
    assert 3 == len(sub_result.universe_set)
    assert {sampleB.id} == set(sub_result.sample_set)

    df = sub_result.toframe()
    assert ['SampleB'] == df['Sample Name'].tolist()
    assert {"dataframe(ids_col=[Sample ID]) AND isa('SampleB')"} == set(df['Query'].tolist())

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

    assert isinstance(query_result, TreeSamplesQuery)
    assert query_result.tree is not None

    assert {'SampleB', 'SampleC', 'genome'} == set(query_result.tree.get_leaf_names())


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
    query_result = query_result.join(data_frame=metadata_df, sample_ids_column='Sample ID')
    assert 2 == len(query_result)
    assert 3 == len(query_result.universe_set)
    assert isinstance(query_result, TreeSamplesQuery)

    df = query_result.toframe()

    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status', 'Color'] == df.columns.tolist()
    assert {'reference:839:C:G AND mutation_tree(genome) AND dataframe(ids_col=[Sample ID])'} == set(
        df['Query'].tolist())

    df = df.sort_values(['Sample Name'])
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert [sampleB.id, sampleC.id] == df['Sample ID'].tolist()
    assert ['green', 'blue'] == df['Color'].tolist()

    # I should still be able to perform within queries since I have a tree attached
    query_result = query_result.isin(['SampleB'], kind='mrca')

    assert 1 == len(query_result)
    assert 3 == len(query_result.universe_set)
    assert ['SampleB'] == query_result.tolist()


def test_within_constructed_tree(loaded_database_connection: DataIndexConnection):
    query_result = query(loaded_database_connection).hasa(
        'reference:839:C:G', kind='mutation').build_tree(
        kind='mutation', scope='genome', include_reference=True, extra_params='--seed 42 -m GTR')
    assert 2 == len(query_result)

    # subs/site
    df = query_result.isin('SampleC', kind='distance', distance=0.005, units='substitutions/site').toframe().sort_values('Sample Name')
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {'reference:839:C:G AND mutation_tree(genome) AND within(0.005 substitutions/site of SampleC)'
            } == set(df['Query'].tolist())

    # subs/site using within
    df = query_result.within('SampleC', distance=0.005, units='substitutions/site').toframe().sort_values('Sample Name')
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {'reference:839:C:G AND mutation_tree(genome) AND within(0.005 substitutions/site of SampleC)'
            } == set(df['Query'].tolist())

    # subs
    df = query_result.isin('SampleC', kind='distance', distance=26, units='substitutions').toframe().sort_values('Sample Name')
    assert 2 == len(df)
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {'reference:839:C:G AND mutation_tree(genome) AND within(26 substitutions of SampleC)'
            } == set(df['Query'].tolist())

    # should not include reference genome
    df = query_result.isin('SampleC', kind='distance', distance=100, units='substitutions').toframe().sort_values('Sample Name')
    assert 2 == len(df)
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {'reference:839:C:G AND mutation_tree(genome) AND within(100 substitutions of SampleC)'
            } == set(df['Query'].tolist())

    # should have only query sample
    df = query_result.isin('SampleC', kind='distance', distance=1, units='substitutions').toframe().sort_values('Sample Name')
    assert 1 == len(df)
    assert ['SampleC'] == df['Sample Name'].tolist()
    assert {'reference:839:C:G AND mutation_tree(genome) AND within(1 substitutions of SampleC)'
            } == set(df['Query'].tolist())

    # mrca
    df = query_result.isin(['SampleB', 'SampleC'], kind='mrca').toframe().sort_values('Sample Name')
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleB', 'SampleC'] == df['Sample Name'].tolist()
    assert {"reference:839:C:G AND mutation_tree(genome) AND within(mrca of ['SampleB', 'SampleC'])"
            } == set(df['Query'].tolist())

    # Sample Names isin()
    df = query_result.isin(['SampleA', 'SampleC']).toframe().sort_values('Sample Name')
    assert 2 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleA', 'SampleC'] == df['Sample Name'].tolist()
    assert {"reference:839:C:G AND mutation_tree(genome) AND isin(['SampleA', 'SampleC'])"
            } == set(df['Query'].tolist())

    # Sample Names isa()
    df = query_result.isa('SampleA').toframe().sort_values('Sample Name')
    assert 1 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
    assert ['SampleA'] == df['Sample Name'].tolist()
    assert {"reference:839:C:G AND mutation_tree(genome) AND isa('SampleA')"
            } == set(df['Query'].tolist())


def test_within_constructed_tree_larger_tree(loaded_database_connection: DataIndexConnection):
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

    # mrca A and B
    df = query_result.isin(['SampleA', 'SampleC'], kind='mrca').toframe().sort_values('Sample Name')
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
