import pytest
from tempfile import TemporaryDirectory
from pathlib import Path

from storage.connector.DataIndexConnection import DataIndexConnection
from storage.variant.model.QueryFeatureMutation import QueryFeatureMutation
from storage.variant.model.QueryFeatureMLST import QueryFeatureMLST
from storage.api.query import query, connect
from storage.variant.model.db import Sample
from storage.variant.SampleSet import AllSampleSet


def test_connect():
    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)

        connection = connect(database_connection='sqlite:///:memory:', database_dir=tmp_file)
        assert connection is not None
        assert connection.reference_service is not None
        assert connection.filesystem_storage.variation_dir.parent == tmp_file


def test_initialized_query(loaded_database_connection: DataIndexConnection):
    initial_query = query(loaded_database_connection)

    assert len(initial_query) == 2**32
    assert isinstance(initial_query.sample_set, AllSampleSet)


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

    query_result = query(loaded_database_connection).has_mutation(
        'reference:839:C:G').has_mutation(
        'reference:5061:G:A')
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
        .has_allele('lmonocytogenes:abcZ:1') \
        .has_allele('lmonocytogenes:lhkA:4')
    assert 1 == len(query_result)
    assert {sample1.id} == set(query_result.sample_set)


@pytest.mark.skip
def test_query_chained_mlst_nucleotide(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sample1 = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    query_result = query(loaded_database_connection) \
        .has_mutation('reference:839:C:G') \
        .has_allele('lmonocytogenes:cat:12')
    assert 1 == len(query_result)
    assert {sample1.id} == set(query_result.sample_set)


def test_query_single_mutation_dataframe(loaded_database_connection: DataIndexConnection):
    db = loaded_database_connection.database
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()

    df = query(loaded_database_connection).has_mutation('reference:839:C:G').toframe()

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
        .has_allele('lmonocytogenes:abcZ:1') \
        .has_allele('lmonocytogenes:lhkA:4').toframe()

    assert 1 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()

    df = df.sort_values(['Sample Name'])
    assert ['CFSAN002349'] == df['Sample Name'].tolist()
    assert [sample1.id] == df['Sample ID'].tolist()
    assert {'lmonocytogenes:abcZ:1 AND lmonocytogenes:lhkA:4'} == set(df['Query'].tolist())


def test_query_single_mutation_no_results_toframe(loaded_database_connection: DataIndexConnection):
    df = query(loaded_database_connection).has_mutation('reference:1:1:A').toframe()
    assert 0 == len(df)
    assert ['Query', 'Sample Name', 'Sample ID', 'Status'] == df.columns.tolist()
