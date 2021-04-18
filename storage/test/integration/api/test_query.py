from tempfile import TemporaryDirectory
from pathlib import Path

from storage.connector.DataIndexConnection import DataIndexConnection
from storage.variant.model.QueryFeatureMutation import QueryFeatureMutation
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
