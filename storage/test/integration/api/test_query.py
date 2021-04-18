from tempfile import TemporaryDirectory
from pathlib import Path

from storage.connector.DataIndexConnection import DataIndexConnection
from storage.variant.model.QueryFeatureMutation import QueryFeatureMutation
from storage.api.query import query, connect


def test_connect():
    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)

        connection = connect(database_connection='sqlite:///:memory:', database_dir=tmp_file)
        assert connection is not None
        assert connection.reference_service is not None
        assert connection.filesystem_storage.variation_dir.parent == tmp_file


def test_initialized_query(loaded_database_connection: DataIndexConnection):
    assert query(loaded_database_connection).sample_set is None


def test_query_single_mutation(loaded_database_connection: DataIndexConnection):
    query_result = query(loaded_database_connection).has(QueryFeatureMutation('reference:5061:G:A'))
    assert 1 == len(query_result)
