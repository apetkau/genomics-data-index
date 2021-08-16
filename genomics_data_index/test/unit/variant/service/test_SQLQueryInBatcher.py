from genomics_data_index.storage.service import SQLQueryInBatcher


def test_sql_query_in_batcher():
    in_data = ['A', 'B', 'C', 'D', 'E']

    # Test batch size 1
    batcher = SQLQueryInBatcher(in_data=in_data, batch_size=1)
    results = batcher.process(lambda in_batch: {x: True for x in in_batch})
    assert isinstance(results, dict)
    assert 5 == len(results)
    assert {'A', 'B', 'C', 'D', 'E'} == set(results.keys())

    # Test batch size 2
    batcher = SQLQueryInBatcher(in_data=in_data, batch_size=2)
    results = batcher.process(lambda in_batch: {x: True for x in in_batch})
    assert isinstance(results, dict)
    assert 5 == len(results)
    assert {'A', 'B', 'C', 'D', 'E'} == set(results.keys())

    # Test batch size 5
    batcher = SQLQueryInBatcher(in_data=in_data, batch_size=5)
    results = batcher.process(lambda in_batch: {x: True for x in in_batch})
    assert isinstance(results, dict)
    assert 5 == len(results)
    assert {'A', 'B', 'C', 'D', 'E'} == set(results.keys())

    # Test batch size 6
    batcher = SQLQueryInBatcher(in_data=in_data, batch_size=6)
    results = batcher.process(lambda in_batch: {x: True for x in in_batch})
    assert isinstance(results, dict)
    assert 5 == len(results)
    assert {'A', 'B', 'C', 'D', 'E'} == set(results.keys())
