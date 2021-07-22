import pytest

from genomics_data_index.storage.util.ListSliceIter import ListSliceIter


def test_slice_list():
    # Case 3 elements
    data = [1, 2, 3]
    assert [[1], [2], [3]] == list(ListSliceIter(data, slice_size=1).islice())
    assert [[1, 2], [3]] == list(ListSliceIter(data, slice_size=2).islice())
    assert [[1, 2, 3]] == list(ListSliceIter(data, slice_size=3).islice())
    assert [[1, 2, 3]] == list(ListSliceIter(data, slice_size=4).islice())

    # Trying slice < 1
    with pytest.raises(Exception) as execinfo:
        ListSliceIter(data, slice_size=0)
    assert 'must be a positive integer' in str(execinfo.value)

    # Case empty data
    data = []
    assert [] == list(ListSliceIter(data, slice_size=1).islice())
    assert [] == list(ListSliceIter(data, slice_size=2).islice())

    # Case 5 elements and strings
    data = ['A', 'B', 'C', 'D', 'E']
    assert [['A'], ['B'], ['C'], ['D'], ['E']] == list(ListSliceIter(data, slice_size=1).islice())
    assert [['A', 'B'], ['C', 'D'], ['E']] == list(ListSliceIter(data, slice_size=2).islice())
    assert [['A', 'B', 'C'], ['D', 'E']] == list(ListSliceIter(data, slice_size=3).islice())
    assert [['A', 'B', 'C', 'D'], ['E']] == list(ListSliceIter(data, slice_size=4).islice())
    assert [['A', 'B', 'C', 'D', 'E']] == list(ListSliceIter(data, slice_size=5).islice())
    assert [['A', 'B', 'C', 'D', 'E']] == list(ListSliceIter(data, slice_size=6).islice())
