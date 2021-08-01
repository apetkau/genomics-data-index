import math

import pytest

from genomics_data_index.storage.SampleSet import SampleSet


def test_create_sample_set_from_list():
    sample_set = SampleSet(sample_ids=[1])

    assert 1 in sample_set
    assert 2 not in sample_set
    assert 1 == len(sample_set)
    assert {1} == set(sample_set)
    assert not sample_set.is_empty()


def test_create_sample_set_from_set():
    assert {1} == set(SampleSet(sample_ids={1}))
    assert {1, 10} == set(SampleSet(sample_ids={1, 10}))
    assert set() == set(SampleSet(sample_ids=set()))


def test_create_sample_set_from_list_multiple():
    sample_set = SampleSet(sample_ids=[1, 10])

    assert 1 in sample_set
    assert 2 not in sample_set
    assert 10 in sample_set
    assert 11 not in sample_set
    assert 2 == len(sample_set)
    assert {1, 10} == set(sample_set)


def test_create_sample_set_from_list_empty():
    sample_set = SampleSet(sample_ids=[])

    assert 1 not in sample_set
    assert 0 == len(sample_set)
    assert set() == set(sample_set)
    assert sample_set.is_empty()


def test_serialize_deserialize():
    sample_set = SampleSet(sample_ids=[1, 3, 10])

    assert {1, 3, 10} == set(sample_set)

    sample_set_bytes = sample_set.get_bytes()
    assert isinstance(sample_set_bytes, bytes)

    sample_set_deserialize = SampleSet.from_bytes(sample_set_bytes)
    assert {1, 3, 10} == set(sample_set_deserialize)


def test_intersect_sample_set():
    sample_set1 = SampleSet(sample_ids=[1, 3, 10])
    sample_set2 = SampleSet(sample_ids=[3, 10, 20])
    sample_set_empty = SampleSet(sample_ids=[])
    sample_set_non_overlap = SampleSet(sample_ids=[50, 100])

    intersection = sample_set1.intersection(sample_set2)
    assert isinstance(intersection, SampleSet)
    assert {3, 10} == set(intersection)

    intersection = sample_set2.intersection(sample_set1)
    assert isinstance(intersection, SampleSet)
    assert {3, 10} == set(intersection)

    intersection = sample_set1.intersection(sample_set_empty)
    assert isinstance(intersection, SampleSet)
    assert set() == set(intersection)

    intersection = sample_set1.intersection(sample_set_non_overlap)
    assert isinstance(intersection, SampleSet)
    assert set() == set(intersection)


def test_union_sample_set():
    sample_set1 = SampleSet(sample_ids=[1, 3, 10])
    sample_set2 = SampleSet(sample_ids=[3, 10, 20])
    sample_set_empty = SampleSet(sample_ids=[])
    sample_set_non_overlap = SampleSet(sample_ids=[50, 100])

    union = sample_set1.union(sample_set2)
    assert isinstance(union, SampleSet)
    assert {1, 3, 10, 20} == set(union)

    union = sample_set2.union(sample_set1)
    assert isinstance(union, SampleSet)
    assert {1, 3, 10, 20} == set(union)

    union = sample_set1.union(sample_set_empty)
    assert isinstance(union, SampleSet)
    assert {1, 3, 10} == set(union)

    union = sample_set1.union(sample_set_non_overlap)
    assert isinstance(union, SampleSet)
    assert {1, 3, 10, 50, 100} == set(union)


def test_jaccard_index_sample_set():
    sample_set1 = SampleSet(sample_ids=[1, 3, 10])
    sample_set2 = SampleSet(sample_ids=[3, 10, 20])
    sample_set_empty = SampleSet(sample_ids=[])
    sample_set_non_overlap = SampleSet(sample_ids=[50, 100])
    sample_set_3 = SampleSet(sample_ids=[1, 3, 10, 20, 30])

    assert math.isclose(1.0, sample_set1.jaccard_index(sample_set1), rel_tol=1e-3)

    assert math.isclose(0.5, sample_set1.jaccard_index(sample_set2), rel_tol=1e-3)
    assert math.isclose(0.5, sample_set2.jaccard_index(sample_set1), rel_tol=1e-3)

    assert math.isclose(0, sample_set1.jaccard_index(sample_set_empty), rel_tol=1e-3)
    assert math.isclose(0, sample_set_empty.jaccard_index(sample_set1), rel_tol=1e-3)

    assert math.isclose(0, sample_set1.jaccard_index(sample_set_non_overlap), rel_tol=1e-3)
    assert math.isclose(0, sample_set_non_overlap.jaccard_index(sample_set1), rel_tol=1e-3)

    assert math.isclose(0.6, sample_set1.jaccard_index(sample_set_3), rel_tol=1e-3)
    assert math.isclose(0.6, sample_set_3.jaccard_index(sample_set1), rel_tol=1e-3)


def test_jaccard_index_sample_ids_set():
    sample_set1 = SampleSet(sample_ids=[1, 3, 10])
    sample_set2 = {3, 10, 20}
    sample_set_empty = set()
    sample_set_non_overlap = {50, 100}
    sample_set_3 = {1, 3, 10, 20, 30}

    assert math.isclose(1.0, sample_set1.jaccard_index(sample_set1), rel_tol=1e-3)
    assert math.isclose(0.5, sample_set1.jaccard_index(sample_set2), rel_tol=1e-3)
    assert math.isclose(0, sample_set1.jaccard_index(sample_set_empty), rel_tol=1e-3)
    assert math.isclose(0, sample_set1.jaccard_index(sample_set_non_overlap), rel_tol=1e-3)
    assert math.isclose(0.6, sample_set1.jaccard_index(sample_set_3), rel_tol=1e-3)


def test_minus_sample_set():
    set1 = SampleSet(sample_ids=[1, 3, 10])
    set2 = SampleSet(sample_ids=[3, 10, 20])
    set3 = SampleSet(sample_ids=[40, 60, 80, 100])

    assert {1} == set(set1.minus(set2))
    assert {20} == set(set2.minus(set1))
    assert set() == set(set1.minus(set1))
    assert {1, 3, 10} == set(set1.minus(set3))
    assert {40, 60, 80, 100} == set(set3.minus(set1))
    assert {1, 3, 10} == set(set1.minus(SampleSet.create_empty()))


def test_intersect_python_set():
    sample_set1 = SampleSet(sample_ids=[1, 3, 10])

    assert {3, 10} == set(sample_set1.intersection({3, 10, 20}))


def test_union_python_set():
    sample_set1 = SampleSet(sample_ids=[1, 3, 10])

    assert {1, 3, 10, 20} == set(sample_set1.union({3, 10, 20}))
    assert {1, 3, 10} == set(sample_set1.union(set()))
    assert {1, 3, 10, 50, 100} == set(sample_set1.union({50, 100}))


def test_create_empty_sample_set():
    empty_set = SampleSet.create_empty()
    assert len(empty_set) == 0
    assert set() == set(empty_set)
    assert empty_set.is_empty()


def test_create_all_sample_set():
    all_set = SampleSet.create_all()
    assert len(all_set) == 2 ** 32
    assert 0 in all_set
    assert 1 in all_set
    assert 2 ** 32 - 1 in all_set
    assert 2 ** 32 not in all_set

    other_set = SampleSet(sample_ids=[1, 3, 10])

    assert all_set.intersection(other_set) == other_set
    assert other_set.intersection(all_set) == other_set
    assert isinstance(all_set.intersection(other_set), SampleSet)
    assert isinstance(other_set.intersection(all_set), SampleSet)


def test_minus_2():
    sample_set1 = SampleSet(sample_ids=[1, 3, 10])
    sample_set1_identical = SampleSet(sample_ids=[1, 3, 10])
    sample_set2 = SampleSet(sample_ids=[3, 10, 20])
    sample_set3 = SampleSet(sample_ids=[3])
    sample_set4 = SampleSet(sample_ids=[1])
    sample_set_empty = SampleSet(sample_ids=[])
    sample_set_non_overlap = SampleSet(sample_ids=[50, 100])

    subtract = sample_set1.minus(sample_set2)
    assert isinstance(subtract, SampleSet)
    assert {1} == set(subtract)

    subtract = sample_set2.minus(sample_set1)
    assert isinstance(subtract, SampleSet)
    assert {20} == set(subtract)

    subtract = sample_set1.minus(sample_set_empty)
    assert isinstance(subtract, SampleSet)
    assert {1, 3, 10} == set(subtract)

    subtract = sample_set_empty.minus(sample_set1)
    assert isinstance(subtract, SampleSet)
    assert set() == set(subtract)

    subtract = sample_set1.minus(sample_set_non_overlap)
    assert isinstance(subtract, SampleSet)
    assert {1, 3, 10} == set(subtract)

    subtract = sample_set1.minus(sample_set1_identical)
    assert isinstance(subtract, SampleSet)
    assert set() == set(subtract)

    subtract = sample_set1.minus(sample_set3)
    assert isinstance(subtract, SampleSet)
    assert {1, 10} == set(subtract)

    subtract = sample_set1.minus([sample_set3, sample_set4])
    assert isinstance(subtract, SampleSet)
    assert {10} == set(subtract)

    # Pass in a set of integers
    subtract = sample_set1.minus({1, 10})
    assert isinstance(subtract, SampleSet)
    assert {3} == set(subtract)

    # Test case of all set
    all_set = SampleSet.create_all()
    subtract = sample_set1.minus(all_set)
    assert isinstance(subtract, SampleSet)
    assert set() == set(subtract)

    # Test case of all set other
    with pytest.raises(Exception) as execinfo:
        all_set.minus(sample_set1)
    assert 'Cannot subtract anything from AllSampleSet' in str(execinfo.value)
