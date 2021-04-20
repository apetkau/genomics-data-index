from storage.variant.SampleSet import SampleSet


def test_create_sample_set_from_list():
    sample_set = SampleSet(sample_ids=[1])

    assert 1 in sample_set
    assert 2 not in sample_set
    assert 1 == len(sample_set)
    assert {1} == set(sample_set)
    assert not sample_set.is_empty()


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

    intersection = sample_set1.intersection(sample_set2)
    assert isinstance(intersection, SampleSet)
    assert {3, 10} == set(intersection)


def test_complement_sample_set():
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
