from storage.variant.SampleSet import SampleSet


def test_create_sample_set_from_list():
    sample_set = SampleSet(sample_ids=[1])

    assert 1 in sample_set
    assert 2 not in sample_set
    assert 1 == len(sample_set)
    assert {1} == set(sample_set)


def test_create_sample_set_from_list_multiple():
    sample_set = SampleSet(sample_ids=[1,10])

    assert 1 in sample_set
    assert 2 not in sample_set
    assert 10 in sample_set
    assert 11 not in sample_set
    assert 2 == len(sample_set)
    assert {1,10} == set(sample_set)


def test_create_sample_set_from_list_empty():
    sample_set = SampleSet(sample_ids=[])

    assert 1 not in sample_set
    assert 0 == len(sample_set)
    assert set() == set(sample_set)


def test_serialize_deserialize():
    sample_set = SampleSet(sample_ids=[1,3,10])

    assert {1,3,10} == set(sample_set)

    sample_set_bytes = sample_set.get_bytes()
    assert isinstance(sample_set_bytes, bytes)

    sample_set_deserialize = SampleSet.from_bytes(sample_set_bytes)
    assert {1,3,10} == set(sample_set_deserialize)


def test_intersect_sample_set():
    sample_set1 = SampleSet(sample_ids=[1, 3, 10])
    sample_set2 = SampleSet(sample_ids=[3, 10, 20])

    assert {3, 10} == set(sample_set1.intersection(sample_set2))


def test_intersect_python_set():
    sample_set1 = SampleSet(sample_ids=[1, 3, 10])

    assert {3, 10} == set(sample_set1.intersection({3, 10, 20}))
