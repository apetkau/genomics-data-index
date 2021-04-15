import pytest

from storage.variant.model.QueryFeatureMutation import QueryFeatureMutation


def test_to_spdi_position():
    assert 'seq:1:1:T' == QueryFeatureMutation.to_spdi('seq', 1, 1, 'T')
    assert 'seq:2:1:T' == QueryFeatureMutation.to_spdi('seq', 2, 1, 'T')
    assert 'seq:0:1:T' == QueryFeatureMutation.to_spdi('seq', 0, 1, 'T')

    with pytest.raises(Exception) as execinfo:
        QueryFeatureMutation.to_spdi('seq', -1, 1, 'T')
    assert 'Position must be non-negative' in str(execinfo.value)


def test_to_spdi_deletion():
    assert 'seq:1:1:T' == QueryFeatureMutation.to_spdi('seq', 1, 1, 'T')
    assert 'seq:1:2:T' == QueryFeatureMutation.to_spdi('seq', 1, 2, 'T')
    assert 'seq:1:1:T' == QueryFeatureMutation.to_spdi('seq', 1, 'A', 'T')
    assert 'seq:1:2:T' == QueryFeatureMutation.to_spdi('seq', 1, 'AT', 'T')
    assert 'seq:1:0:T' == QueryFeatureMutation.to_spdi('seq', 1, '', 'T')
    assert 'seq:1:1:T' == QueryFeatureMutation.to_spdi('seq', 1, '1', 'T')


def test_to_spdi_alt():
    assert 'seq:1:1:T' == QueryFeatureMutation.to_spdi('seq', 1, 1, 'T')
    assert 'seq:1:1:AT' == QueryFeatureMutation.to_spdi('seq', 1, 1, 'AT')
    assert 'seq:1:1:' == QueryFeatureMutation.to_spdi('seq', 1, 1, '')


def test_from_spdi_position():
    assert ('seq', 1, 1, 'T') == QueryFeatureMutation.from_spdi('seq:1:A:T')
    assert ('seq', 2, 1, 'T') == QueryFeatureMutation.from_spdi('seq:2:A:T')
    assert ('seq', 0, 1, 'T') == QueryFeatureMutation.from_spdi('seq:0:A:T')

    with pytest.raises(Exception) as execinfo:
        QueryFeatureMutation.from_spdi('seq:-1:A:T')
    assert 'Position must be non-negative' in str(execinfo.value)


def test_from_spdi_deletion():
    assert ('seq', 1, 1, 'T') == QueryFeatureMutation.from_spdi('seq:1:A:T')
    assert ('seq', 1, 1, 'T') == QueryFeatureMutation.from_spdi('seq:1:1:T')
    assert ('seq', 1, 2, 'T') == QueryFeatureMutation.from_spdi('seq:1:2:T')
    assert ('seq', 1, 2, 'T') == QueryFeatureMutation.from_spdi('seq:1:AT:T')
    assert ('seq', 1, 0, 'T') == QueryFeatureMutation.from_spdi('seq:1:0:T')
    assert ('seq', 1, 0, 'T') == QueryFeatureMutation.from_spdi('seq:1::T')


def test_from_spdi_insertion():
    assert ('seq', 1, 1, 'T') == QueryFeatureMutation.from_spdi('seq:1:A:T')
    assert ('seq', 2, 1, 'TA') == QueryFeatureMutation.from_spdi('seq:2:A:TA')
    assert ('seq', 0, 1, '') == QueryFeatureMutation.from_spdi('seq:0:A:')
