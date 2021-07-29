import pytest

from genomics_data_index.storage.model.QueryFeatureMLST import QueryFeatureMLST


def test_create_from_id():
    f = QueryFeatureMLST.create_from_id('mlst:ecoli:adk:100')
    assert 'mlst:ecoli:adk:100' == f.id
    assert 'ecoli' == f.scheme
    assert 'ecoli' == f.scope
    assert 'adk' == f.locus
    assert '100' == f.allele
    assert not f.is_unknown()

    f = QueryFeatureMLST.create_from_id('mlst:ecoli:adk:?')
    assert 'mlst:ecoli:adk:?' == f.id
    assert 'ecoli' == f.scheme
    assert 'ecoli' == f.scope
    assert 'adk' == f.locus
    assert '?' == f.allele
    assert f.is_unknown()

    # Test to make sure error if prefix is incorrect
    with pytest.raises(Exception) as execinfo:
        QueryFeatureMLST.create_from_id('ecoli:adk:100')

    assert 'mlst_id=[ecoli:adk:100] must start with [mlst:]' in str(execinfo)