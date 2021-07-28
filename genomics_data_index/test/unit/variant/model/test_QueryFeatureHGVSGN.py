import pytest

from genomics_data_index.storage.model.QueryFeatureHGVSGN import QueryFeatureHGVSGN


def test_create_from_id():
    f = QueryFeatureHGVSGN.create_from_id('hgvs_gn:NC_011083:murF:c.497C>A')
    assert 'hgvs_gn:NC_011083:murF:c.497C>A' == f.id
    assert 'NC_011083' == f.sequence
    assert 'NC_011083' == f.scope
    assert 'murF' == f.gene
    assert 'c.497C>A' == f.variant
    assert f.has_id()
    assert f.has_gene()
    assert f.is_nucleotide()
    assert not f.is_protein()

    f = QueryFeatureHGVSGN.create_from_id('hgvs_gn:NC_011083:murF:p.Ala166Glu')
    assert 'hgvs_gn:NC_011083:murF:p.Ala166Glu' == f.id
    assert 'NC_011083' == f.sequence
    assert 'NC_011083' == f.scope
    assert 'murF' == f.gene
    assert 'p.Ala166Glu' == f.variant
    assert f.has_id()
    assert f.has_gene()
    assert not f.is_nucleotide()
    assert f.is_protein()

    # Test to make sure error if prefix is incorrect
    with pytest.raises(Exception) as execinfo:
        QueryFeatureHGVSGN.create_from_id('hgvs:NC_011083:SEHA_RS01180:c.497C>A')

    assert 'hgvs_id=[hgvs:NC_011083:SEHA_RS01180:c.497C>A] must start with [hgvs_gn:]' in str(execinfo)
