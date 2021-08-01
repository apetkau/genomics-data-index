import pytest

from genomics_data_index.storage.model.QueryFeatureHGVSGN import QueryFeatureHGVSGN


def test_create_from_id():
    f = QueryFeatureHGVSGN.create_from_id('hgvs_gn:NC_011083:murF:c.497C>A')
    assert 'hgvs_gn:NC_011083:murF:c.497C>A' == f.id
    assert 'NC_011083:murF:c.497C>A' == f.id_no_prefix
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
    assert 'NC_011083:murF:p.Ala166Glu' == f.id_no_prefix
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


def test_create():
    # Test create where hgvs sequence id begins with 'c.'
    f = QueryFeatureHGVSGN.create(sequence_name='NC_011083', gene='murF', variant='c.497C>A')
    assert 'hgvs_gn:NC_011083:murF:c.497C>A' == f.id
    assert 'NC_011083:murF:c.497C>A' == f.id_no_prefix
    assert 'NC_011083' == f.sequence
    assert 'NC_011083' == f.scope
    assert 'murF' == f.gene
    assert 'c.497C>A' == f.variant
    assert f.has_id()
    assert f.has_gene()
    assert f.is_nucleotide()
    assert not f.is_protein()

    # Test create where hgvs sequence id begins with 'n.'
    f = QueryFeatureHGVSGN.create(sequence_name='NC_011083', gene='gene', variant='n.1031571T>C')
    assert 'hgvs_gn:NC_011083:n.1031571T>C' == f.id
    assert 'NC_011083:n.1031571T>C' == f.id_no_prefix
    assert 'NC_011083' == f.sequence
    assert 'NC_011083' == f.scope
    assert f.gene is None
    assert 'n.1031571T>C' == f.variant
    assert f.has_id()
    assert not f.has_gene()
    assert f.is_nucleotide()
    assert not f.is_protein()
