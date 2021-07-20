from typing import cast

from genomics_data_index.storage.model.QueryFeatureHGVS import QueryFeatureHGVS


def test_create_from_id():
    f = QueryFeatureHGVS.create_from_id('hgvs:NC_011083:SEHA_RS01180:c.497C>A')
    assert 'hgvs:NC_011083:SEHA_RS01180:c.497C>A' == f.id
    assert 'NC_011083' == f.sequence
    assert 'NC_011083' == f.scope
    assert 'SEHA_RS01180' == f.gene
    assert 'c.497C>A' == f.variant
    assert f.has_id()
    assert f.has_gene()
    assert f.is_nucleotide()
    assert not f.is_protein()

    f = QueryFeatureHGVS.create_from_id('hgvs:NC_011083:SEHA_RS01180:p.Ala166Glu')
    assert 'hgvs:NC_011083:SEHA_RS01180:p.Ala166Glu' == f.id
    assert 'NC_011083' == f.sequence
    assert 'NC_011083' == f.scope
    assert 'SEHA_RS01180' == f.gene
    assert 'p.Ala166Glu' == f.variant
    assert f.has_id()
    assert f.has_gene()
    assert not f.is_nucleotide()
    assert f.is_protein()


def test_create():
    f = QueryFeatureHGVS.create('NC_011083', 'SEHA_RS01180', 'c.497C>A')
    assert 'hgvs:NC_011083:SEHA_RS01180:c.497C>A' == f.id
    assert 'NC_011083' == f.sequence
    assert 'NC_011083' == f.scope
    assert 'SEHA_RS01180' == f.gene
    assert 'c.497C>A' == f.variant
    assert f.has_id()
    assert f.has_gene()
    assert f.is_nucleotide()
    assert not f.is_protein()

    f = QueryFeatureHGVS.create('NC_011083', 'SEHA_RS01180', 'p.Ala166Glu')
    assert 'hgvs:NC_011083:SEHA_RS01180:p.Ala166Glu' == f.id
    assert 'NC_011083' == f.sequence
    assert 'NC_011083' == f.scope
    assert 'SEHA_RS01180' == f.gene
    assert 'p.Ala166Glu' == f.variant
    assert f.has_id()
    assert f.has_gene()
    assert not f.is_nucleotide()
    assert f.is_protein()

    f = QueryFeatureHGVS.create('NC_011083', None, 'n.1031571T>C')
    assert 'hgvs:NC_011083:n.1031571T>C' == f.id
    assert 'NC_011083' == f.sequence
    assert 'NC_011083' == f.scope
    assert f.gene is None
    assert 'n.1031571T>C' == f.variant
    assert f.has_id()
    assert not f.has_gene()
    assert f.is_nucleotide()
    assert not f.is_protein()

    # Don't include gene id if variant starts with 'n.' (coordinates with respect to sequence
    f = QueryFeatureHGVS.create('NC_011083', 'gene', 'n.1031571T>C')
    assert 'hgvs:NC_011083:n.1031571T>C' == f.id
    assert 'NC_011083' == f.sequence
    assert 'NC_011083' == f.scope
    assert f.gene is None
    assert 'n.1031571T>C' == f.variant
    assert f.has_id()
    assert not f.has_gene()
    assert f.is_nucleotide()
    assert not f.is_protein()

    # None
    f = QueryFeatureHGVS.create('NC_011083', 'gene', None)
    assert f.id is None
    assert 'NC_011083' == f.sequence
    assert 'NC_011083' == f.scope
    assert f.gene == 'gene'
    assert f.variant is None
    assert f.has_gene()
    assert not f.has_id()
    assert not f.is_nucleotide()
    assert not f.is_protein()

    # None2
    f = QueryFeatureHGVS.create('NC_011083', None, None)
    assert f.id is None
    assert 'NC_011083' == f.sequence
    assert 'NC_011083' == f.scope
    assert f.gene is None
    assert f.variant is None
    assert not f.has_gene()
    assert not f.has_id()
    assert not f.is_nucleotide()
    assert not f.is_protein()


def test_unknown():
    f = QueryFeatureHGVS.create_from_id('hgvs:NC_011083:SEHA_RS01180:c.497C>A')
    assert not f.is_unknown()
    assert f.id == 'hgvs:NC_011083:SEHA_RS01180:c.497C>A'
    assert f.id_minus_unknown == 'hgvs:NC_011083:SEHA_RS01180:c.497C>A'

    u = f.to_unknown()
    assert u.is_unknown()
    assert u.id == 'hgvs:NC_011083:SEHA_RS01180:c.497C>A:unknown'
    assert cast(QueryFeatureHGVS, u).id_minus_unknown == 'hgvs:NC_011083:SEHA_RS01180:c.497C>A'
    assert not f.is_unknown()

    ue = f.to_unknown_explode()
    assert 1 == len(ue)
    assert ue[0].is_unknown()
    assert ue[0].id == 'hgvs:NC_011083:SEHA_RS01180:c.497C>A:unknown'
    assert cast(QueryFeatureHGVS, ue[0]).id_minus_unknown == 'hgvs:NC_011083:SEHA_RS01180:c.497C>A'

    uu = u.to_unknown()
    assert uu.is_unknown()
    assert uu.id == 'hgvs:NC_011083:SEHA_RS01180:c.497C>A:unknown'
    assert cast(QueryFeatureHGVS, uu).id_minus_unknown == 'hgvs:NC_011083:SEHA_RS01180:c.497C>A'
