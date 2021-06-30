from genomics_data_index.storage.model.QueryFeatureHGVS import QueryFeatureHGVS


def test_parse_hgvs_id():
    f = QueryFeatureHGVS('hgvs:NC_011083:murF:c.497C>A')
    assert 'hgvs:NC_011083:murF:c.497C>A' == f.id
    assert 'NC_011083' == f.reference
    assert 'NC_011083' == f.scope
    assert 'murF' == f.gene
    assert 'c.497C>A' == f.mutation
    assert f.has_gene()
    assert f.is_nucleotide()
    assert not f.is_protein()

    f = QueryFeatureHGVS('hgvs:NC_011083:murF:p.Ala166Glu')
    assert 'hgvs:NC_011083:murF:p.Ala166Glu' == f.id
    assert 'NC_011083' == f.reference
    assert 'NC_011083' == f.scope
    assert 'murF' == f.gene
    assert 'p.Ala166Glu' == f.mutation
    assert f.has_gene()
    assert not f.is_nucleotide()
    assert f.is_protein()
