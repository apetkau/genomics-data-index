import pytest

from genomics_data_index.storage.model.QueryFeatureFactory import QueryFeatureFactory
from genomics_data_index.storage.model.QueryFeatureHGVS import QueryFeatureHGVS


query_feature_factory = QueryFeatureFactory.instance()


def test_create_hgvs():
    f = query_feature_factory.create_feature('hgvs:NC_011083:SEHA_RS01180:c.497C>A')
    assert isinstance(f, QueryFeatureHGVS)
    assert 'hgvs:NC_011083:SEHA_RS01180:c.497C>A' == f.id
    assert 'NC_011083' == f.sequence
    assert 'NC_011083' == f.scope
    assert 'SEHA_RS01180' == f.gene
    assert 'c.497C>A' == f.variant
    assert f.has_id()
    assert f.has_gene()
    assert f.is_nucleotide()
    assert not f.is_protein()

    f = query_feature_factory.create_feature('hgvs:NC_011083:SEHA_RS01180:p.Ala166Glu')
    assert isinstance(f, QueryFeatureHGVS)
    assert 'hgvs:NC_011083:SEHA_RS01180:p.Ala166Glu' == f.id
    assert 'NC_011083' == f.sequence
    assert 'NC_011083' == f.scope
    assert 'SEHA_RS01180' == f.gene
    assert 'p.Ala166Glu' == f.variant
    assert f.has_id()
    assert f.has_gene()
    assert not f.is_nucleotide()
    assert f.is_protein()


def test_create_unknown():
    with pytest.raises(Exception) as execinfo:
        query_feature_factory.create_feature('unknown:feature')

    assert 'Unknown feature type for feature=[unknown:feature]' in str(execinfo)
