import pytest

from genomics_data_index.storage.model.QueryFeatureFactory import QueryFeatureFactory
from genomics_data_index.storage.model.QueryFeatureHGVS import QueryFeatureHGVS
from genomics_data_index.storage.model.QueryFeatureHGVSGN import QueryFeatureHGVSGN
from genomics_data_index.storage.model.QueryFeatureMLST import QueryFeatureMLST
from genomics_data_index.storage.model.QueryFeatureMutationSPDI import QueryFeatureMutationSPDI

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


def test_create_hgvs_gn():
    f = query_feature_factory.create_feature('hgvs_gn:NC_011083:murF:c.497C>A')
    assert isinstance(f, QueryFeatureHGVSGN)
    assert 'hgvs_gn:NC_011083:murF:c.497C>A' == f.id
    assert 'NC_011083' == f.sequence
    assert 'NC_011083' == f.scope
    assert 'murF' == f.gene
    assert 'c.497C>A' == f.variant
    assert f.has_id()
    assert f.has_gene()
    assert f.is_nucleotide()
    assert not f.is_protein()

    f = query_feature_factory.create_feature('hgvs_gn:NC_011083:murF:p.Ala166Glu')
    assert isinstance(f, QueryFeatureHGVSGN)
    assert 'hgvs_gn:NC_011083:murF:p.Ala166Glu' == f.id
    assert 'NC_011083' == f.sequence
    assert 'NC_011083' == f.scope
    assert 'murF' == f.gene
    assert 'p.Ala166Glu' == f.variant
    assert f.has_id()
    assert f.has_gene()
    assert not f.is_nucleotide()
    assert f.is_protein()


def test_create_spdi():
    # Test SNV
    f = query_feature_factory.create_feature('reference:10:A:C')
    assert isinstance(f, QueryFeatureMutationSPDI)
    assert 'reference:10:A:C' == f.id
    assert 'reference' == f.sequence
    assert 10 == f.position
    assert 'A' == f.deletion
    assert 'C' == f.insertion
    assert 1 == f.deletion_length()
    assert f.has_deletion_sequence()
    assert not f.is_unknown()

    # Test complex
    f = query_feature_factory.create_feature('reference:10:ATT:CCC')
    assert isinstance(f, QueryFeatureMutationSPDI)
    assert 'reference:10:ATT:CCC' == f.id
    assert 'reference' == f.sequence
    assert 10 == f.position
    assert 'ATT' == f.deletion
    assert 'CCC' == f.insertion
    assert 3 == f.deletion_length()
    assert f.has_deletion_sequence()
    assert not f.is_unknown()

    # Test complex with integer deletion
    f = query_feature_factory.create_feature('reference:10:3:CCC')
    assert isinstance(f, QueryFeatureMutationSPDI)
    assert 'reference:10:3:CCC' == f.id
    assert 'reference' == f.sequence
    assert 10 == f.position
    assert 3 == f.deletion
    assert 'CCC' == f.insertion
    assert 3 == f.deletion_length()
    assert not f.has_deletion_sequence()
    assert not f.is_unknown()

    # Test unknown integer deletion
    f = query_feature_factory.create_feature('reference:10:1:?')
    assert isinstance(f, QueryFeatureMutationSPDI)
    assert 'reference:10:1:?' == f.id
    assert 'reference' == f.sequence
    assert 10 == f.position
    assert 1 == f.deletion
    assert '?' == f.insertion
    assert 1 == f.deletion_length()
    assert not f.has_deletion_sequence()
    assert f.is_unknown()

    # Test unknown, deletion sequence
    f = query_feature_factory.create_feature('reference:10:A:?')
    assert isinstance(f, QueryFeatureMutationSPDI)
    assert 'reference:10:A:?' == f.id
    assert 'reference' == f.sequence
    assert 10 == f.position
    assert 'A' == f.deletion
    assert '?' == f.insertion
    assert 1 == f.deletion_length()
    assert f.has_deletion_sequence()
    assert f.is_unknown()


def test_create_mlst():
    f = query_feature_factory.create_feature('mlst:ecoli:adk:100')
    assert isinstance(f, QueryFeatureMLST)
    assert 'mlst:ecoli:adk:100' == f.id
    assert 'ecoli:adk:100' == f.id_no_prefix
    assert 'ecoli' == f.scheme
    assert 'ecoli' == f.scope
    assert 'adk' == f.locus
    assert '100' == f.allele
    assert not f.is_unknown()

    f = query_feature_factory.create_feature('mlst:ecoli:adk:?')
    assert isinstance(f, QueryFeatureMLST)
    assert 'mlst:ecoli:adk:?' == f.id
    assert 'ecoli:adk:?' == f.id_no_prefix
    assert 'ecoli' == f.scheme
    assert 'ecoli' == f.scope
    assert 'adk' == f.locus
    assert '?' == f.allele
    assert f.is_unknown()


def test_create_unknown():
    with pytest.raises(Exception) as execinfo:
        query_feature_factory.create_feature('unknown:feature')

    assert 'Unknown feature type for feature=[unknown:feature]' in str(execinfo)
