import pytest

from genomics_data_index.storage.model.NucleotideMutationTranslater import NucleotideMutationTranslater
from genomics_data_index.storage.model.QueryFeatureMutationSPDI import QueryFeatureMutationSPDI


def test_to_spdi_position():
    assert 'seq:1:1:T' == NucleotideMutationTranslater.to_spdi('seq', 1, 1, 'T')
    assert 'seq:2:1:T' == NucleotideMutationTranslater.to_spdi('seq', 2, 1, 'T')
    assert 'seq:0:1:T' == NucleotideMutationTranslater.to_spdi('seq', 0, 1, 'T')

    with pytest.raises(Exception) as execinfo:
        NucleotideMutationTranslater.to_spdi('seq', -1, 1, 'T')
    assert 'Position must be non-negative' in str(execinfo.value)


def test_to_spdi_deletion():
    assert 'seq:1:1:T' == NucleotideMutationTranslater.to_spdi('seq', 1, 1, 'T')
    assert 'seq:1:2:T' == NucleotideMutationTranslater.to_spdi('seq', 1, 2, 'T')
    assert 'seq:1:1:T' == NucleotideMutationTranslater.to_spdi('seq', 1, 'A', 'T')
    assert 'seq:1:2:T' == NucleotideMutationTranslater.to_spdi('seq', 1, 'AT', 'T')
    assert 'seq:1:0:T' == NucleotideMutationTranslater.to_spdi('seq', 1, '', 'T')
    assert 'seq:1:1:T' == NucleotideMutationTranslater.to_spdi('seq', 1, '1', 'T')


def test_to_spdi_deletion_no_convert_deletion():
    assert 'seq:1:1:T' == NucleotideMutationTranslater.to_spdi('seq', 1, 1, 'T', convert_deletion=False)
    assert 'seq:1:2:T' == NucleotideMutationTranslater.to_spdi('seq', 1, 2, 'T', convert_deletion=False)
    assert 'seq:1:A:T' == NucleotideMutationTranslater.to_spdi('seq', 1, 'A', 'T', convert_deletion=False)
    assert 'seq:1:AT:T' == NucleotideMutationTranslater.to_spdi('seq', 1, 'AT', 'T', convert_deletion=False)
    assert 'seq:1::T' == NucleotideMutationTranslater.to_spdi('seq', 1, '', 'T', convert_deletion=False)
    assert 'seq:1:1:T' == NucleotideMutationTranslater.to_spdi('seq', 1, '1', 'T', convert_deletion=False)


def test_to_spdi_alt():
    assert 'seq:1:1:T' == NucleotideMutationTranslater.to_spdi('seq', 1, 1, 'T')
    assert 'seq:1:1:AT' == NucleotideMutationTranslater.to_spdi('seq', 1, 1, 'AT')
    assert 'seq:1:1:' == NucleotideMutationTranslater.to_spdi('seq', 1, 1, '')


def test_from_spdi_position():
    assert ('seq', 1, 1, 'T') == NucleotideMutationTranslater.from_spdi('seq:1:A:T')
    assert ('seq', 2, 1, 'T') == NucleotideMutationTranslater.from_spdi('seq:2:A:T')
    assert ('seq', 0, 1, 'T') == NucleotideMutationTranslater.from_spdi('seq:0:A:T')

    with pytest.raises(Exception) as execinfo:
        NucleotideMutationTranslater.from_spdi('seq:-1:A:T')
    assert 'Position must be non-negative' in str(execinfo.value)


def test_from_spdi_deletion():
    assert ('seq', 1, 1, 'T') == NucleotideMutationTranslater.from_spdi('seq:1:A:T')
    assert ('seq', 1, 1, 'T') == NucleotideMutationTranslater.from_spdi('seq:1:1:T')
    assert ('seq', 1, 2, 'T') == NucleotideMutationTranslater.from_spdi('seq:1:2:T')
    assert ('seq', 1, 2, 'T') == NucleotideMutationTranslater.from_spdi('seq:1:AT:T')
    assert ('seq', 1, 0, 'T') == NucleotideMutationTranslater.from_spdi('seq:1:0:T')
    assert ('seq', 1, 0, 'T') == NucleotideMutationTranslater.from_spdi('seq:1::T')


def test_from_spdi_deletion_no_convert_deletion():
    assert ('seq', 1, 'A', 'T') == NucleotideMutationTranslater.from_spdi('seq:1:A:T', convert_deletion=False)
    assert ('seq', 1, 1, 'T') == NucleotideMutationTranslater.from_spdi('seq:1:1:T', convert_deletion=False)
    assert ('seq', 1, 2, 'T') == NucleotideMutationTranslater.from_spdi('seq:1:2:T', convert_deletion=False)
    assert ('seq', 1, 'AT', 'T') == NucleotideMutationTranslater.from_spdi('seq:1:AT:T', convert_deletion=False)
    assert ('seq', 1, 0, 'T') == NucleotideMutationTranslater.from_spdi('seq:1:0:T', convert_deletion=False)

    with pytest.raises(Exception) as execinfo:
        NucleotideMutationTranslater.from_spdi('seq:1::T', convert_deletion=False)
    assert 'deletion=[] but convert_deletion is False' in str(execinfo.value)


def test_from_spdi_insertion():
    assert ('seq', 1, 1, 'T') == NucleotideMutationTranslater.from_spdi('seq:1:A:T')
    assert ('seq', 2, 1, 'TA') == NucleotideMutationTranslater.from_spdi('seq:2:A:TA')
    assert ('seq', 0, 1, '') == NucleotideMutationTranslater.from_spdi('seq:0:A:')


def test_to_db_feature():
    assert 'seq:1:1:T' == NucleotideMutationTranslater.to_db_feature(QueryFeatureMutationSPDI('seq:1:A:T')).id
    assert 'seq:1:1:T' == NucleotideMutationTranslater.to_db_feature(QueryFeatureMutationSPDI('seq:1:1:T')).id
    assert 'seq:1:2:T' == NucleotideMutationTranslater.to_db_feature(QueryFeatureMutationSPDI('seq:1:AT:T')).id
    assert 'seq:1:2:T' == NucleotideMutationTranslater.to_db_feature(QueryFeatureMutationSPDI('seq:1:2:T')).id
    assert 'seq:10:2:TA' == NucleotideMutationTranslater.to_db_feature(QueryFeatureMutationSPDI('seq:10:AT:TA')).id


def test_to_hgvs_id():
    assert 'hgvs:NC_011083:SEHA_RS01180:p.Ala166Glu' == NucleotideMutationTranslater.to_hgvs_id('NC_011083',
                                                                                                'SEHA_RS01180',
                                                                                                'p.Ala166Glu')
    assert 'hgvs:NC_011083:n.1031571T>C' == NucleotideMutationTranslater.to_hgvs_id('NC_011083', None, 'n.1031571T>C')
    assert 'hgvs:NC_011083:n.1031571T>C' == NucleotideMutationTranslater.to_hgvs_id('NC_011083', 'gene_name',
                                                                                    'n.1031571T>C')
