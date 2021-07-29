import pytest

from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.db import NucleotideVariantsSamples, MLSTAllelesSamples


def test_update_sample_ids_both_overlap():
    s1 = SampleSet([1, 2])
    v1 = MLSTAllelesSamples(sla='ecoli:abc:1', sample_ids=s1)

    s2 = SampleSet([2, 3])
    v2 = MLSTAllelesSamples(sla='ecoli:abc:1', sample_ids=s2)

    v1.update_sample_ids(v2)

    # v1 should only have sample_ids updated
    assert v1.id == 'ecoli:abc:1'
    assert v1.sla == 'ecoli:abc:1'
    assert set(v1.sample_ids) == {1, 2, 3}

    # v2 should not be changed
    assert v2.id == 'ecoli:abc:1'
    assert set(v2.sample_ids) == {2, 3}


def test_update_sample_ids_disjoint():
    s1 = SampleSet([1, 2])
    v1 = MLSTAllelesSamples(sla='ecoli:abc:1', sample_ids=s1)

    s2 = SampleSet([4, 10])
    v2 = MLSTAllelesSamples(sla='ecoli:abc:1', sample_ids=s2)

    v1.update_sample_ids(v2)

    # v1 should only have sample_ids updated
    assert v1.id == 'ecoli:abc:1'
    assert v1.sla == 'ecoli:abc:1'
    assert set(v1.sample_ids) == {1, 2, 4, 10}

    # v2 should not be changed
    assert v2.id == 'ecoli:abc:1'
    assert set(v2.sample_ids) == {4, 10}


def test_update_sample_ids_feature_mismatch():
    s1 = SampleSet([1, 2])
    v1 = MLSTAllelesSamples(sla='ecoli:abc:1', sample_ids=s1)

    s2 = SampleSet([2, 3])
    v2 = NucleotideVariantsSamples(spdi='ref:10:1:A', var_type='SNP', sample_ids=s2)

    with pytest.raises(Exception) as execinfo:
        v1.update_sample_ids(v2)

    assert 'Cannot merge other' in str(execinfo.value)
    assert 'since it is not of type' in str(execinfo.value)


def test_update_sample_ids_feature_id_mismatch():
    s1 = SampleSet([1, 2])
    v1 = MLSTAllelesSamples(sla='ecoli:abc:1', sample_ids=s1)

    s2 = SampleSet([2, 3])
    v2 = MLSTAllelesSamples(sla='ecoli-n:abc:1', sample_ids=s2)

    with pytest.raises(Exception) as execinfo:
        v1.update_sample_ids(v2)

    assert 'Cannot merge other' in str(execinfo.value)
    assert 'since identifiers are not equal' in str(execinfo.value)


def test_query_id():
    s1 = SampleSet([1, 2])
    v1 = MLSTAllelesSamples(sla='ecoli:abc:1', sample_ids=s1)

    assert v1.query_id == 'mlst:ecoli:abc:1'
