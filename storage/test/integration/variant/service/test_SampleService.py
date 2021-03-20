import pytest

from storage.variant.SampleSet import SampleSet


def test_samples_with_variants(sample_service, variation_service):
    samples_with_variants = sample_service.get_samples_with_variants('genome')
    assert {'SampleA', 'SampleB', 'SampleC'} == {sample.name for sample in samples_with_variants}


def test_samples_with_variants_on_sequence(sample_service, variation_service):
    samples_with_variants = sample_service.get_samples_with_variants_on_sequence('reference')
    assert {'SampleA', 'SampleB', 'SampleC'} == {sample.name for sample in samples_with_variants}


def test_samples_associated_with_reference(sample_service, variation_service):
    samples_reference = sample_service.get_samples_associated_with_reference('genome')
    assert {'SampleA', 'SampleB', 'SampleC'} == {sample.name for sample in samples_reference}


def test_find_sample_name_ids(sample_service, variation_service):
    sample_ids_list = sample_service.find_sample_name_ids(['SampleA', 'SampleB', 'SampleC'])
    assert {'SampleA': 1, 'SampleB': 2, 'SampleC': 3} == sample_ids_list


def test_count_samples_associated_with_reference(sample_service, variation_service):
    assert 3 == sample_service.count_samples_associated_with_reference('genome')


def test_count_samples_associated_with_reference_empty(sample_service, variation_service):
    assert 0 == sample_service.count_samples_associated_with_reference('no_exist')


def test_samples_associated_with_reference_empty(sample_service, variation_service):
    samples_with_variants = sample_service.get_samples_associated_with_reference('no_exist')
    assert set() == {sample.name for sample in samples_with_variants}


def test_samples_with_variants_empty(sample_service, variation_service):
    samples_with_variants = sample_service.get_samples_with_variants('does_not_exist')
    assert set() == set(samples_with_variants)


def test_samples_which_exists(sample_service, variation_service):
    assert {'SampleA', 'SampleB', 'SampleC'} == set(sample_service.which_exists(['SampleA', 'SampleB',
                                                                                 'SampleC', 'not_exist']))


def test_samples_which_exists_only_one(sample_service, variation_service):
    assert {'SampleA'} == set(sample_service.which_exists(['SampleA']))


def test_samples_which_exists_none(sample_service, variation_service):
    assert [] == sample_service.which_exists(['not_exist'])


def test_samples_which_exists_none2(sample_service, variation_service):
    assert [] == sample_service.which_exists([])


def test_get_samples(sample_service, variation_service):
    assert {'SampleA', 'SampleB', 'SampleC'} == {sample.name for sample in sample_service.get_samples()}


def test_get_samples_empty(sample_service):
    assert [] == sample_service.get_samples()


def test_find_samples_by_ids(sample_service, variation_service):
    sample_name_ids = sample_service.find_sample_name_ids(['SampleA', 'SampleB', 'SampleC'])
    sample_set = SampleSet(sample_ids=[sample_name_ids['SampleA']])

    samples = sample_service.find_samples_by_ids(sample_set)
    assert {'SampleA'} == {s.name for s in samples}


def test_find_samples_by_ids_2_samples(sample_service, variation_service):
    sample_name_ids = sample_service.find_sample_name_ids(['SampleA', 'SampleB', 'SampleC'])
    sample_set = SampleSet(sample_ids=[sample_name_ids['SampleA'], sample_name_ids['SampleC']])

    samples = sample_service.find_samples_by_ids(sample_set)
    assert {'SampleA', 'SampleC'} == {s.name for s in samples}


def test_find_samples_by_ids_3_samples(sample_service, variation_service):
    sample_name_ids = sample_service.find_sample_name_ids(['SampleA', 'SampleB', 'SampleC'])
    sample_set = SampleSet(sample_ids=[sample_name_ids['SampleA'], sample_name_ids['SampleB'],
                                       sample_name_ids['SampleC']])

    samples = sample_service.find_samples_by_ids(sample_set)
    assert {'SampleA', 'SampleB', 'SampleC'} == {s.name for s in samples}


def test_find_samples_by_variation_ids(sample_service, variation_service):
    variant_samples = sample_service.find_samples_by_variation_ids(['reference:5061:G:A'])

    assert 'reference:5061:G:A' in variant_samples
    assert {'SampleB'} == {s.name for s in variant_samples['reference:5061:G:A']}
    assert {2} == {s.id for s in variant_samples['reference:5061:G:A']}
