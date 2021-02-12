def test_samples_with_variants(sample_service, variation_service):
    samples_with_variants = sample_service.get_samples_with_variants('genome')
    assert {'SampleA', 'SampleB', 'SampleC'} == {sample.name for sample in samples_with_variants}


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
