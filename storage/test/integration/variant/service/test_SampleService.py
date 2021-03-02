def test_samples_with_variants(sample_service, variation_service):
    samples_with_variants = sample_service.get_samples_with_variants('genome')
    assert {'SampleA', 'SampleB', 'SampleC'} == {sample.name for sample in samples_with_variants}


def test_samples_with_variants_on_sequence(sample_service, variation_service):
    samples_with_variants = sample_service.get_samples_with_variants_on_sequence('reference')
    assert {'SampleA', 'SampleB', 'SampleC'} == {sample.name for sample in samples_with_variants}


def test_samples_associated_with_sequence(sample_service, variation_service):
    samples_with_variants = sample_service.get_samples_associated_with_sequence('reference')
    assert {'SampleA', 'SampleB', 'SampleC'} == {sample.name for sample in samples_with_variants}


def test_samples_associated_with_sequence_empty(sample_service, variation_service):
    samples_with_variants = sample_service.get_samples_associated_with_sequence('no_exist')
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


def test_find_samples_by_variation_ids(sample_service, variation_service):
    variant_samples = sample_service.find_samples_by_variation_ids(['reference:5061:G:A'])

    assert 'reference:5061:G:A' in variant_samples
    assert {'SampleB'} == {s.name for s in variant_samples['reference:5061:G:A']}
    assert {2} == {s.id for s in variant_samples['reference:5061:G:A']}
