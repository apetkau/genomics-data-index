def test_sample_sequence(sample_sequence_service):
    sample_sequences = sample_sequence_service.get_sample_sequences(reference_name='genome',
                                                                    samples=['SampleA', 'SampleB', 'SampleC'])

    assert {'reference'} == set(sample_sequences.keys()), 'Incorrect number of sequences'

    sample_names = {ss.sample.name for ss in sample_sequences['reference']}
    assert {'SampleA', 'SampleB', 'SampleC'} == sample_names


def test_missing_in_sequence(sample_sequence_service):
    assert sample_sequence_service.missing_in_sequence('SampleA', 'reference', positions=[1])
    assert sample_sequence_service.missing_in_sequence('SampleA', 'reference', positions=[1, 2])

    assert not sample_sequence_service.missing_in_sequence('SampleA', 'reference', positions=[500])
    assert sample_sequence_service.missing_in_sequence('SampleA', 'reference', positions=[1, 500])

    assert not sample_sequence_service.missing_in_sequence('SampleA', 'reference', positions=[])
    assert not sample_sequence_service.missing_in_sequence('SampleA', 'reference', positions=None)
