def test_sample_sequence(sample_sequence_service):
    sample_sequences = sample_sequence_service.get_sample_sequences(reference_name='genome',
                                                                    samples=['SampleA', 'SampleB', 'SampleC'])

    assert {'reference'} == set(sample_sequences.keys()), 'Incorrect number of sequences'

    sample_names = {ss.sample.name for ss in sample_sequences['reference']}
    assert {'SampleA', 'SampleB', 'SampleC'} == sample_names
