from storage.api.GenomicDataStore import GenomicDataStore


def test_summaries_loaded_data(loaded_database_genomic_data_store: GenomicDataStore):
    gds = loaded_database_genomic_data_store

    # Samples
    assert 9 == gds.count_samples()
    assert {'2014C-3598', '2014C-3599', '2014D-0067', '2014D-0068',
            'CFSAN002349', 'CFSAN023463',
            'SampleA', 'SampleB', 'SampleC'} == set(gds.sample_names())

    # References
    assert 1 == gds.count_references()
    assert ['genome'] == gds.reference_names()

    # Mutations
    assert 112 == gds.count_mutations('genome')
    ms = gds.mutations_summary('genome', id_type='spdi')
    print(ms)
    assert 112 == len(ms)
    assert 'Mutation' == ms.index.name
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Count'] == list(ms.columns)
    assert ['reference', 839, 1, 'G', 2] == ms.loc['reference:839:1:G'].values.tolist()
    assert ['reference', 866, 9, 'G', 1] == ms.loc['reference:866:9:G'].values.tolist()
    assert ['reference', 1048, 1, 'G', 1] == ms.loc['reference:1048:1:G'].values.tolist()
    assert ['reference', 3897, 5, 'G', 2] == ms.loc['reference:3897:5:G'].values.tolist()

    ms = gds.mutations_summary('genome', id_type='spdi_ref')
    print(ms)
    assert 112 == len(ms)
    assert 'Mutation' == ms.index.name
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Count'] == list(ms.columns)
    assert ['reference', 839, 'C', 'G', 2] == ms.loc['reference:839:C:G'].values.tolist()
    assert ['reference', 866, 'GCCAGATCC', 'G', 1] == ms.loc['reference:866:GCCAGATCC:G'].values.tolist()
    assert ['reference', 1048, 'C', 'G', 1] == ms.loc['reference:1048:C:G'].values.tolist()
    assert ['reference', 3897, 'GCGCA', 'G', 2] == ms.loc['reference:3897:GCGCA:G'].values.tolist()
