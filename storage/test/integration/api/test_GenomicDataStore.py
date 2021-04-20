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
    assert 112 == len(ms)
    assert 2 == ms.loc['reference:839:1:G']
    assert 1 == ms.loc['reference:866:9:G']
    assert 1 == ms.loc['reference:1048:1:G']
    assert 2 == ms.loc['reference:3897:5:G']

    ms = gds.mutations_summary('genome', id_type='spdi_ref')
    assert 112 == len(ms)
    assert 2 == ms.loc['reference:839:C:G']
    assert 1 == ms.loc['reference:866:GCCAGATCC:G']
    assert 1 == ms.loc['reference:1048:C:G']
    assert 2 == ms.loc['reference:3897:GCGCA:G']
