from storage.api.GenomicDataStore import GenomicDataStore


def test_summaries_loaded_data(loaded_database_genomic_data_store: GenomicDataStore):
    gds = loaded_database_genomic_data_store

    assert 9 == gds.count_samples()
    assert {'2014C-3598', '2014C-3599', '2014D-0067', '2014D-0068',
            'CFSAN002349', 'CFSAN023463',
            'SampleA', 'SampleB', 'SampleC'} == set(gds.sample_names())
    assert 112 == gds.count_mutations('genome')

    ms = gds.mutations_summary('genome')
    print(ms)
    assert 112 == len(ms)
    assert 2 == ms.loc['reference:839:1:G']
    assert 1 == ms.loc['reference:866:9:G']
    assert 1 == ms.loc['reference:1048:1:G']
    assert 2 == ms.loc['reference:3897:5:G']
