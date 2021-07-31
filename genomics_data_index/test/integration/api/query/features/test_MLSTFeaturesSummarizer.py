from genomics_data_index.api.query.GenomicsDataIndex import GenomicsDataIndex
from genomics_data_index.api.query.features.MLSTFeaturesSummarizer import MLSTFeaturesSummarizer
from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.db import Sample


def test_summary_all(loaded_database_genomic_data_store: GenomicsDataIndex):
    db = loaded_database_genomic_data_store.connection.database
    all_sample_ids = {s.id for s in db.get_session().query(Sample).all()}
    assert 9 == len(all_sample_ids)

    mlst_summarizier = MLSTFeaturesSummarizer(connection=loaded_database_genomic_data_store.connection)

    present_set = SampleSet(all_sample_ids)
    unknown_set = SampleSet.create_empty()
    absent_set = SampleSet.create_empty()

    summary_df = mlst_summarizier.summary(present_samples=present_set,
                                          unknown_samples=unknown_set,
                                          absent_samples=absent_set,
                                          selection='all')
    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    assert 24 == len(summary_df)
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Total', 'Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'abcZ', '1', 5, 9, 55] == summary_df.loc['mlst:lmonocytogenes:abcZ:1'].tolist()
    assert ['lmonocytogenes', 'bglA', '51', 3, 9, 33] == summary_df.loc['mlst:lmonocytogenes:bglA:51'].tolist()
    assert ['lmonocytogenes', 'bglA', '52', 2, 9, 22] == summary_df.loc['mlst:lmonocytogenes:bglA:52'].tolist()
    assert ['ecoli', 'adk', '100', 2, 9, 22] == summary_df.loc['mlst:ecoli:adk:100'].tolist()
    assert ['ecoli', 'recA', '7', 2, 9, 22] == summary_df.loc['mlst:ecoli:recA:7'].tolist()
    assert ['campylobacter', 'uncA', '6', 1, 9, 11] == summary_df.loc['mlst:campylobacter:uncA:6'].tolist()

    # Test unique (should give me identical results since nothing is absent from the selection)
    # summary_df = mlst_summarizier.summary(present_samples=present_set,
    #                                       unknown_samples=unknown_set,
    #                                       absent_samples=absent_set,
    #                                       selection='unique')
    # summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    # assert 24 == len(summary_df)
    # assert 'MLST Feature' == summary_df.index.name
    # assert ['Scheme', 'Locus', 'Allele', 'Count', 'Total', 'Percent'] == list(summary_df.columns)
    # assert ['lmonocytogenes', 'abcZ', '1', 5, 9, 55] == summary_df.loc['mlst:lmonocytogenes:abcZ:1'].tolist()


def test_summary_selections(loaded_database_genomic_data_store: GenomicsDataIndex):
    db = loaded_database_genomic_data_store.connection.database
    all_sample_ids = {s.id for s in db.get_session().query(Sample).all()}
    assert 9 == len(all_sample_ids)
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sample_CFSAN002349 = db.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()
    sample_2014D_0067 = db.get_session().query(Sample).filter(Sample.name == '2014D-0067').one()

    mlst_summarizier = MLSTFeaturesSummarizer(connection=loaded_database_genomic_data_store.connection)

    # Test only single sample features
    present_set = SampleSet([sampleA.id])
    unknown_set = SampleSet.create_empty()
    absent_set = SampleSet(all_sample_ids - {sampleA.id})

    summary_df = mlst_summarizier.summary(present_samples=present_set,
                                          unknown_samples=unknown_set,
                                          absent_samples=absent_set,
                                          selection='all')

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    assert 7 == len(summary_df)
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Total', 'Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'abcZ', '1', 1, 1, 100] == summary_df.loc['mlst:lmonocytogenes:abcZ:1'].tolist()
    assert ['lmonocytogenes', 'bglA', '51', 1, 1, 100] == summary_df.loc['mlst:lmonocytogenes:bglA:51'].tolist()

    # Test two samples
    present_set = SampleSet([sampleA.id, sampleB.id])
    unknown_set = SampleSet.create_empty()
    absent_set = SampleSet(all_sample_ids - {sampleA.id, sampleB.id})

    summary_df = mlst_summarizier.summary(present_samples=present_set,
                                          unknown_samples=unknown_set,
                                          absent_samples=absent_set,
                                          selection='all')

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    assert 8 == len(summary_df)
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Total', 'Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'abcZ', '1', 2, 2, 100] == summary_df.loc['mlst:lmonocytogenes:abcZ:1'].tolist()
    assert ['lmonocytogenes', 'bglA', '51', 1, 2, 50] == summary_df.loc['mlst:lmonocytogenes:bglA:51'].tolist()
    assert ['lmonocytogenes', 'bglA', '52', 1, 2, 50] == summary_df.loc['mlst:lmonocytogenes:bglA:52'].tolist()
    assert ['lmonocytogenes', 'ldh', '5', 1, 2, 50] == summary_df.loc['mlst:lmonocytogenes:ldh:5'].tolist()

    # Test three samples
    present_set = SampleSet([sampleA.id, sampleB.id, sample_CFSAN002349.id])
    unknown_set = SampleSet.create_empty()
    absent_set = SampleSet(all_sample_ids - {sampleA.id, sampleB.id, sample_CFSAN002349.id})

    summary_df = mlst_summarizier.summary(present_samples=present_set,
                                          unknown_samples=unknown_set,
                                          absent_samples=absent_set,
                                          selection='all')

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    assert 9 == len(summary_df)
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Total', 'Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'abcZ', '1', 3, 3, 100] == summary_df.loc['mlst:lmonocytogenes:abcZ:1'].tolist()
    assert ['lmonocytogenes', 'bglA', '51', 2, 3, 66] == summary_df.loc['mlst:lmonocytogenes:bglA:51'].tolist()
    assert ['lmonocytogenes', 'bglA', '52', 1, 3, 33] == summary_df.loc['mlst:lmonocytogenes:bglA:52'].tolist()
    assert ['lmonocytogenes', 'ldh', '5', 2, 3, 66] == summary_df.loc['mlst:lmonocytogenes:ldh:5'].tolist()
    assert ['lmonocytogenes', 'lhkA', '5', 2, 3, 66] == summary_df.loc['mlst:lmonocytogenes:lhkA:5'].tolist()
    assert ['lmonocytogenes', 'lhkA', '4', 1, 3, 33] == summary_df.loc['mlst:lmonocytogenes:lhkA:4'].tolist()

    # Test multiple schemes
    present_set = SampleSet([sampleA.id, sampleB.id, sample_CFSAN002349.id, sample_2014D_0067.id])
    unknown_set = SampleSet.create_empty()
    absent_set = SampleSet(all_sample_ids - {sampleA.id, sampleB.id, sample_CFSAN002349.id, sample_2014D_0067.id})

    summary_df = mlst_summarizier.summary(present_samples=present_set,
                                          unknown_samples=unknown_set,
                                          absent_samples=absent_set,
                                          selection='all')

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    assert 15 == len(summary_df)
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Total', 'Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'abcZ', '1', 3, 4, 75] == summary_df.loc['mlst:lmonocytogenes:abcZ:1'].tolist()
    assert ['lmonocytogenes', 'bglA', '51', 2, 4, 50] == summary_df.loc['mlst:lmonocytogenes:bglA:51'].tolist()
    assert ['lmonocytogenes', 'bglA', '52', 1, 4, 25] == summary_df.loc['mlst:lmonocytogenes:bglA:52'].tolist()
    assert ['lmonocytogenes', 'ldh', '5', 2, 4, 50] == summary_df.loc['mlst:lmonocytogenes:ldh:5'].tolist()
    assert ['lmonocytogenes', 'lhkA', '5', 2, 4, 50] == summary_df.loc['mlst:lmonocytogenes:lhkA:5'].tolist()
    assert ['lmonocytogenes', 'lhkA', '4', 1, 4, 25] == summary_df.loc['mlst:lmonocytogenes:lhkA:4'].tolist()
    assert ['campylobacter', 'aspA', '2', 1, 4, 25] == summary_df.loc['mlst:campylobacter:aspA:2'].tolist()
    assert ['campylobacter', 'glyA', '3', 1, 4, 25] == summary_df.loc['mlst:campylobacter:glyA:3'].tolist()
    assert 6 == len(summary_df[summary_df['Scheme'] == 'campylobacter']) # Missing one feature since it's unknown
