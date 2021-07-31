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

    summary_df = mlst_summarizier.summary(present_set)
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


def test_unique_summary(loaded_database_genomic_data_store: GenomicsDataIndex):
    db = loaded_database_genomic_data_store.connection.database
    all_sample_ids = {s.id for s in db.get_session().query(Sample).all()}
    sample_CFSAN002349 = db.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()
    sample_2014D_0068 = db.get_session().query(Sample).filter(Sample.name == '2014D-0068').one()
    sample_2014D_0067 = db.get_session().query(Sample).filter(Sample.name == '2014D-0067').one()
    sample_2014C_3598 = db.get_session().query(Sample).filter(Sample.name == '2014C-3598').one()
    sample_2014C_3599 = db.get_session().query(Sample).filter(Sample.name == '2014C-3599').one()
    assert 9 == len(all_sample_ids)

    mlst_summarizier = MLSTFeaturesSummarizer(connection=loaded_database_genomic_data_store.connection)

    # Test unique on all (should give me identical results to all since nothing is absent from the selection)
    present_set = SampleSet(all_sample_ids)
    complement_set = SampleSet.create_empty()
    summary_df = mlst_summarizier.unique_summary(present_set, other_set=complement_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    assert 24 == len(summary_df)
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Total', 'Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'abcZ', '1', 5, 9, 55] == summary_df.loc['mlst:lmonocytogenes:abcZ:1'].tolist()

    # Test unique on single sample (only a singl feature)
    present_set = SampleSet({sample_CFSAN002349.id})
    complement_set = SampleSet(all_sample_ids - {sample_CFSAN002349.id})
    summary_df = mlst_summarizier.unique_summary(present_set, other_set=complement_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    assert 1 == len(summary_df)
    print(summary_df)
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Total', 'Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'lhkA', '4', 1, 1, 100] == summary_df.loc['mlst:lmonocytogenes:lhkA:4'].tolist()

    # Test unique on two samples
    present_set = SampleSet({sample_CFSAN002349.id, sampleC.id})
    complement_set = SampleSet(all_sample_ids - {sample_CFSAN002349.id, sampleC.id})
    summary_df = mlst_summarizier.unique_summary(present_set, other_set=complement_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    assert 2 == len(summary_df)
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Total', 'Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'lhkA', '4', 1, 2, 50] == summary_df.loc['mlst:lmonocytogenes:lhkA:4'].tolist()
    assert ['lmonocytogenes', 'cat', '12', 1, 2, 50] == summary_df.loc['mlst:lmonocytogenes:cat:12'].tolist()

    # Test unique within a scheme
    present_set = SampleSet({sample_2014C_3598.id, sample_2014C_3599.id})
    complement_set = SampleSet(all_sample_ids - {sample_2014C_3598.id, sample_2014C_3599.id})
    summary_df = mlst_summarizier.unique_summary(present_set, other_set=complement_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    assert 7 == len(summary_df)
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Total', 'Percent'] == list(summary_df.columns)
    assert ['ecoli', 'adk', '100', 2, 2, 100] == summary_df.loc['mlst:ecoli:adk:100'].tolist()
    assert ['ecoli', 'fumC', '23', 2, 2, 100] == summary_df.loc['mlst:ecoli:fumC:23'].tolist()
    assert ['ecoli', 'gyrB', '68', 2, 2, 100] == summary_df.loc['mlst:ecoli:gyrB:68'].tolist()
    assert ['ecoli', 'icd', '45', 2, 2, 100] == summary_df.loc['mlst:ecoli:icd:45'].tolist()
    assert ['ecoli', 'mdh', '1', 2, 2, 100] == summary_df.loc['mlst:ecoli:mdh:1'].tolist()
    assert ['ecoli', 'purA', '35', 2, 2, 100] == summary_df.loc['mlst:ecoli:purA:35'].tolist()
    assert ['ecoli', 'recA', '7', 2, 2, 100] == summary_df.loc['mlst:ecoli:recA:7'].tolist()

    # Test unique across schemes
    present_set = SampleSet({sample_CFSAN002349.id, sample_2014D_0068.id})
    complement_set = SampleSet(all_sample_ids - {sample_CFSAN002349.id, sample_2014D_0068.id})
    summary_df = mlst_summarizier.unique_summary(present_set, other_set=complement_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    assert 2 == len(summary_df)
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Total', 'Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'lhkA', '4', 1, 2, 50] == summary_df.loc['mlst:lmonocytogenes:lhkA:4'].tolist()
    assert ['campylobacter', 'uncA', '6', 1, 2, 50] == summary_df.loc['mlst:campylobacter:uncA:6'].tolist()

    # Test unique only unknown
    mlst_summarizier = MLSTFeaturesSummarizer(connection=loaded_database_genomic_data_store.connection,
                                              include_present=False, include_unknown=True)

    present_set = SampleSet({sampleB.id, sample_2014D_0067.id})
    complement_set = SampleSet(all_sample_ids - {sampleB.id, sample_2014D_0067.id})
    summary_df = mlst_summarizier.unique_summary(present_set, other_set=complement_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    assert 2 == len(summary_df)
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Total', 'Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'ldh', '?', 1, 2, 50] == summary_df.loc['mlst:lmonocytogenes:ldh:?'].tolist()
    assert ['campylobacter', 'uncA', '?', 1, 2, 50] == summary_df.loc['mlst:campylobacter:uncA:?'].tolist()

    # Test unique only unknown, restricted to specific scheme
    mlst_summarizier = MLSTFeaturesSummarizer(connection=loaded_database_genomic_data_store.connection,
                                              include_present=False, include_unknown=True,
                                              scheme='lmonocytogenes')

    present_set = SampleSet({sampleB.id, sample_2014D_0067.id})
    complement_set = SampleSet(all_sample_ids - {sampleB.id, sample_2014D_0067.id})
    summary_df = mlst_summarizier.unique_summary(present_set, other_set=complement_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    assert 1 == len(summary_df)
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Total', 'Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'ldh', '?', 1, 2, 50] == summary_df.loc['mlst:lmonocytogenes:ldh:?'].tolist()


def test_summary_selections(loaded_database_genomic_data_store: GenomicsDataIndex):
    db = loaded_database_genomic_data_store.connection.database
    all_sample_ids = {s.id for s in db.get_session().query(Sample).all()}
    assert 9 == len(all_sample_ids)
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sample_CFSAN002349 = db.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()
    sample_2014D_0067 = db.get_session().query(Sample).filter(Sample.name == '2014D-0067').one()

    mlst_summarizer = MLSTFeaturesSummarizer(connection=loaded_database_genomic_data_store.connection)

    # Test only single sample features
    present_set = SampleSet([sampleA.id])

    summary_df = mlst_summarizer.summary(present_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    assert 7 == len(summary_df)
    assert {'lmonocytogenes'} == set(summary_df['Scheme'].tolist())
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Total', 'Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'abcZ', '1', 1, 1, 100] == summary_df.loc['mlst:lmonocytogenes:abcZ:1'].tolist()
    assert ['lmonocytogenes', 'bglA', '51', 1, 1, 100] == summary_df.loc['mlst:lmonocytogenes:bglA:51'].tolist()

    # Test two samples
    present_set = SampleSet([sampleA.id, sampleB.id])
    summary_df = mlst_summarizer.summary(present_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    assert 8 == len(summary_df)
    assert {'lmonocytogenes'} == set(summary_df['Scheme'].tolist())
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Total', 'Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'abcZ', '1', 2, 2, 100] == summary_df.loc['mlst:lmonocytogenes:abcZ:1'].tolist()
    assert ['lmonocytogenes', 'bglA', '51', 1, 2, 50] == summary_df.loc['mlst:lmonocytogenes:bglA:51'].tolist()
    assert ['lmonocytogenes', 'bglA', '52', 1, 2, 50] == summary_df.loc['mlst:lmonocytogenes:bglA:52'].tolist()
    assert ['lmonocytogenes', 'ldh', '5', 1, 2, 50] == summary_df.loc['mlst:lmonocytogenes:ldh:5'].tolist()

    # Test three samples
    present_set = SampleSet([sampleA.id, sampleB.id, sample_CFSAN002349.id])
    summary_df = mlst_summarizer.summary(present_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    assert 9 == len(summary_df)
    assert {'lmonocytogenes'} == set(summary_df['Scheme'].tolist())
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
    summary_df = mlst_summarizer.summary(present_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    assert 15 == len(summary_df)
    assert {'lmonocytogenes', 'campylobacter'} == set(summary_df['Scheme'].tolist())
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
    assert 6 == len(summary_df[summary_df['Scheme'] == 'campylobacter'])  # Missing one feature since it's unknown

    # Test multiple schemes sample set but summarize for only a particular scheme
    present_set = SampleSet([sampleA.id, sampleB.id, sample_CFSAN002349.id, sample_2014D_0067.id])
    mlst_summarizer = MLSTFeaturesSummarizer(connection=loaded_database_genomic_data_store.connection,
                                             scheme='lmonocytogenes')
    summary_df = mlst_summarizer.summary(present_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    assert 9 == len(summary_df)
    assert {'lmonocytogenes'} == set(summary_df['Scheme'].tolist())  # Only results for one scheme
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Total', 'Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'abcZ', '1', 3, 4, 75] == summary_df.loc['mlst:lmonocytogenes:abcZ:1'].tolist()
    assert ['lmonocytogenes', 'bglA', '51', 2, 4, 50] == summary_df.loc['mlst:lmonocytogenes:bglA:51'].tolist()
    assert ['lmonocytogenes', 'bglA', '52', 1, 4, 25] == summary_df.loc['mlst:lmonocytogenes:bglA:52'].tolist()
    assert ['lmonocytogenes', 'ldh', '5', 2, 4, 50] == summary_df.loc['mlst:lmonocytogenes:ldh:5'].tolist()
    assert ['lmonocytogenes', 'lhkA', '5', 2, 4, 50] == summary_df.loc['mlst:lmonocytogenes:lhkA:5'].tolist()
    assert ['lmonocytogenes', 'lhkA', '4', 1, 4, 25] == summary_df.loc['mlst:lmonocytogenes:lhkA:4'].tolist()

    # Test multiple schemes sample set but summarize for only a particular scheme/locus
    present_set = SampleSet([sampleA.id, sampleB.id, sample_CFSAN002349.id, sample_2014D_0067.id])
    mlst_summarizer = MLSTFeaturesSummarizer(connection=loaded_database_genomic_data_store.connection,
                                             scheme='lmonocytogenes', locus='bglA')
    summary_df = mlst_summarizer.summary(present_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    assert 2 == len(summary_df)
    assert {'lmonocytogenes'} == set(summary_df['Scheme'].tolist())  # Only results for one scheme
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Total', 'Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'bglA', '51', 2, 4, 50] == summary_df.loc['mlst:lmonocytogenes:bglA:51'].tolist()
    assert ['lmonocytogenes', 'bglA', '52', 1, 4, 25] == summary_df.loc['mlst:lmonocytogenes:bglA:52'].tolist()

    # Test multiple schemes, include unknown
    present_set = SampleSet([sampleA.id, sampleB.id, sample_CFSAN002349.id, sample_2014D_0067.id])
    mlst_summarizer = MLSTFeaturesSummarizer(connection=loaded_database_genomic_data_store.connection,
                                             include_unknown=True)
    summary_df = mlst_summarizer.summary(present_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    assert 17 == len(summary_df)
    assert {'lmonocytogenes', 'campylobacter'} == set(summary_df['Scheme'].tolist())
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Total', 'Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'abcZ', '1', 3, 4, 75] == summary_df.loc['mlst:lmonocytogenes:abcZ:1'].tolist()
    assert ['lmonocytogenes', 'bglA', '51', 2, 4, 50] == summary_df.loc['mlst:lmonocytogenes:bglA:51'].tolist()
    assert ['lmonocytogenes', 'bglA', '52', 1, 4, 25] == summary_df.loc['mlst:lmonocytogenes:bglA:52'].tolist()
    assert ['lmonocytogenes', 'ldh', '5', 2, 4, 50] == summary_df.loc['mlst:lmonocytogenes:ldh:5'].tolist()
    assert ['lmonocytogenes', 'ldh', '?', 1, 4, 25] == summary_df.loc['mlst:lmonocytogenes:ldh:?'].tolist()
    assert ['lmonocytogenes', 'lhkA', '5', 2, 4, 50] == summary_df.loc['mlst:lmonocytogenes:lhkA:5'].tolist()
    assert ['lmonocytogenes', 'lhkA', '4', 1, 4, 25] == summary_df.loc['mlst:lmonocytogenes:lhkA:4'].tolist()
    assert ['campylobacter', 'aspA', '2', 1, 4, 25] == summary_df.loc['mlst:campylobacter:aspA:2'].tolist()
    assert ['campylobacter', 'glyA', '3', 1, 4, 25] == summary_df.loc['mlst:campylobacter:glyA:3'].tolist()
    assert ['campylobacter', 'uncA', '?', 1, 4, 25] == summary_df.loc['mlst:campylobacter:uncA:?'].tolist()

    # Test multiple schemes, only unknown
    present_set = SampleSet([sampleA.id, sampleB.id, sample_CFSAN002349.id, sample_2014D_0067.id])
    mlst_summarizer = MLSTFeaturesSummarizer(connection=loaded_database_genomic_data_store.connection,
                                             include_present=False, include_unknown=True)
    summary_df = mlst_summarizer.summary(present_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    assert 2 == len(summary_df)
    assert {'lmonocytogenes', 'campylobacter'} == set(summary_df['Scheme'].tolist())
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Total', 'Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'ldh', '?', 1, 4, 25] == summary_df.loc['mlst:lmonocytogenes:ldh:?'].tolist()
    assert ['campylobacter', 'uncA', '?', 1, 4, 25] == summary_df.loc['mlst:campylobacter:uncA:?'].tolist()
