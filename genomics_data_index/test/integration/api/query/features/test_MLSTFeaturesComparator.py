from genomics_data_index.api.query.GenomicsDataIndex import GenomicsDataIndex
from genomics_data_index.api.query.features.MLSTFeaturesComparator import MLSTFeaturesComparator
from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.db import Sample


def test_summary_all(loaded_database_genomic_data_store: GenomicsDataIndex):
    db = loaded_database_genomic_data_store.connection.database
    all_sample_ids = {s.id for s in db.get_session().query(Sample).all()}
    assert 9 == len(all_sample_ids)

    mlst_summarizier = MLSTFeaturesComparator(connection=loaded_database_genomic_data_store.connection,
                                              include_unknown_samples=True)

    present_set = SampleSet(all_sample_ids)

    summary_df = mlst_summarizier.summary(present_set)
    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    summary_df['Unknown Percent'] = summary_df['Unknown Percent'].astype(int)
    summary_df['Present and Unknown Percent'] = summary_df['Present and Unknown Percent'].astype(int)
    assert 24 == len(summary_df)
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Unknown Count', 'Present and Unknown Count',
            'Total', 'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'abcZ', '1', 5, 0, 5, 9, 55, 0, 55] == summary_df.loc[
        'mlst:lmonocytogenes:abcZ:1'].tolist()
    assert ['lmonocytogenes', 'bglA', '51', 3, 0, 3, 9, 33, 0, 33] == summary_df.loc[
        'mlst:lmonocytogenes:bglA:51'].tolist()
    assert ['lmonocytogenes', 'bglA', '52', 2, 0, 2, 9, 22, 0, 22] == summary_df.loc[
        'mlst:lmonocytogenes:bglA:52'].tolist()
    assert ['ecoli', 'adk', '100', 2, 0, 2, 9, 22, 0, 22] == summary_df.loc['mlst:ecoli:adk:100'].tolist()
    assert ['ecoli', 'recA', '7', 2, 0, 2, 9, 22, 0, 22] == summary_df.loc['mlst:ecoli:recA:7'].tolist()
    assert ['campylobacter', 'uncA', '6', 1, 1, 2, 9, 11, 11, 22] == summary_df.loc[
        'mlst:campylobacter:uncA:6'].tolist()


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

    mlst_summarizier = MLSTFeaturesComparator(connection=loaded_database_genomic_data_store.connection)

    # Test unique on all (should give me identical results to all since nothing is absent from the selection)
    present_set = SampleSet(all_sample_ids)
    complement_set = SampleSet.create_empty()
    summary_df = mlst_summarizier.unique_summary(present_set, other_set=complement_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    summary_df['Unknown Percent'] = summary_df['Unknown Percent'].astype(int)
    summary_df['Present and Unknown Percent'] = summary_df['Present and Unknown Percent'].astype(int)
    assert 24 == len(summary_df)
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Unknown Count', 'Present and Unknown Count',
            'Total', 'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'abcZ', '1', 5, 0, 5, 9, 55, 0, 55] == summary_df.loc[
        'mlst:lmonocytogenes:abcZ:1'].tolist()

    # Test unique on single sample (only a single feature)
    present_set = SampleSet({sample_CFSAN002349.id})
    complement_set = SampleSet(all_sample_ids - {sample_CFSAN002349.id})
    summary_df = mlst_summarizier.unique_summary(present_set, other_set=complement_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    summary_df['Unknown Percent'] = summary_df['Unknown Percent'].astype(int)
    summary_df['Present and Unknown Percent'] = summary_df['Present and Unknown Percent'].astype(int)
    assert 1 == len(summary_df)
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Unknown Count', 'Present and Unknown Count',
            'Total', 'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'lhkA', '4', 1, 0, 1, 1, 100, 0, 100] == summary_df.loc[
        'mlst:lmonocytogenes:lhkA:4'].tolist()

    # Test unique on two samples
    present_set = SampleSet({sample_CFSAN002349.id, sampleC.id})
    complement_set = SampleSet(all_sample_ids - {sample_CFSAN002349.id, sampleC.id})
    summary_df = mlst_summarizier.unique_summary(present_set, other_set=complement_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    summary_df['Unknown Percent'] = summary_df['Unknown Percent'].astype(int)
    summary_df['Present and Unknown Percent'] = summary_df['Present and Unknown Percent'].astype(int)
    assert 2 == len(summary_df)
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Unknown Count', 'Present and Unknown Count',
            'Total', 'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'lhkA', '4', 1, 0, 1, 2, 50, 0, 50] == summary_df.loc[
        'mlst:lmonocytogenes:lhkA:4'].tolist()
    assert ['lmonocytogenes', 'cat', '12', 1, 0, 1, 2, 50, 0, 50] == summary_df.loc[
        'mlst:lmonocytogenes:cat:12'].tolist()

    # Test unique within a scheme
    present_set = SampleSet({sample_2014C_3598.id, sample_2014C_3599.id})
    complement_set = SampleSet(all_sample_ids - {sample_2014C_3598.id, sample_2014C_3599.id})
    summary_df = mlst_summarizier.unique_summary(present_set, other_set=complement_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    summary_df['Unknown Percent'] = summary_df['Unknown Percent'].astype(int)
    summary_df['Present and Unknown Percent'] = summary_df['Present and Unknown Percent'].astype(int)
    assert 7 == len(summary_df)
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Unknown Count', 'Present and Unknown Count',
            'Total', 'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(summary_df.columns)
    assert ['ecoli', 'adk', '100', 2, 0, 2, 2, 100, 0, 100] == summary_df.loc['mlst:ecoli:adk:100'].tolist()
    assert ['ecoli', 'fumC', '23', 2, 0, 2, 2, 100, 0, 100] == summary_df.loc['mlst:ecoli:fumC:23'].tolist()
    assert ['ecoli', 'gyrB', '68', 2, 0, 2, 2, 100, 0, 100] == summary_df.loc['mlst:ecoli:gyrB:68'].tolist()
    assert ['ecoli', 'icd', '45', 2, 0, 2, 2, 100, 0, 100] == summary_df.loc['mlst:ecoli:icd:45'].tolist()
    assert ['ecoli', 'mdh', '1', 2, 0, 2, 2, 100, 0, 100] == summary_df.loc['mlst:ecoli:mdh:1'].tolist()
    assert ['ecoli', 'purA', '35', 2, 0, 2, 2, 100, 0, 100] == summary_df.loc['mlst:ecoli:purA:35'].tolist()
    assert ['ecoli', 'recA', '7', 2, 0, 2, 2, 100, 0, 100] == summary_df.loc['mlst:ecoli:recA:7'].tolist()

    # Test unique across schemes
    present_set = SampleSet({sample_CFSAN002349.id, sample_2014D_0068.id})
    complement_set = SampleSet(all_sample_ids - {sample_CFSAN002349.id, sample_2014D_0068.id})
    summary_df = mlst_summarizier.unique_summary(present_set, other_set=complement_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    summary_df['Unknown Percent'] = summary_df['Unknown Percent'].astype(int)
    summary_df['Present and Unknown Percent'] = summary_df['Present and Unknown Percent'].astype(int)
    assert 2 == len(summary_df)
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Unknown Count', 'Present and Unknown Count',
            'Total', 'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'lhkA', '4', 1, 0, 1, 2, 50, 0, 50] == summary_df.loc[
        'mlst:lmonocytogenes:lhkA:4'].tolist()
    assert ['campylobacter', 'uncA', '6', 1, 0, 1, 2, 50, 0, 50] == summary_df.loc['mlst:campylobacter:uncA:6'].tolist()

    # Test unique across schemes, include unknown features when no samples present
    mlst_summarizier = MLSTFeaturesComparator(connection=loaded_database_genomic_data_store.connection,
                                              include_unknown_no_present_samples=True)

    present_set = SampleSet({sample_CFSAN002349.id, sample_2014D_0068.id})
    complement_set = SampleSet(all_sample_ids - {sample_CFSAN002349.id, sample_2014D_0068.id})
    summary_df = mlst_summarizier.unique_summary(present_set, other_set=complement_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    summary_df['Unknown Percent'] = summary_df['Unknown Percent'].astype(int)
    summary_df['Present and Unknown Percent'] = summary_df['Present and Unknown Percent'].astype(int)
    assert 1 == len(summary_df)
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Unknown Count', 'Present and Unknown Count',
            'Total', 'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'lhkA', '4', 1, 0, 1, 2, 50, 0, 50] == summary_df.loc[
        'mlst:lmonocytogenes:lhkA:4'].tolist()

    # Test unique only unknown
    mlst_summarizier = MLSTFeaturesComparator(connection=loaded_database_genomic_data_store.connection,
                                              include_present=False, include_unknown=True)

    present_set = SampleSet({sampleB.id, sample_2014D_0067.id})
    complement_set = SampleSet(all_sample_ids - {sampleB.id, sample_2014D_0067.id})
    summary_df = mlst_summarizier.unique_summary(present_set, other_set=complement_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    summary_df['Unknown Count'] = summary_df['Unknown Count'].fillna(-1).astype(int)
    summary_df['Present and Unknown Count'] = summary_df['Present and Unknown Count'].fillna(-1).astype(int)
    summary_df['Unknown Percent'] = summary_df['Unknown Percent'].fillna(-1).astype(int)
    summary_df['Present and Unknown Percent'] = summary_df['Present and Unknown Percent'].fillna(-1).astype(int)
    assert 2 == len(summary_df)
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Unknown Count', 'Present and Unknown Count',
            'Total', 'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'ldh', '?', 1, -1, -1, 2, 50, -1, -1] == summary_df.loc[
        'mlst:lmonocytogenes:ldh:?'].tolist()
    assert ['campylobacter', 'uncA', '?', 1, -1, -1, 2, 50, -1, -1] == summary_df.loc[
        'mlst:campylobacter:uncA:?'].tolist()

    # Test unique only unknown, restricted to specific scheme
    mlst_summarizier = MLSTFeaturesComparator(connection=loaded_database_genomic_data_store.connection,
                                              include_present=False, include_unknown=True,
                                              scheme='lmonocytogenes')

    present_set = SampleSet({sampleB.id, sample_2014D_0067.id})
    complement_set = SampleSet(all_sample_ids - {sampleB.id, sample_2014D_0067.id})
    summary_df = mlst_summarizier.unique_summary(present_set, other_set=complement_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    summary_df['Unknown Count'] = summary_df['Unknown Count'].fillna(-1).astype(int)
    summary_df['Present and Unknown Count'] = summary_df['Present and Unknown Count'].fillna(-1).astype(int)
    summary_df['Unknown Percent'] = summary_df['Unknown Percent'].fillna(-1).astype(int)
    summary_df['Present and Unknown Percent'] = summary_df['Present and Unknown Percent'].fillna(-1).astype(int)
    assert 1 == len(summary_df)
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Unknown Count', 'Present and Unknown Count',
            'Total', 'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'ldh', '?', 1, -1, -1, 2, 50, -1, -1] == summary_df.loc[
        'mlst:lmonocytogenes:ldh:?'].tolist()


def test_summary_selections(loaded_database_genomic_data_store: GenomicsDataIndex):
    db = loaded_database_genomic_data_store.connection.database
    all_sample_ids = {s.id for s in db.get_session().query(Sample).all()}
    assert 9 == len(all_sample_ids)
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sample_CFSAN002349 = db.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()
    sample_2014D_0067 = db.get_session().query(Sample).filter(Sample.name == '2014D-0067').one()

    mlst_summarizer = MLSTFeaturesComparator(connection=loaded_database_genomic_data_store.connection,
                                             include_unknown_samples=True)

    # Test only single sample features
    present_set = SampleSet([sampleA.id])

    summary_df = mlst_summarizer.summary(present_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    summary_df['Unknown Percent'] = summary_df['Unknown Percent'].astype(int)
    summary_df['Present and Unknown Percent'] = summary_df['Present and Unknown Percent'].astype(int)
    assert 7 == len(summary_df)
    assert {'lmonocytogenes'} == set(summary_df['Scheme'].tolist())
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Unknown Count', 'Present and Unknown Count',
            'Total', 'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'abcZ', '1', 1, 0, 1, 1, 100, 0, 100] == summary_df.loc[
        'mlst:lmonocytogenes:abcZ:1'].tolist()
    assert ['lmonocytogenes', 'bglA', '51', 1, 0, 1, 1, 100, 0, 100] == summary_df.loc[
        'mlst:lmonocytogenes:bglA:51'].tolist()

    # Test two samples
    present_set = SampleSet([sampleA.id, sampleB.id])
    summary_df = mlst_summarizer.summary(present_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    summary_df['Unknown Percent'] = summary_df['Unknown Percent'].astype(int)
    summary_df['Present and Unknown Percent'] = summary_df['Present and Unknown Percent'].astype(int)
    assert 8 == len(summary_df)
    assert {'lmonocytogenes'} == set(summary_df['Scheme'].tolist())
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Unknown Count', 'Present and Unknown Count',
            'Total', 'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'abcZ', '1', 2, 0, 2, 2, 100, 0, 100] == summary_df.loc[
        'mlst:lmonocytogenes:abcZ:1'].tolist()
    assert ['lmonocytogenes', 'bglA', '51', 1, 0, 1, 2, 50, 0, 50] == summary_df.loc[
        'mlst:lmonocytogenes:bglA:51'].tolist()
    assert ['lmonocytogenes', 'bglA', '52', 1, 0, 1, 2, 50, 0, 50] == summary_df.loc[
        'mlst:lmonocytogenes:bglA:52'].tolist()
    assert ['lmonocytogenes', 'ldh', '5', 1, 1, 2, 2, 50, 50, 100] == summary_df.loc[
        'mlst:lmonocytogenes:ldh:5'].tolist()

    # Test three samples
    present_set = SampleSet([sampleA.id, sampleB.id, sample_CFSAN002349.id])
    summary_df = mlst_summarizer.summary(present_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    summary_df['Unknown Percent'] = summary_df['Unknown Percent'].astype(int)
    summary_df['Present and Unknown Percent'] = summary_df['Present and Unknown Percent'].astype(int)
    assert 9 == len(summary_df)
    assert {'lmonocytogenes'} == set(summary_df['Scheme'].tolist())
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Unknown Count', 'Present and Unknown Count',
            'Total', 'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'abcZ', '1', 3, 0, 3, 3, 100, 0, 100] == summary_df.loc[
        'mlst:lmonocytogenes:abcZ:1'].tolist()
    assert ['lmonocytogenes', 'bglA', '51', 2, 0, 2, 3, 66, 0, 66] == summary_df.loc[
        'mlst:lmonocytogenes:bglA:51'].tolist()
    assert ['lmonocytogenes', 'bglA', '52', 1, 0, 1, 3, 33, 0, 33] == summary_df.loc[
        'mlst:lmonocytogenes:bglA:52'].tolist()
    assert ['lmonocytogenes', 'ldh', '5', 2, 1, 3, 3, 66, 33, 100] == summary_df.loc[
        'mlst:lmonocytogenes:ldh:5'].tolist()
    assert ['lmonocytogenes', 'lhkA', '5', 2, 0, 2, 3, 66, 0, 66] == summary_df.loc[
        'mlst:lmonocytogenes:lhkA:5'].tolist()
    assert ['lmonocytogenes', 'lhkA', '4', 1, 0, 1, 3, 33, 0, 33] == summary_df.loc[
        'mlst:lmonocytogenes:lhkA:4'].tolist()

    # Test multiple schemes
    present_set = SampleSet([sampleA.id, sampleB.id, sample_CFSAN002349.id, sample_2014D_0067.id])
    summary_df = mlst_summarizer.summary(present_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    summary_df['Unknown Percent'] = summary_df['Unknown Percent'].astype(int)
    summary_df['Present and Unknown Percent'] = summary_df['Present and Unknown Percent'].astype(int)
    assert 15 == len(summary_df)
    assert {'lmonocytogenes', 'campylobacter'} == set(summary_df['Scheme'].tolist())
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Unknown Count', 'Present and Unknown Count',
            'Total', 'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'abcZ', '1', 3, 0, 3, 4, 75, 0, 75] == summary_df.loc[
        'mlst:lmonocytogenes:abcZ:1'].tolist()
    assert ['lmonocytogenes', 'bglA', '51', 2, 0, 2, 4, 50, 0, 50] == summary_df.loc[
        'mlst:lmonocytogenes:bglA:51'].tolist()
    assert ['lmonocytogenes', 'bglA', '52', 1, 0, 1, 4, 25, 0, 25] == summary_df.loc[
        'mlst:lmonocytogenes:bglA:52'].tolist()
    assert ['lmonocytogenes', 'ldh', '5', 2, 1, 3, 4, 50, 25, 75] == summary_df.loc[
        'mlst:lmonocytogenes:ldh:5'].tolist()
    assert ['lmonocytogenes', 'lhkA', '5', 2, 0, 2, 4, 50, 0, 50] == summary_df.loc[
        'mlst:lmonocytogenes:lhkA:5'].tolist()
    assert ['lmonocytogenes', 'lhkA', '4', 1, 0, 1, 4, 25, 0, 25] == summary_df.loc[
        'mlst:lmonocytogenes:lhkA:4'].tolist()
    assert ['campylobacter', 'aspA', '2', 1, 0, 1, 4, 25, 0, 25] == summary_df.loc['mlst:campylobacter:aspA:2'].tolist()
    assert ['campylobacter', 'glyA', '3', 1, 0, 1, 4, 25, 0, 25] == summary_df.loc['mlst:campylobacter:glyA:3'].tolist()
    assert 6 == len(summary_df[summary_df['Scheme'] == 'campylobacter'])

    # Test multiple schemes, include unknown no present samples
    mlst_summarizer = MLSTFeaturesComparator(connection=loaded_database_genomic_data_store.connection,
                                             include_unknown_no_present_samples=True)
    present_set = SampleSet([sampleA.id, sampleB.id, sample_CFSAN002349.id, sample_2014D_0067.id])
    summary_df = mlst_summarizer.summary(present_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    summary_df['Unknown Percent'] = summary_df['Unknown Percent'].astype(int)
    summary_df['Present and Unknown Percent'] = summary_df['Present and Unknown Percent'].astype(int)
    assert 16 == len(summary_df)
    assert {'lmonocytogenes', 'campylobacter'} == set(summary_df['Scheme'].tolist())
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Unknown Count', 'Present and Unknown Count',
            'Total', 'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'abcZ', '1', 3, 0, 3, 4, 75, 0, 75] == summary_df.loc[
        'mlst:lmonocytogenes:abcZ:1'].tolist()
    assert ['lmonocytogenes', 'bglA', '51', 2, 0, 2, 4, 50, 0, 50] == summary_df.loc[
        'mlst:lmonocytogenes:bglA:51'].tolist()
    assert ['lmonocytogenes', 'bglA', '52', 1, 0, 1, 4, 25, 0, 25] == summary_df.loc[
        'mlst:lmonocytogenes:bglA:52'].tolist()
    assert ['lmonocytogenes', 'ldh', '5', 2, 1, 3, 4, 50, 25, 75] == summary_df.loc[
        'mlst:lmonocytogenes:ldh:5'].tolist()
    assert ['lmonocytogenes', 'lhkA', '5', 2, 0, 2, 4, 50, 0, 50] == summary_df.loc[
        'mlst:lmonocytogenes:lhkA:5'].tolist()
    assert ['lmonocytogenes', 'lhkA', '4', 1, 0, 1, 4, 25, 0, 25] == summary_df.loc[
        'mlst:lmonocytogenes:lhkA:4'].tolist()
    assert ['campylobacter', 'aspA', '2', 1, 0, 1, 4, 25, 0, 25] == summary_df.loc['mlst:campylobacter:aspA:2'].tolist()
    assert ['campylobacter', 'glyA', '3', 1, 0, 1, 4, 25, 0, 25] == summary_df.loc['mlst:campylobacter:glyA:3'].tolist()
    assert ['campylobacter', 'uncA', '6', 0, 1, 1, 4, 0, 25, 25] == summary_df.loc['mlst:campylobacter:uncA:6'].tolist()
    assert 7 == len(summary_df[summary_df['Scheme'] == 'campylobacter'])

    # Test multiple schemes sample set but summarize for only a particular scheme
    present_set = SampleSet([sampleA.id, sampleB.id, sample_CFSAN002349.id, sample_2014D_0067.id])
    mlst_summarizer = MLSTFeaturesComparator(connection=loaded_database_genomic_data_store.connection,
                                             scheme='lmonocytogenes',
                                             include_unknown_samples=True)
    summary_df = mlst_summarizer.summary(present_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    summary_df['Unknown Percent'] = summary_df['Unknown Percent'].astype(int)
    summary_df['Present and Unknown Percent'] = summary_df['Present and Unknown Percent'].astype(int)
    assert 9 == len(summary_df)
    assert {'lmonocytogenes'} == set(summary_df['Scheme'].tolist())  # Only results for one scheme
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Unknown Count', 'Present and Unknown Count',
            'Total', 'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'abcZ', '1', 3, 0, 3, 4, 75, 0, 75] == summary_df.loc[
        'mlst:lmonocytogenes:abcZ:1'].tolist()
    assert ['lmonocytogenes', 'bglA', '51', 2, 0, 2, 4, 50, 0, 50] == summary_df.loc[
        'mlst:lmonocytogenes:bglA:51'].tolist()
    assert ['lmonocytogenes', 'bglA', '52', 1, 0, 1, 4, 25, 0, 25] == summary_df.loc[
        'mlst:lmonocytogenes:bglA:52'].tolist()
    assert ['lmonocytogenes', 'ldh', '5', 2, 1, 3, 4, 50, 25, 75] == summary_df.loc[
        'mlst:lmonocytogenes:ldh:5'].tolist()
    assert ['lmonocytogenes', 'lhkA', '5', 2, 0, 2, 4, 50, 0, 50] == summary_df.loc[
        'mlst:lmonocytogenes:lhkA:5'].tolist()
    assert ['lmonocytogenes', 'lhkA', '4', 1, 0, 1, 4, 25, 0, 25] == summary_df.loc[
        'mlst:lmonocytogenes:lhkA:4'].tolist()

    # Test multiple schemes sample set but summarize for only a particular scheme/locus
    present_set = SampleSet([sampleA.id, sampleB.id, sample_CFSAN002349.id, sample_2014D_0067.id])
    mlst_summarizer = MLSTFeaturesComparator(connection=loaded_database_genomic_data_store.connection,
                                             scheme='lmonocytogenes', locus='bglA',
                                             include_unknown_samples=True)
    summary_df = mlst_summarizer.summary(present_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    summary_df['Unknown Percent'] = summary_df['Unknown Percent'].astype(int)
    summary_df['Present and Unknown Percent'] = summary_df['Present and Unknown Percent'].astype(int)
    assert 2 == len(summary_df)
    assert {'lmonocytogenes'} == set(summary_df['Scheme'].tolist())  # Only results for one scheme
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Unknown Count', 'Present and Unknown Count',
            'Total', 'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'bglA', '51', 2, 0, 2, 4, 50, 0, 50] == summary_df.loc[
        'mlst:lmonocytogenes:bglA:51'].tolist()
    assert ['lmonocytogenes', 'bglA', '52', 1, 0, 1, 4, 25, 0, 25] == summary_df.loc[
        'mlst:lmonocytogenes:bglA:52'].tolist()

    # Test multiple schemes, include unknown, include unknown no present samples
    present_set = SampleSet([sampleA.id, sampleB.id, sample_CFSAN002349.id, sample_2014D_0067.id])
    mlst_summarizer = MLSTFeaturesComparator(connection=loaded_database_genomic_data_store.connection,
                                             include_unknown=True, include_unknown_no_present_samples=True)
    summary_df = mlst_summarizer.summary(present_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    summary_df['Unknown Count'] = summary_df['Unknown Count'].fillna(-1).astype(int)
    summary_df['Present and Unknown Count'] = summary_df['Present and Unknown Count'].fillna(-1).astype(int)
    summary_df['Unknown Percent'] = summary_df['Unknown Percent'].fillna(-1).astype(int)
    summary_df['Present and Unknown Percent'] = summary_df['Present and Unknown Percent'].fillna(-1).astype(int)
    assert 18 == len(summary_df)
    assert {'lmonocytogenes', 'campylobacter'} == set(summary_df['Scheme'].tolist())
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Unknown Count', 'Present and Unknown Count',
            'Total', 'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'abcZ', '1', 3, 0, 3, 4, 75, 0, 75] == summary_df.loc[
        'mlst:lmonocytogenes:abcZ:1'].tolist()
    assert ['lmonocytogenes', 'bglA', '51', 2, 0, 2, 4, 50, 0, 50] == summary_df.loc[
        'mlst:lmonocytogenes:bglA:51'].tolist()
    assert ['lmonocytogenes', 'bglA', '52', 1, 0, 1, 4, 25, 0, 25] == summary_df.loc[
        'mlst:lmonocytogenes:bglA:52'].tolist()
    assert ['lmonocytogenes', 'ldh', '5', 2, 1, 3, 4, 50, 25, 75] == summary_df.loc[
        'mlst:lmonocytogenes:ldh:5'].tolist()
    assert ['lmonocytogenes', 'ldh', '?', 1, -1, -1, 4, 25, -1, -1] == summary_df.loc[
        'mlst:lmonocytogenes:ldh:?'].tolist()
    assert ['lmonocytogenes', 'lhkA', '5', 2, 0, 2, 4, 50, 0, 50] == summary_df.loc[
        'mlst:lmonocytogenes:lhkA:5'].tolist()
    assert ['lmonocytogenes', 'lhkA', '4', 1, 0, 1, 4, 25, 0, 25] == summary_df.loc[
        'mlst:lmonocytogenes:lhkA:4'].tolist()
    assert ['campylobacter', 'aspA', '2', 1, 0, 1, 4, 25, 0, 25] == summary_df.loc['mlst:campylobacter:aspA:2'].tolist()
    assert ['campylobacter', 'glyA', '3', 1, 0, 1, 4, 25, 0, 25] == summary_df.loc['mlst:campylobacter:glyA:3'].tolist()
    assert ['campylobacter', 'uncA', '6', 0, 1, 1, 4, 0, 25, 25] == summary_df.loc['mlst:campylobacter:uncA:6'].tolist()
    assert ['campylobacter', 'uncA', '?', 1, -1, -1, 4, 25, -1, -1] == summary_df.loc[
        'mlst:campylobacter:uncA:?'].tolist()

    # Test multiple schemes, include unknown, no include unknown no present samples
    present_set = SampleSet([sampleA.id, sampleB.id, sample_CFSAN002349.id, sample_2014D_0067.id])
    mlst_summarizer = MLSTFeaturesComparator(connection=loaded_database_genomic_data_store.connection,
                                             include_unknown=True, include_unknown_no_present_samples=False)
    summary_df = mlst_summarizer.summary(present_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    summary_df['Unknown Count'] = summary_df['Unknown Count'].fillna(-1).astype(int)
    summary_df['Present and Unknown Count'] = summary_df['Present and Unknown Count'].fillna(-1).astype(int)
    summary_df['Unknown Percent'] = summary_df['Unknown Percent'].fillna(-1).astype(int)
    summary_df['Present and Unknown Percent'] = summary_df['Present and Unknown Percent'].fillna(-1).astype(int)
    assert 17 == len(summary_df)
    assert {'lmonocytogenes', 'campylobacter'} == set(summary_df['Scheme'].tolist())
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Unknown Count', 'Present and Unknown Count',
            'Total', 'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'abcZ', '1', 3, 0, 3, 4, 75, 0, 75] == summary_df.loc[
        'mlst:lmonocytogenes:abcZ:1'].tolist()
    assert ['lmonocytogenes', 'bglA', '51', 2, 0, 2, 4, 50, 0, 50] == summary_df.loc[
        'mlst:lmonocytogenes:bglA:51'].tolist()
    assert ['lmonocytogenes', 'bglA', '52', 1, 0, 1, 4, 25, 0, 25] == summary_df.loc[
        'mlst:lmonocytogenes:bglA:52'].tolist()
    assert ['lmonocytogenes', 'ldh', '5', 2, 1, 3, 4, 50, 25, 75] == summary_df.loc[
        'mlst:lmonocytogenes:ldh:5'].tolist()
    assert ['lmonocytogenes', 'ldh', '?', 1, -1, -1, 4, 25, -1, -1] == summary_df.loc[
        'mlst:lmonocytogenes:ldh:?'].tolist()
    assert ['lmonocytogenes', 'lhkA', '5', 2, 0, 2, 4, 50, 0, 50] == summary_df.loc[
        'mlst:lmonocytogenes:lhkA:5'].tolist()
    assert ['lmonocytogenes', 'lhkA', '4', 1, 0, 1, 4, 25, 0, 25] == summary_df.loc[
        'mlst:lmonocytogenes:lhkA:4'].tolist()
    assert ['campylobacter', 'aspA', '2', 1, 0, 1, 4, 25, 0, 25] == summary_df.loc['mlst:campylobacter:aspA:2'].tolist()
    assert ['campylobacter', 'glyA', '3', 1, 0, 1, 4, 25, 0, 25] == summary_df.loc['mlst:campylobacter:glyA:3'].tolist()
    assert ['campylobacter', 'uncA', '?', 1, -1, -1, 4, 25, -1, -1] == summary_df.loc[
        'mlst:campylobacter:uncA:?'].tolist()
    assert 'mlst:campylobacter:uncA:6' not in summary_df

    # Test multiple schemes, only unknown
    present_set = SampleSet([sampleA.id, sampleB.id, sample_CFSAN002349.id, sample_2014D_0067.id])
    mlst_summarizer = MLSTFeaturesComparator(connection=loaded_database_genomic_data_store.connection,
                                             include_present=False, include_unknown=True,
                                             include_unknown_samples=True)
    summary_df = mlst_summarizer.summary(present_set)

    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    summary_df['Unknown Count'] = summary_df['Unknown Count'].fillna(-1)
    summary_df['Unknown Count'] = summary_df['Unknown Count'].astype(int)
    summary_df['Present and Unknown Count'] = summary_df['Present and Unknown Count'].fillna(-1)
    summary_df['Present and Unknown Count'] = summary_df['Present and Unknown Count'].astype(int)
    summary_df['Unknown Percent'] = summary_df['Unknown Percent'].fillna(-1)
    summary_df['Unknown Percent'] = summary_df['Unknown Percent'].astype(int)
    summary_df['Present and Unknown Percent'] = summary_df['Present and Unknown Percent'].fillna(-1)
    summary_df['Present and Unknown Percent'] = summary_df['Present and Unknown Percent'].astype(int)
    assert 2 == len(summary_df)
    assert {'lmonocytogenes', 'campylobacter'} == set(summary_df['Scheme'].tolist())
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Unknown Count', 'Present and Unknown Count',
            'Total', 'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'ldh', '?', 1, -1, -1, 4, 25, -1, -1] == summary_df.loc[
        'mlst:lmonocytogenes:ldh:?'].tolist()
    assert ['campylobacter', 'uncA', '?', 1, -1, -1, 4, 25, -1, -1] == summary_df.loc[
        'mlst:campylobacter:uncA:?'].tolist()


def test_features_comparison(loaded_database_genomic_data_store: GenomicsDataIndex):
    db = loaded_database_genomic_data_store.connection.database
    all_sample_ids = {s.id for s in db.get_session().query(Sample).all()}
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()
    sample_CFSAN002349 = db.get_session().query(Sample).filter(Sample.name == 'CFSAN002349').one()
    sample_CFSAN023463 = db.get_session().query(Sample).filter(Sample.name == 'CFSAN023463').one()
    lmonocytogenes = {sampleA.id, sampleB.id, sampleC.id, sample_CFSAN002349.id, sample_CFSAN023463.id}
    assert 9 == len(all_sample_ids)

    present_set = SampleSet(all_sample_ids)
    mlst_summarizer = MLSTFeaturesComparator(
        connection=loaded_database_genomic_data_store.connection)

    # Test single category of all
    sample_categories = [present_set]
    comparison_df = mlst_summarizer.features_comparison(selected_samples=present_set,
                                                        sample_categories=sample_categories,
                                                        category_prefixes=['All'],
                                                        unit='count')
    assert 24 == len(comparison_df)
    assert 'MLST Feature' == comparison_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Total', 'All_count', 'All_Unknown count',
            'All_Present and Unknown count', 'All_total'] == list(comparison_df.columns)
    assert {9} == set(comparison_df['Total'].tolist())
    assert {9} == set(comparison_df['All_total'].tolist())
    assert 5 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'All_count']
    assert 3 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'All_count']
    assert 2 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'All_count']
    assert 2 == comparison_df.loc['mlst:ecoli:adk:100', 'All_count']
    assert 2 == comparison_df.loc['mlst:ecoli:recA:7', 'All_count']
    assert 1 == comparison_df.loc['mlst:campylobacter:uncA:6', 'All_count']

    assert 0 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'All_Unknown count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'All_Unknown count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'All_Unknown count']
    assert 1 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'All_Unknown count']
    assert 1 == comparison_df.loc['mlst:campylobacter:uncA:6', 'All_Unknown count']

    assert 5 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'All_Present and Unknown count']
    assert 3 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'All_Present and Unknown count']
    assert 2 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'All_Present and Unknown count']
    assert 5 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'All_Present and Unknown count']
    assert 2 == comparison_df.loc['mlst:campylobacter:uncA:6', 'All_Present and Unknown count']

    # Test two categories: one of lmonocytogenes and one of the rest
    sample_categories = [SampleSet(lmonocytogenes), SampleSet(all_sample_ids - lmonocytogenes)]
    comparison_df = mlst_summarizer.features_comparison(selected_samples=present_set,
                                                        sample_categories=sample_categories,
                                                        category_prefixes=['lmonocytogenes', 'other'],
                                                        unit='count')
    assert 24 == len(comparison_df)
    assert 'MLST Feature' == comparison_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Total',
            'lmonocytogenes_count', 'other_count',
            'lmonocytogenes_Unknown count', 'other_Unknown count',
            'lmonocytogenes_Present and Unknown count', 'other_Present and Unknown count',
            'lmonocytogenes_total', 'other_total'] == list(comparison_df.columns)
    assert {9} == set(comparison_df['Total'].tolist())
    assert {5} == set(comparison_df['lmonocytogenes_total'].tolist())
    assert {4} == set(comparison_df['other_total'].tolist())
    assert 5 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'lmonocytogenes_count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'other_count']
    assert 3 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'lmonocytogenes_count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'other_count']
    assert 2 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'lmonocytogenes_count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'other_count']
    assert 0 == comparison_df.loc['mlst:ecoli:adk:100', 'lmonocytogenes_count']
    assert 2 == comparison_df.loc['mlst:ecoli:adk:100', 'other_count']
    assert 0 == comparison_df.loc['mlst:ecoli:recA:7', 'lmonocytogenes_count']
    assert 2 == comparison_df.loc['mlst:ecoli:recA:7', 'other_count']
    assert 4 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'lmonocytogenes_count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'other_count']
    assert 0 == comparison_df.loc['mlst:campylobacter:uncA:6', 'lmonocytogenes_count']
    assert 1 == comparison_df.loc['mlst:campylobacter:uncA:6', 'other_count']

    assert 0 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'lmonocytogenes_Unknown count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'other_Unknown count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'lmonocytogenes_Unknown count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'other_Unknown count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'lmonocytogenes_Unknown count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'other_Unknown count']
    assert 0 == comparison_df.loc['mlst:ecoli:adk:100', 'lmonocytogenes_Unknown count']
    assert 0 == comparison_df.loc['mlst:ecoli:adk:100', 'other_Unknown count']
    assert 0 == comparison_df.loc['mlst:ecoli:recA:7', 'lmonocytogenes_Unknown count']
    assert 0 == comparison_df.loc['mlst:ecoli:recA:7', 'other_Unknown count']
    assert 1 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'lmonocytogenes_Unknown count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'other_Unknown count']
    assert 0 == comparison_df.loc['mlst:campylobacter:uncA:6', 'lmonocytogenes_Unknown count']
    assert 1 == comparison_df.loc['mlst:campylobacter:uncA:6', 'other_Unknown count']

    assert 5 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'lmonocytogenes_Present and Unknown count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'other_Present and Unknown count']
    assert 3 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'lmonocytogenes_Present and Unknown count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'other_Present and Unknown count']
    assert 2 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'lmonocytogenes_Present and Unknown count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'other_Present and Unknown count']
    assert 0 == comparison_df.loc['mlst:ecoli:adk:100', 'lmonocytogenes_Present and Unknown count']
    assert 2 == comparison_df.loc['mlst:ecoli:adk:100', 'other_Present and Unknown count']
    assert 0 == comparison_df.loc['mlst:ecoli:recA:7', 'lmonocytogenes_Present and Unknown count']
    assert 2 == comparison_df.loc['mlst:ecoli:recA:7', 'other_Present and Unknown count']
    assert 5 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'lmonocytogenes_Present and Unknown count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'other_Present and Unknown count']
    assert 0 == comparison_df.loc['mlst:campylobacter:uncA:6', 'lmonocytogenes_Present and Unknown count']
    assert 2 == comparison_df.loc['mlst:campylobacter:uncA:6', 'other_Present and Unknown count']

    # Test two categories percent: one of lmonocytogenes and one of the rest
    sample_categories = [SampleSet(lmonocytogenes), SampleSet(all_sample_ids - lmonocytogenes)]
    comparison_df = mlst_summarizer.features_comparison(selected_samples=present_set,
                                                        sample_categories=sample_categories,
                                                        category_prefixes=['lmonocytogenes', 'other'],
                                                        unit='percent')
    assert 24 == len(comparison_df)
    assert 'MLST Feature' == comparison_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Total',
            'lmonocytogenes_percent', 'other_percent',
            'lmonocytogenes_Unknown percent', 'other_Unknown percent',
            'lmonocytogenes_Present and Unknown percent', 'other_Present and Unknown percent',
            'lmonocytogenes_total', 'other_total'] == list(comparison_df.columns)
    comparison_df['lmonocytogenes_percent'] = comparison_df['lmonocytogenes_percent'].astype(
        int)  # Convert to int for easier comparison
    comparison_df['other_percent'] = comparison_df['other_percent'].astype(int)  # Convert to int for easier comparison
    assert {9} == set(comparison_df['Total'].tolist())
    assert {5} == set(comparison_df['lmonocytogenes_total'].tolist())
    assert {4} == set(comparison_df['other_total'].tolist())
    assert 100 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'lmonocytogenes_percent']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'other_percent']
    assert 60 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'lmonocytogenes_percent']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'other_percent']
    assert 40 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'lmonocytogenes_percent']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'other_percent']
    assert 0 == comparison_df.loc['mlst:ecoli:adk:100', 'lmonocytogenes_percent']
    assert 50 == comparison_df.loc['mlst:ecoli:adk:100', 'other_percent']
    assert 0 == comparison_df.loc['mlst:ecoli:recA:7', 'lmonocytogenes_percent']
    assert 50 == comparison_df.loc['mlst:ecoli:recA:7', 'other_percent']
    assert 80 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'lmonocytogenes_percent']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'other_percent']
    assert 0 == comparison_df.loc['mlst:campylobacter:uncA:6', 'lmonocytogenes_percent']
    assert 25 == comparison_df.loc['mlst:campylobacter:uncA:6', 'other_percent']

    assert 0 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'lmonocytogenes_Unknown percent']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'other_Unknown percent']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'lmonocytogenes_Unknown percent']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'other_Unknown percent']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'lmonocytogenes_Unknown percent']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'other_Unknown percent']
    assert 0 == comparison_df.loc['mlst:ecoli:adk:100', 'lmonocytogenes_Unknown percent']
    assert 0 == comparison_df.loc['mlst:ecoli:adk:100', 'other_Unknown percent']
    assert 0 == comparison_df.loc['mlst:ecoli:recA:7', 'lmonocytogenes_Unknown percent']
    assert 0 == comparison_df.loc['mlst:ecoli:recA:7', 'other_Unknown percent']
    assert 20 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'lmonocytogenes_Unknown percent']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'other_Unknown percent']
    assert 0 == comparison_df.loc['mlst:campylobacter:uncA:6', 'lmonocytogenes_Unknown percent']
    assert 25 == comparison_df.loc['mlst:campylobacter:uncA:6', 'other_Unknown percent']

    assert 100 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'lmonocytogenes_Present and Unknown percent']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'other_Present and Unknown percent']
    assert 60 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'lmonocytogenes_Present and Unknown percent']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'other_Present and Unknown percent']
    assert 40 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'lmonocytogenes_Present and Unknown percent']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'other_Present and Unknown percent']
    assert 0 == comparison_df.loc['mlst:ecoli:adk:100', 'lmonocytogenes_Present and Unknown percent']
    assert 50 == comparison_df.loc['mlst:ecoli:adk:100', 'other_Present and Unknown percent']
    assert 0 == comparison_df.loc['mlst:ecoli:recA:7', 'lmonocytogenes_Present and Unknown percent']
    assert 50 == comparison_df.loc['mlst:ecoli:recA:7', 'other_Present and Unknown percent']
    assert 100 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'lmonocytogenes_Present and Unknown percent']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'other_Present and Unknown percent']
    assert 0 == comparison_df.loc['mlst:campylobacter:uncA:6', 'lmonocytogenes_Present and Unknown percent']
    assert 50 == comparison_df.loc['mlst:campylobacter:uncA:6', 'other_Present and Unknown percent']

    # Test two categories proportion: one of lmonocytogenes and one of the rest
    sample_categories = [SampleSet(lmonocytogenes), SampleSet(all_sample_ids - lmonocytogenes)]
    comparison_df = mlst_summarizer.features_comparison(selected_samples=present_set,
                                                        sample_categories=sample_categories,
                                                        category_prefixes=['lmonocytogenes', 'other'],
                                                        unit='proportion')
    assert 24 == len(comparison_df)
    assert 'MLST Feature' == comparison_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Total',
            'lmonocytogenes_proportion', 'other_proportion',
            'lmonocytogenes_Unknown proportion', 'other_Unknown proportion',
            'lmonocytogenes_Present and Unknown proportion', 'other_Present and Unknown proportion',
            'lmonocytogenes_total', 'other_total'] == list(comparison_df.columns)
    # Convert to percent as int for easier comparison
    comparison_df['lmonocytogenes_proportion'] = (comparison_df['lmonocytogenes_proportion'] * 100).astype(int)
    comparison_df['lmonocytogenes_Unknown proportion'] = (
                comparison_df['lmonocytogenes_Unknown proportion'] * 100).astype(int)
    comparison_df['lmonocytogenes_Present and Unknown proportion'] = (
                comparison_df['lmonocytogenes_Present and Unknown proportion'] * 100).astype(int)
    comparison_df['other_proportion'] = (comparison_df['other_proportion'] * 100).astype(int)
    comparison_df['other_Unknown proportion'] = (comparison_df['other_Unknown proportion'] * 100).astype(int)
    comparison_df['other_Present and Unknown proportion'] = (
                comparison_df['other_Present and Unknown proportion'] * 100).astype(int)
    assert {9} == set(comparison_df['Total'].tolist())
    assert {5} == set(comparison_df['lmonocytogenes_total'].tolist())
    assert {4} == set(comparison_df['other_total'].tolist())
    assert 100 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'lmonocytogenes_proportion']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'other_proportion']
    assert 60 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'lmonocytogenes_proportion']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'other_proportion']
    assert 40 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'lmonocytogenes_proportion']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'other_proportion']
    assert 0 == comparison_df.loc['mlst:ecoli:adk:100', 'lmonocytogenes_proportion']
    assert 50 == comparison_df.loc['mlst:ecoli:adk:100', 'other_proportion']
    assert 0 == comparison_df.loc['mlst:ecoli:recA:7', 'lmonocytogenes_proportion']
    assert 50 == comparison_df.loc['mlst:ecoli:recA:7', 'other_proportion']
    assert 80 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'lmonocytogenes_proportion']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'other_proportion']
    assert 0 == comparison_df.loc['mlst:campylobacter:uncA:6', 'lmonocytogenes_proportion']
    assert 25 == comparison_df.loc['mlst:campylobacter:uncA:6', 'other_proportion']

    assert 0 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'lmonocytogenes_Unknown proportion']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'other_Unknown proportion']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'lmonocytogenes_Unknown proportion']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'other_Unknown proportion']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'lmonocytogenes_Unknown proportion']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'other_Unknown proportion']
    assert 0 == comparison_df.loc['mlst:ecoli:adk:100', 'lmonocytogenes_Unknown proportion']
    assert 0 == comparison_df.loc['mlst:ecoli:adk:100', 'other_Unknown proportion']
    assert 0 == comparison_df.loc['mlst:ecoli:recA:7', 'lmonocytogenes_Unknown proportion']
    assert 0 == comparison_df.loc['mlst:ecoli:recA:7', 'other_Unknown proportion']
    assert 20 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'lmonocytogenes_Unknown proportion']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'other_Unknown proportion']
    assert 0 == comparison_df.loc['mlst:campylobacter:uncA:6', 'lmonocytogenes_Unknown proportion']
    assert 25 == comparison_df.loc['mlst:campylobacter:uncA:6', 'other_Unknown proportion']

    assert 100 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'lmonocytogenes_Present and Unknown proportion']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'other_Present and Unknown proportion']
    assert 60 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'lmonocytogenes_Present and Unknown proportion']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'other_Present and Unknown proportion']
    assert 40 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'lmonocytogenes_Present and Unknown proportion']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'other_Present and Unknown proportion']
    assert 0 == comparison_df.loc['mlst:ecoli:adk:100', 'lmonocytogenes_Present and Unknown proportion']
    assert 50 == comparison_df.loc['mlst:ecoli:adk:100', 'other_Present and Unknown proportion']
    assert 0 == comparison_df.loc['mlst:ecoli:recA:7', 'lmonocytogenes_Present and Unknown proportion']
    assert 50 == comparison_df.loc['mlst:ecoli:recA:7', 'other_Present and Unknown proportion']
    assert 100 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'lmonocytogenes_Present and Unknown proportion']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'other_Present and Unknown proportion']
    assert 0 == comparison_df.loc['mlst:campylobacter:uncA:6', 'lmonocytogenes_Present and Unknown proportion']
    assert 50 == comparison_df.loc['mlst:campylobacter:uncA:6', 'other_Present and Unknown proportion']

    # Test two categories: one of lmonocytogenes and one of the rest threshold below
    sample_categories = [SampleSet(lmonocytogenes), SampleSet(all_sample_ids - lmonocytogenes)]
    comparison_df = mlst_summarizer.features_comparison(selected_samples=present_set,
                                                        sample_categories=sample_categories,
                                                        category_prefixes=['lmonocytogenes', 'other'],
                                                        category_samples_threshold=4,
                                                        unit='count')
    assert 24 == len(comparison_df)
    assert 'MLST Feature' == comparison_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Total',
            'lmonocytogenes_count', 'other_count',
            'lmonocytogenes_Unknown count', 'other_Unknown count',
            'lmonocytogenes_Present and Unknown count', 'other_Present and Unknown count',
            'lmonocytogenes_total', 'other_total'] == list(comparison_df.columns)
    assert {9} == set(comparison_df['Total'].tolist())
    assert {5} == set(comparison_df['lmonocytogenes_total'].tolist())
    assert {4} == set(comparison_df['other_total'].tolist())
    assert 5 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'lmonocytogenes_count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'other_count']
    assert 3 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'lmonocytogenes_count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'other_count']
    assert 2 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'lmonocytogenes_count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'other_count']
    assert 0 == comparison_df.loc['mlst:ecoli:adk:100', 'lmonocytogenes_count']
    assert 2 == comparison_df.loc['mlst:ecoli:adk:100', 'other_count']
    assert 0 == comparison_df.loc['mlst:ecoli:recA:7', 'lmonocytogenes_count']
    assert 2 == comparison_df.loc['mlst:ecoli:recA:7', 'other_count']
    assert 4 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'lmonocytogenes_count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'other_count']
    assert 0 == comparison_df.loc['mlst:campylobacter:uncA:6', 'lmonocytogenes_count']
    assert 1 == comparison_df.loc['mlst:campylobacter:uncA:6', 'other_count']

    assert 0 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'lmonocytogenes_Unknown count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'other_Unknown count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'lmonocytogenes_Unknown count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'other_Unknown count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'lmonocytogenes_Unknown count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'other_Unknown count']
    assert 0 == comparison_df.loc['mlst:ecoli:adk:100', 'lmonocytogenes_Unknown count']
    assert 0 == comparison_df.loc['mlst:ecoli:adk:100', 'other_Unknown count']
    assert 0 == comparison_df.loc['mlst:ecoli:recA:7', 'lmonocytogenes_Unknown count']
    assert 0 == comparison_df.loc['mlst:ecoli:recA:7', 'other_Unknown count']
    assert 1 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'lmonocytogenes_Unknown count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'other_Unknown count']
    assert 0 == comparison_df.loc['mlst:campylobacter:uncA:6', 'lmonocytogenes_Unknown count']
    assert 1 == comparison_df.loc['mlst:campylobacter:uncA:6', 'other_Unknown count']

    assert 5 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'lmonocytogenes_Present and Unknown count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'other_Present and Unknown count']
    assert 3 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'lmonocytogenes_Present and Unknown count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'other_Present and Unknown count']
    assert 2 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'lmonocytogenes_Present and Unknown count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'other_Present and Unknown count']
    assert 0 == comparison_df.loc['mlst:ecoli:adk:100', 'lmonocytogenes_Present and Unknown count']
    assert 2 == comparison_df.loc['mlst:ecoli:adk:100', 'other_Present and Unknown count']
    assert 0 == comparison_df.loc['mlst:ecoli:recA:7', 'lmonocytogenes_Present and Unknown count']
    assert 2 == comparison_df.loc['mlst:ecoli:recA:7', 'other_Present and Unknown count']
    assert 5 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'lmonocytogenes_Present and Unknown count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'other_Present and Unknown count']
    assert 0 == comparison_df.loc['mlst:campylobacter:uncA:6', 'lmonocytogenes_Present and Unknown count']
    assert 2 == comparison_df.loc['mlst:campylobacter:uncA:6', 'other_Present and Unknown count']

    # Test two categories: one of lmonocytogenes and one of the rest threshold above
    sample_categories = [SampleSet(lmonocytogenes), SampleSet(all_sample_ids - lmonocytogenes)]
    comparison_df = mlst_summarizer.features_comparison(selected_samples=present_set,
                                                        sample_categories=sample_categories,
                                                        category_prefixes=['lmonocytogenes', 'other'],
                                                        category_samples_threshold=5,
                                                        unit='count')
    assert 24 == len(comparison_df)
    assert 'MLST Feature' == comparison_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Total',
            'lmonocytogenes_count', 'lmonocytogenes_Unknown count', 'lmonocytogenes_Present and Unknown count',
            'lmonocytogenes_total'] == list(comparison_df.columns)
    assert {9} == set(comparison_df['Total'].tolist())
    assert {5} == set(comparison_df['lmonocytogenes_total'].tolist())
    assert 5 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'lmonocytogenes_count']
    assert 3 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'lmonocytogenes_count']
    assert 2 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'lmonocytogenes_count']
    assert 0 == comparison_df.loc['mlst:ecoli:adk:100', 'lmonocytogenes_count']
    assert 0 == comparison_df.loc['mlst:ecoli:recA:7', 'lmonocytogenes_count']
    assert 4 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'lmonocytogenes_count']
    assert 0 == comparison_df.loc['mlst:campylobacter:uncA:6', 'lmonocytogenes_count']

    assert 0 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'lmonocytogenes_Unknown count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'lmonocytogenes_Unknown count']
    assert 0 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'lmonocytogenes_Unknown count']
    assert 0 == comparison_df.loc['mlst:ecoli:adk:100', 'lmonocytogenes_Unknown count']
    assert 0 == comparison_df.loc['mlst:ecoli:recA:7', 'lmonocytogenes_Unknown count']
    assert 1 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'lmonocytogenes_Unknown count']
    assert 0 == comparison_df.loc['mlst:campylobacter:uncA:6', 'lmonocytogenes_Unknown count']

    assert 5 == comparison_df.loc['mlst:lmonocytogenes:abcZ:1', 'lmonocytogenes_Present and Unknown count']
    assert 3 == comparison_df.loc['mlst:lmonocytogenes:bglA:51', 'lmonocytogenes_Present and Unknown count']
    assert 2 == comparison_df.loc['mlst:lmonocytogenes:bglA:52', 'lmonocytogenes_Present and Unknown count']
    assert 0 == comparison_df.loc['mlst:ecoli:adk:100', 'lmonocytogenes_Present and Unknown count']
    assert 0 == comparison_df.loc['mlst:ecoli:recA:7', 'lmonocytogenes_Present and Unknown count']
    assert 5 == comparison_df.loc['mlst:lmonocytogenes:ldh:5', 'lmonocytogenes_Present and Unknown count']
    assert 0 == comparison_df.loc['mlst:campylobacter:uncA:6', 'lmonocytogenes_Present and Unknown count']
