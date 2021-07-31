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

