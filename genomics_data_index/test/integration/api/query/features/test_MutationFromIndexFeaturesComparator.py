from pathlib import Path
from typing import List, Union

import pandas as pd

from genomics_data_index.api.query.GenomicsDataIndex import GenomicsDataIndex
from genomics_data_index.api.query.features.MutationFeaturesFromIndexComparator import \
    MutationFeaturesFromIndexComparator
from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.db import Sample
from genomics_data_index.test.integration import snippy_all_dataframes


def read_expected_snippy_df(snippy_mutations: Union[Path, List[Path]], total: int,
                            mutations_not_in: Union[Path, List[Path]] = None) -> pd.DataFrame:
    if isinstance(snippy_mutations, Path):
        snippy_mutations = [snippy_mutations]

    if isinstance(mutations_not_in, Path):
        mutations_not_in = [mutations_not_in]

    snippy_dfs = [pd.read_csv(p, sep='\t') for p in snippy_mutations]
    expected_df = pd.concat(snippy_dfs)
    expected_df = expected_df.groupby('Mutation').agg({
        'Sequence': 'first',
        'Position': 'first',
        'Deletion': 'first',
        'Insertion': 'first',
        'Type': 'first',
        'Mutation': 'count',
    }).rename(columns={'Mutation': 'Count'}).sort_index()
    expected_df['Unknown Count'] = '<NA>'
    expected_df['Present and Unknown Count'] = '<NA>'
    expected_df['Total'] = total
    expected_df['Percent'] = 100 * (expected_df['Count'] / expected_df['Total'])
    expected_df['Unknown Percent'] = '<NA>'
    expected_df['Present and Unknown Percent'] = '<NA>'

    if mutations_not_in is not None:
        notin_dfs = pd.concat([pd.read_csv(p, sep='\t') for p in mutations_not_in])
        expected_df = expected_df.loc[~expected_df.index.isin(list(notin_dfs['Mutation']))]

    return expected_df


def test_summary_all(loaded_database_genomic_data_store: GenomicsDataIndex):
    db = loaded_database_genomic_data_store.connection.database
    all_sample_ids = {s.id for s in db.get_session().query(Sample).all()}

    expected_df = read_expected_snippy_df(list(snippy_all_dataframes.values()), total=9)

    present_set = SampleSet(all_sample_ids)
    mutations_summarizer = MutationFeaturesFromIndexComparator(connection=loaded_database_genomic_data_store.connection,
                                                               include_unknown_samples=False,
                                                               ignore_annotations=True)

    mutations_df = mutations_summarizer.summary(present_set)
    mutations_df['Percent'] = mutations_df['Percent'].astype(int)  # Convert to int for easier comparison
    mutations_df['Unknown Percent'] = mutations_df['Unknown Percent'].fillna('<NA>')
    mutations_df['Unknown Count'] = mutations_df['Unknown Count'].fillna('<NA>')
    mutations_df['Present and Unknown Percent'] = mutations_df['Present and Unknown Percent'].fillna('<NA>')
    mutations_df['Present and Unknown Count'] = mutations_df['Present and Unknown Count'].fillna('<NA>')
    mutations_df = mutations_df.sort_index()

    assert len(expected_df) == len(mutations_df)
    assert list(expected_df.columns) == list(mutations_df.columns)
    assert list(expected_df.index) == list(mutations_df.index)
    assert list(expected_df['Deletion']) == list(mutations_df['Deletion'])
    assert list(expected_df['Count']) == list(mutations_df['Count'])
    assert list(expected_df['Unknown Count']) == list(mutations_df['Unknown Count'])
    assert list(expected_df['Present and Unknown Count']) == list(mutations_df['Present and Unknown Count'])
    assert list(expected_df['Total']) == list(mutations_df['Total'])
    assert list(expected_df['Type']) == list(mutations_df['Type'])
    assert 22 == mutations_df.loc['reference:619:G:C', 'Percent']
    assert '<NA>' == mutations_df.loc['reference:619:G:C', 'Unknown Percent']
    assert '<NA>' == mutations_df.loc['reference:619:G:C', 'Present and Unknown Percent']
    assert 11 == mutations_df.loc['reference:461:AAAT:G', 'Percent']
    assert '<NA>' == mutations_df.loc['reference:461:AAAT:G', 'Unknown Percent']
    assert '<NA>' == mutations_df.loc['reference:461:AAAT:G', 'Present and Unknown Percent']

    # Test with unknown features
    mutations_summarizer = MutationFeaturesFromIndexComparator(connection=loaded_database_genomic_data_store.connection,
                                                               include_unknown_samples=False,
                                                               ignore_annotations=True, include_unknown=True)
    mutations_df = mutations_summarizer.summary(present_set)
    mutations_df['Percent'] = mutations_df['Percent'].astype(int)  # Convert to int for easier comparison
    mutations_df['Unknown Percent'] = mutations_df['Unknown Percent'].fillna('<NA>')
    mutations_df['Unknown Count'] = mutations_df['Unknown Count'].fillna('<NA>')
    mutations_df = mutations_df.sort_index()

    assert 112 + 440 == len(mutations_df)
    assert list(expected_df.columns) == list(mutations_df.columns)
    assert 2 == mutations_df.loc['reference:619:G:C', 'Count']
    assert 'SNP' == mutations_df.loc['reference:619:G:C', 'Type']
    assert 2 == mutations_df.loc['reference:3063:A:ATGCAGC', 'Count']
    assert 'INDEL' == mutations_df.loc['reference:3063:A:ATGCAGC', 'Type']
    assert 1 == mutations_df.loc['reference:1984:GTGATTG:TTGA', 'Count']
    assert 'OTHER' == mutations_df.loc['reference:1984:GTGATTG:TTGA', 'Type']
    assert 1 == mutations_df.loc['reference:866:GCCAGATCC:G', 'Count']
    assert 'INDEL' == mutations_df.loc['reference:866:GCCAGATCC:G', 'Type']
    assert 3 == mutations_df.loc['reference:90:T:?', 'Count']
    assert 'UNKNOWN_MISSING' == mutations_df.loc['reference:90:T:?', 'Type']
    assert 2 == mutations_df.loc['reference:190:A:?', 'Count']
    assert 'UNKNOWN_MISSING' == mutations_df.loc['reference:190:A:?', 'Type']
    assert 1 == mutations_df.loc['reference:210:C:?', 'Count']
    assert 'UNKNOWN_MISSING' == mutations_df.loc['reference:210:C:?', 'Type']

    # Test only include unknown features
    mutations_summarizer = MutationFeaturesFromIndexComparator(connection=loaded_database_genomic_data_store.connection,
                                                               ignore_annotations=True, include_unknown=True,
                                                               include_unknown_samples=False,
                                                               include_present=False)
    mutations_df = mutations_summarizer.summary(present_set)
    mutations_df['Percent'] = mutations_df['Percent'].astype(int)  # Convert to int for easier comparison
    mutations_df['Unknown Percent'] = mutations_df['Unknown Percent'].fillna('<NA>')
    mutations_df['Unknown Count'] = mutations_df['Unknown Count'].fillna('<NA>')
    mutations_df['Present and Unknown Percent'] = mutations_df['Present and Unknown Percent'].fillna('<NA>')
    mutations_df['Present and Unknown Count'] = mutations_df['Present and Unknown Count'].fillna('<NA>')
    mutations_df = mutations_df.sort_index()

    assert 440 == len(mutations_df)
    assert list(expected_df.columns) == list(mutations_df.columns)
    assert 3 == mutations_df.loc['reference:90:T:?', 'Count']
    assert 'UNKNOWN_MISSING' == mutations_df.loc['reference:90:T:?', 'Type']
    assert 2 == mutations_df.loc['reference:190:A:?', 'Count']
    assert 'UNKNOWN_MISSING' == mutations_df.loc['reference:190:A:?', 'Type']
    assert 1 == mutations_df.loc['reference:210:C:?', 'Count']
    assert 'UNKNOWN_MISSING' == mutations_df.loc['reference:210:C:?', 'Type']
    assert 'reference:619:G:C' not in mutations_df

    # Test with different id type where deletion is int instead of sequence
    mutations_summarizer = MutationFeaturesFromIndexComparator(connection=loaded_database_genomic_data_store.connection,
                                                               include_unknown_samples=False,
                                                               ignore_annotations=True, id_type='spdi')
    mutations_df = mutations_summarizer.summary(present_set)
    mutations_df['Percent'] = mutations_df['Percent'].astype(int)  # Convert to int for easier comparison
    mutations_df['Unknown Percent'] = mutations_df['Unknown Percent'].fillna('<NA>')
    mutations_df['Unknown Count'] = mutations_df['Unknown Count'].fillna('<NA>')
    mutations_df['Present and Unknown Percent'] = mutations_df['Present and Unknown Percent'].fillna('<NA>')
    mutations_df['Present and Unknown Count'] = mutations_df['Present and Unknown Count'].fillna('<NA>')
    mutations_df = mutations_df.sort_index()

    assert 112 == len(mutations_df)
    assert list(expected_df.columns) == list(mutations_df.columns)
    assert 2 == mutations_df.loc['reference:619:1:C', 'Count']
    assert 'SNP' == mutations_df.loc['reference:619:1:C', 'Type']
    assert 2 == mutations_df.loc['reference:3063:1:ATGCAGC', 'Count']
    assert 'INDEL' == mutations_df.loc['reference:3063:1:ATGCAGC', 'Type']
    assert 1 == mutations_df.loc['reference:1984:7:TTGA', 'Count']
    assert 'OTHER' == mutations_df.loc['reference:1984:7:TTGA', 'Type']
    assert 1 == mutations_df.loc['reference:866:9:G', 'Count']
    assert 'INDEL' == mutations_df.loc['reference:866:9:G', 'Type']

    # Test with different id type include unknown features
    mutations_summarizer = MutationFeaturesFromIndexComparator(connection=loaded_database_genomic_data_store.connection,
                                                               ignore_annotations=True, id_type='spdi',
                                                               include_unknown_samples=False,
                                                               include_unknown=True)
    mutations_df = mutations_summarizer.summary(present_set)
    mutations_df['Percent'] = mutations_df['Percent'].astype(int)  # Convert to int for easier comparison
    mutations_df['Unknown Percent'] = mutations_df['Unknown Percent'].fillna('<NA>')
    mutations_df['Unknown Count'] = mutations_df['Unknown Count'].fillna('<NA>')
    mutations_df['Present and Unknown Percent'] = mutations_df['Present and Unknown Percent'].fillna('<NA>')
    mutations_df['Present and Unknown Count'] = mutations_df['Present and Unknown Count'].fillna('<NA>')
    mutations_df = mutations_df.sort_index()

    assert 112 + 440 == len(mutations_df)
    assert list(expected_df.columns) == list(mutations_df.columns)
    assert 2 == mutations_df.loc['reference:619:1:C', 'Count']
    assert 'SNP' == mutations_df.loc['reference:619:1:C', 'Type']
    assert 2 == mutations_df.loc['reference:3063:1:ATGCAGC', 'Count']
    assert 'INDEL' == mutations_df.loc['reference:3063:1:ATGCAGC', 'Type']
    assert 1 == mutations_df.loc['reference:1984:7:TTGA', 'Count']
    assert 'OTHER' == mutations_df.loc['reference:1984:7:TTGA', 'Type']
    assert 1 == mutations_df.loc['reference:866:9:G', 'Count']
    assert 'INDEL' == mutations_df.loc['reference:866:9:G', 'Type']
    assert 3 == mutations_df.loc['reference:90:1:?', 'Count']
    assert 'UNKNOWN_MISSING' == mutations_df.loc['reference:90:1:?', 'Type']
    assert 2 == mutations_df.loc['reference:190:1:?', 'Count']
    assert 'UNKNOWN_MISSING' == mutations_df.loc['reference:190:1:?', 'Type']
    assert 1 == mutations_df.loc['reference:210:1:?', 'Count']
    assert 'UNKNOWN_MISSING' == mutations_df.loc['reference:210:1:?', 'Type']


def test_summary_unique(loaded_database_genomic_data_store: GenomicsDataIndex):
    db = loaded_database_genomic_data_store.connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()
    all_sample_ids = {s.id for s in db.get_session().query(Sample).all()}

    mutations_summarizer = MutationFeaturesFromIndexComparator(connection=loaded_database_genomic_data_store.connection,
                                                               include_unknown_samples=False,
                                                               ignore_annotations=True)

    # Unique to A
    present_set = SampleSet({sampleA.id})
    other_set = SampleSet(all_sample_ids - {sampleA.id})
    mutations_df = mutations_summarizer.unique_summary(present_set, other_set=other_set).sort_index()

    expected_df = read_expected_snippy_df(snippy_all_dataframes['SampleA'], total=1)

    mutations_df['Percent'] = mutations_df['Percent'].astype(int)  # Convert to int for easier comparison

    assert len(expected_df) == len(mutations_df)
    assert 46 == len(mutations_df)  # Check length against independently generated length
    assert list(expected_df.index) == list(mutations_df.index)
    assert list(expected_df['Count']) == list(mutations_df['Count'])
    assert list(expected_df['Total']) == list(mutations_df['Total'])
    assert 100 == mutations_df.loc['reference:3656:CATT:C', 'Percent']
    assert 100 == mutations_df.loc['reference:461:AAAT:G', 'Percent']

    # Unique to B
    present_set = SampleSet({sampleB.id})
    other_set = SampleSet(all_sample_ids - {sampleB.id})
    mutations_df = mutations_summarizer.unique_summary(present_set, other_set=other_set).sort_index()

    expected_df = read_expected_snippy_df(snippy_all_dataframes['SampleB'],
                                          mutations_not_in=[snippy_all_dataframes['SampleA'],
                                                            snippy_all_dataframes['SampleC']],
                                          total=1)

    mutations_df['Percent'] = mutations_df['Percent'].astype(int)  # Convert to int for easier comparison

    assert len(expected_df) == len(mutations_df)
    assert list(expected_df.index) == list(mutations_df.index)
    assert list(expected_df['Count']) == list(mutations_df['Count'])
    assert list(expected_df['Total']) == list(mutations_df['Total'])
    assert 100 == mutations_df.loc['reference:349:AAGT:A', 'Percent']

    # Unique to BC
    present_set = SampleSet({sampleB.id, sampleC.id})
    other_set = SampleSet(all_sample_ids - {sampleB.id, sampleC.id})
    mutations_df = mutations_summarizer.unique_summary(present_set, other_set=other_set).sort_index()

    expected_df = read_expected_snippy_df([snippy_all_dataframes['SampleB'],
                                           snippy_all_dataframes['SampleC']],
                                          mutations_not_in=[snippy_all_dataframes['SampleA']],
                                          total=2)

    mutations_df['Percent'] = mutations_df['Percent'].astype(int)  # Convert to int for easier comparison

    assert len(expected_df) == len(mutations_df)
    assert 66 == len(mutations_df)  # Check length against independently generated length
    assert list(expected_df.index) == list(mutations_df.index)
    assert list(expected_df['Count']) == list(mutations_df['Count'])
    assert list(expected_df['Total']) == list(mutations_df['Total'])
    assert 100 == mutations_df.loc['reference:619:G:C', 'Percent']
    assert 50 == mutations_df.loc['reference:866:GCCAGATCC:G', 'Percent']
    assert 50 == mutations_df.loc['reference:349:AAGT:A', 'Percent']


def test_summary_annotations(loaded_database_genomic_data_store_annotations: GenomicsDataIndex):
    db = loaded_database_genomic_data_store_annotations.connection.database

    mutations_summarizer = MutationFeaturesFromIndexComparator(
        connection=loaded_database_genomic_data_store_annotations.connection,
        include_unknown_samples=False,
        ignore_annotations=False)

    sample_sh14_001 = db.get_session().query(Sample).filter(Sample.name == 'SH14-001').one()
    sample_sh14_014 = db.get_session().query(Sample).filter(Sample.name == 'SH14-014').one()
    sample_sh10_014 = db.get_session().query(Sample).filter(Sample.name == 'SH10-014').one()
    three_samples = {sample_sh14_001.id, sample_sh14_014.id, sample_sh10_014.id}

    present_set = SampleSet(three_samples)
    mutations_df = mutations_summarizer.summary(present_set)

    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type',
            'Count', 'Unknown Count', 'Present and Unknown Count', 'Total',
            'Percent', 'Unknown Percent', 'Present and Unknown Percent',
            'Annotation', 'Annotation_Impact',
            'Gene_Name', 'Gene_ID', 'Feature_Type', 'Transcript_BioType',
            'HGVS.c', 'HGVS.p', 'ID_HGVS.c', 'ID_HGVS.p', 'ID_HGVS_GN.c', 'ID_HGVS_GN.p'] == list(mutations_df.columns)
    assert 177 == len(mutations_df)
    mutations_df['Percent'] = mutations_df['Percent'].astype(int)  # easier to compare percents in assert
    mutations_df['Unknown Count'] = mutations_df['Unknown Count'].fillna('<NA>')
    mutations_df['Unknown Percent'] = mutations_df['Unknown Percent'].fillna('<NA>')
    mutations_df['Present and Unknown Count'] = mutations_df['Present and Unknown Count'].fillna('<NA>')
    mutations_df['Present and Unknown Percent'] = mutations_df['Present and Unknown Percent'].fillna('<NA>')

    # missense variant (3/3)
    assert ['NC_011083', 140658, 'C', 'A', 'SNP', 3, '<NA>', '<NA>', 3, 100, '<NA>', '<NA>',
            'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
            'c.497C>A', 'p.Ala166Glu',
            'hgvs:NC_011083:SEHA_RS01180:c.497C>A', 'hgvs:NC_011083:SEHA_RS01180:p.Ala166Glu',
            'hgvs_gn:NC_011083:murF:c.497C>A', 'hgvs_gn:NC_011083:murF:p.Ala166Glu'] == list(
        mutations_df.loc['NC_011083:140658:C:A'])

    # Intergenic variant (1/3)
    assert ['NC_011083', 4555461, 'T', 'TC', 'INDEL', 1, '<NA>', '<NA>', 3, 33, '<NA>', '<NA>',
            'intergenic_region', 'MODIFIER', 'SEHA_RS22510-SEHA_RS26685', 'SEHA_RS22510-SEHA_RS26685',
            'intergenic_region', 'NA',
            'n.4555461_4555462insC', 'NA',
            'hgvs:NC_011083:n.4555461_4555462insC', 'NA',
            'hgvs_gn:NC_011083:n.4555461_4555462insC', 'NA'] == list(
        mutations_df.loc['NC_011083:4555461:T:TC'].fillna('NA'))

    # MNP variant (1/3)
    assert ['NC_011083', 3535698, 'GCC', 'CAT', 'MNP', 2, '<NA>', '<NA>', 3, 66, '<NA>', '<NA>',
            'missense_variant', 'MODERATE', 'oadA', 'SEHA_RS17780',
            'transcript', 'protein_coding',
            'c.544_546delGGCinsATG', 'p.Gly182Met',
            'hgvs:NC_011083:SEHA_RS17780:c.544_546delGGCinsATG', 'hgvs:NC_011083:SEHA_RS17780:p.Gly182Met',
            'hgvs_gn:NC_011083:oadA:c.544_546delGGCinsATG', 'hgvs_gn:NC_011083:oadA:p.Gly182Met'] == list(
        mutations_df.loc['NC_011083:3535698:GCC:CAT'])


def test_summary_annotations_unknown(loaded_database_genomic_data_store_annotations_include_unknown: GenomicsDataIndex):
    db = loaded_database_genomic_data_store_annotations_include_unknown.connection.database

    mutations_summarizer = MutationFeaturesFromIndexComparator(
        connection=loaded_database_genomic_data_store_annotations_include_unknown.connection,
        include_unknown_samples=True,
        ignore_annotations=False)

    sample_sh14_001 = db.get_session().query(Sample).filter(Sample.name == 'SH14-001').one()
    sample_sh14_014 = db.get_session().query(Sample).filter(Sample.name == 'SH14-014').one()
    sample_sh10_014 = db.get_session().query(Sample).filter(Sample.name == 'SH10-014').one()
    three_samples = {sample_sh14_001.id, sample_sh14_014.id, sample_sh10_014.id}

    present_set = SampleSet(three_samples)
    mutations_df = mutations_summarizer.summary(present_set)

    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type',
            'Count', 'Unknown Count', 'Present and Unknown Count', 'Total',
            'Percent', 'Unknown Percent', 'Present and Unknown Percent',
            'Annotation', 'Annotation_Impact',
            'Gene_Name', 'Gene_ID', 'Feature_Type', 'Transcript_BioType',
            'HGVS.c', 'HGVS.p', 'ID_HGVS.c', 'ID_HGVS.p', 'ID_HGVS_GN.c', 'ID_HGVS_GN.p'] == list(mutations_df.columns)
    assert 177 == len(mutations_df)
    mutations_df['Percent'] = mutations_df['Percent'].astype(int)  # easier to compare percents in assert
    mutations_df = mutations_df.fillna('<NA>')
    mutations_df['Unknown Count'] = mutations_df['Unknown Count'].astype(int)
    mutations_df['Unknown Percent'] = mutations_df['Unknown Percent'].astype(int)
    mutations_df['Present and Unknown Count'] = mutations_df['Present and Unknown Count'].astype(int)
    mutations_df['Present and Unknown Percent'] = mutations_df['Present and Unknown Percent'].astype(int)

    # missense variant (3/3, 0/3)
    assert ['NC_011083', 140658, 'C', 'A', 'SNP', 3, 0, 3, 3, 100, 0, 100,
            'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
            'c.497C>A', 'p.Ala166Glu',
            'hgvs:NC_011083:SEHA_RS01180:c.497C>A', 'hgvs:NC_011083:SEHA_RS01180:p.Ala166Glu',
            'hgvs_gn:NC_011083:murF:c.497C>A', 'hgvs_gn:NC_011083:murF:p.Ala166Glu'] == list(
        mutations_df.loc['NC_011083:140658:C:A'])

    # Intergenic variant (1/3, 2/3)
    assert ['NC_011083', 4555461, 'T', 'TC', 'INDEL', 1, 2, 3, 3, 33, 66, 100,
            'intergenic_region', 'MODIFIER', 'SEHA_RS22510-SEHA_RS26685', 'SEHA_RS22510-SEHA_RS26685',
            'intergenic_region', '<NA>',
            'n.4555461_4555462insC', '<NA>',
            'hgvs:NC_011083:n.4555461_4555462insC', '<NA>',
            'hgvs_gn:NC_011083:n.4555461_4555462insC', '<NA>'] == list(
        mutations_df.loc['NC_011083:4555461:T:TC'])

    # MNP variant (2/3, 1/3)
    assert ['NC_011083', 3535698, 'GCC', 'CAT', 'MNP', 2, 1, 3, 3, 66, 33, 100,
            'missense_variant', 'MODERATE', 'oadA', 'SEHA_RS17780',
            'transcript', 'protein_coding',
            'c.544_546delGGCinsATG', 'p.Gly182Met',
            'hgvs:NC_011083:SEHA_RS17780:c.544_546delGGCinsATG', 'hgvs:NC_011083:SEHA_RS17780:p.Gly182Met',
            'hgvs_gn:NC_011083:oadA:c.544_546delGGCinsATG', 'hgvs_gn:NC_011083:oadA:p.Gly182Met'] == list(
        mutations_df.loc['NC_011083:3535698:GCC:CAT'])

    # Long MNP variant (2/3, 1/3)
    assert ['NC_011083', 3535143, 'AATGCCTGCC', 'TATCCCGGCG', 'MNP', 2, 1, 3, 3, 66, 33, 100,
            'synonymous_variant', 'LOW', 'oadA', 'SEHA_RS17780', 'transcript', 'protein_coding',
            'c.1092_1101delGGCAGGCATTinsCGCCGGGATA', 'p.368',
            'hgvs:NC_011083:SEHA_RS17780:c.1092_1101delGGCAGGCATTinsCGCCGGGATA',
            'hgvs:NC_011083:SEHA_RS17780:p.368',
            'hgvs_gn:NC_011083:oadA:c.1092_1101delGGCAGGCATTinsCGCCGGGATA',
            'hgvs_gn:NC_011083:oadA:p.368'] == list(
        mutations_df.loc['NC_011083:3535143:AATGCCTGCC:TATCCCGGCG'])

    # synonymous variant, no unknown (1/3, 0/3)
    assert ['NC_011083', 508378, 'C', 'T', 'SNP', 1, 0, 1, 3, 33, 0, 33,
            'synonymous_variant', 'LOW', 'tgt', 'SEHA_RS02965', 'transcript', 'protein_coding',
            'c.423C>T', 'p.Ile141Ile',
            'hgvs:NC_011083:SEHA_RS02965:c.423C>T', 'hgvs:NC_011083:SEHA_RS02965:p.Ile141Ile',
            'hgvs_gn:NC_011083:tgt:c.423C>T', 'hgvs_gn:NC_011083:tgt:p.Ile141Ile'] == list(
        mutations_df.loc['NC_011083:508378:C:T'])

    # variant where there is an overlap with present and unknown, (present 2/3, unknown 1/3, overlap 3869320)
    assert ['NC_011083', 3869320, 'C', 'A', 'SNP', 1, 1, 2, 3, 33, 33, 66,
            'synonymous_variant', 'LOW', 'yiaK', 'SEHA_RS19360', 'transcript', 'protein_coding',
            'c.591C>A', 'p.Gly197Gly',
            'hgvs:NC_011083:SEHA_RS19360:c.591C>A', 'hgvs:NC_011083:SEHA_RS19360:p.Gly197Gly',
            'hgvs_gn:NC_011083:yiaK:c.591C>A', 'hgvs_gn:NC_011083:yiaK:p.Gly197Gly'] == list(
        mutations_df.loc['NC_011083:3869320:C:A'])

    # deletion where there is an overlap with present and unknown, (present 3/3, unknown 2/3,
    # overlap 1676762 and 1676763)
    assert ['NC_011083', 1676762, 'CA', 'C', 'INDEL', 1, 2, 3, 3, 33, 66, 100,
            'intergenic_region', 'MODIFIER', 'SEHA_RS26130-SEHA_RS08880', 'SEHA_RS26130-SEHA_RS08880',
            'intergenic_region', '<NA>',
            'n.1676763delA', '<NA>',
            'hgvs:NC_011083:n.1676763delA', '<NA>',
            'hgvs_gn:NC_011083:n.1676763delA', '<NA>'] == list(
        mutations_df.loc['NC_011083:1676762:CA:C'])

    # All unknown (should not exist in table)
    assert 'NC_011083:1:A:C' not in mutations_df


def test_summary_no_annotations_unknown(
        loaded_database_genomic_data_store_annotations_include_unknown: GenomicsDataIndex):
    db = loaded_database_genomic_data_store_annotations_include_unknown.connection.database

    mutations_summarizer = MutationFeaturesFromIndexComparator(
        connection=loaded_database_genomic_data_store_annotations_include_unknown.connection,
        include_unknown_samples=True,
        ignore_annotations=True)

    sample_sh14_001 = db.get_session().query(Sample).filter(Sample.name == 'SH14-001').one()
    sample_sh14_014 = db.get_session().query(Sample).filter(Sample.name == 'SH14-014').one()
    sample_sh10_014 = db.get_session().query(Sample).filter(Sample.name == 'SH10-014').one()
    three_samples = {sample_sh14_001.id, sample_sh14_014.id, sample_sh10_014.id}

    present_set = SampleSet(three_samples)
    mutations_df = mutations_summarizer.summary(present_set)

    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type',
            'Count', 'Unknown Count', 'Present and Unknown Count', 'Total',
            'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(mutations_df.columns)
    assert 177 == len(mutations_df)
    mutations_df['Percent'] = mutations_df['Percent'].astype(int)  # easier to compare percents in assert
    mutations_df = mutations_df.fillna('<NA>')
    mutations_df['Unknown Count'] = mutations_df['Unknown Count'].astype(int)
    mutations_df['Unknown Percent'] = mutations_df['Unknown Percent'].astype(int)
    mutations_df['Present and Unknown Count'] = mutations_df['Present and Unknown Count'].astype(int)
    mutations_df['Present and Unknown Percent'] = mutations_df['Present and Unknown Percent'].astype(int)

    # missense variant (3/3, 0/3)
    assert ['NC_011083', 140658, 'C', 'A', 'SNP', 3, 0, 3, 3, 100, 0, 100] == list(
        mutations_df.loc['NC_011083:140658:C:A'])

    # Intergenic variant (1/3, 2/3)
    assert ['NC_011083', 4555461, 'T', 'TC', 'INDEL', 1, 2, 3, 3, 33, 66, 100] == list(
        mutations_df.loc['NC_011083:4555461:T:TC'])


def test_summary_annotations_unknown_column_unknown_rows(
        loaded_database_genomic_data_store_annotations_include_unknown: GenomicsDataIndex):
    db = loaded_database_genomic_data_store_annotations_include_unknown.connection.database

    mutations_summarizer = MutationFeaturesFromIndexComparator(
        connection=loaded_database_genomic_data_store_annotations_include_unknown.connection,
        include_unknown_samples=True,
        include_unknown=True,
        ignore_annotations=False)

    sample_sh14_001 = db.get_session().query(Sample).filter(Sample.name == 'SH14-001').one()
    sample_sh14_014 = db.get_session().query(Sample).filter(Sample.name == 'SH14-014').one()
    sample_sh10_014 = db.get_session().query(Sample).filter(Sample.name == 'SH10-014').one()
    three_samples = {sample_sh14_001.id, sample_sh14_014.id, sample_sh10_014.id}

    present_set = SampleSet(three_samples)
    mutations_df = mutations_summarizer.summary(present_set)

    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type',
            'Count', 'Unknown Count', 'Present and Unknown Count', 'Total',
            'Percent', 'Unknown Percent', 'Present and Unknown Percent',
            'Annotation', 'Annotation_Impact',
            'Gene_Name', 'Gene_ID', 'Feature_Type', 'Transcript_BioType',
            'HGVS.c', 'HGVS.p', 'ID_HGVS.c', 'ID_HGVS.p', 'ID_HGVS_GN.c', 'ID_HGVS_GN.p'] == list(mutations_df.columns)
    assert 648 == len(mutations_df)
    mutations_df['Percent'] = mutations_df['Percent'].astype(int)  # easier to compare percents in assert
    mutations_df['Unknown Count'] = mutations_df['Unknown Count'].fillna(-1)
    mutations_df['Unknown Percent'] = mutations_df['Unknown Percent'].fillna(-1)
    mutations_df['Present and Unknown Count'] = mutations_df['Present and Unknown Count'].fillna(-1)
    mutations_df['Present and Unknown Percent'] = mutations_df['Present and Unknown Percent'].fillna(-1)
    mutations_df = mutations_df.fillna('<NA>')
    mutations_df['Unknown Count'] = mutations_df['Unknown Count'].astype(int)
    mutations_df['Unknown Percent'] = mutations_df['Unknown Percent'].astype(int)
    mutations_df['Present and Unknown Count'] = mutations_df['Present and Unknown Count'].astype(int)
    mutations_df['Present and Unknown Percent'] = mutations_df['Present and Unknown Percent'].astype(int)

    # missense variant (3/3, 0/3)
    assert ['NC_011083', 140658, 'C', 'A', 'SNP', 3, 0, 3, 3, 100, 0, 100,
            'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
            'c.497C>A', 'p.Ala166Glu',
            'hgvs:NC_011083:SEHA_RS01180:c.497C>A', 'hgvs:NC_011083:SEHA_RS01180:p.Ala166Glu',
            'hgvs_gn:NC_011083:murF:c.497C>A', 'hgvs_gn:NC_011083:murF:p.Ala166Glu'] == list(
        mutations_df.loc['NC_011083:140658:C:A'])

    # Intergenic variant (1/3, 2/3)
    assert ['NC_011083', 4555461, 'T', 'TC', 'INDEL', 1, 2, 3, 3, 33, 66, 100,
            'intergenic_region', 'MODIFIER', 'SEHA_RS22510-SEHA_RS26685', 'SEHA_RS22510-SEHA_RS26685',
            'intergenic_region', '<NA>',
            'n.4555461_4555462insC', '<NA>',
            'hgvs:NC_011083:n.4555461_4555462insC', '<NA>',
            'hgvs_gn:NC_011083:n.4555461_4555462insC', '<NA>'] == list(
        mutations_df.loc['NC_011083:4555461:T:TC'])

    # Long deletion
    assert ['NC_011083', 3167187, 'AACCACGACCACGACCACGACCACGACCACGACCACG', 'A', 'INDEL', 2, 0, 2, 3, 66, 0, 66,
            'disruptive_inframe_deletion', 'MODERATE', 'SEHA_RS15905', 'SEHA_RS15905',
            'transcript', 'protein_coding',
            'c.429_464delCGACCACGACCACGACCACGACCACGACCACGACCA', 'p.Asp144_His155del',
            'hgvs:NC_011083:SEHA_RS15905:c.429_464delCGACCACGACCACGACCACGACCACGACCACGACCA',
            'hgvs:NC_011083:SEHA_RS15905:p.Asp144_His155del',
            'hgvs_gn:NC_011083:SEHA_RS15905:c.429_464delCGACCACGACCACGACCACGACCACGACCACGACCA',
            'hgvs_gn:NC_011083:SEHA_RS15905:p.Asp144_His155del'] == list(
        mutations_df.loc['NC_011083:3167187:AACCACGACCACGACCACGACCACGACCACGACCACG:A'])

    # variant where there is an overlap with present and unknown, (present 2/3, unknown 1/3, overlap 3869320)
    assert ['NC_011083', 3869320, 'C', 'A', 'SNP', 1, 1, 2, 3, 33, 33, 66,
            'synonymous_variant', 'LOW', 'yiaK', 'SEHA_RS19360', 'transcript', 'protein_coding',
            'c.591C>A', 'p.Gly197Gly',
            'hgvs:NC_011083:SEHA_RS19360:c.591C>A', 'hgvs:NC_011083:SEHA_RS19360:p.Gly197Gly',
            'hgvs_gn:NC_011083:yiaK:c.591C>A', 'hgvs_gn:NC_011083:yiaK:p.Gly197Gly'] == list(
        mutations_df.loc['NC_011083:3869320:C:A'])

    # deletion where there is an overlap with present and unknown, (present 3/3, unknown 2/3,
    # overlap 1676762 and 1676763)
    assert ['NC_011083', 1676762, 'CA', 'C', 'INDEL', 1, 2, 3, 3, 33, 66, 100,
            'intergenic_region', 'MODIFIER', 'SEHA_RS26130-SEHA_RS08880', 'SEHA_RS26130-SEHA_RS08880',
            'intergenic_region', '<NA>',
            'n.1676763delA', '<NA>',
            'hgvs:NC_011083:n.1676763delA', '<NA>',
            'hgvs_gn:NC_011083:n.1676763delA', '<NA>'] == list(
        mutations_df.loc['NC_011083:1676762:CA:C'])

    # The unknown feature for the above case
    assert ['NC_011083', 3869320, 'C', '?', 'UNKNOWN_MISSING', 1, -1, -1, 3, 33, -1, -1] + 12 * ['<NA>'] == list(
        mutations_df.loc['NC_011083:3869320:C:?'])

    # Another Unknown position
    assert ['NC_011083', 145096, 'A', '?', 'UNKNOWN_MISSING', 1, -1, -1, 3, 33, -1, -1] + 12 * ['<NA>'] == list(
        mutations_df.loc['NC_011083:145096:A:?'])


def test_features_comparison(loaded_database_genomic_data_store: GenomicsDataIndex):
    db = loaded_database_genomic_data_store.connection.database
    sampleA = db.get_session().query(Sample).filter(Sample.name == 'SampleA').one()
    sampleB = db.get_session().query(Sample).filter(Sample.name == 'SampleB').one()
    sampleC = db.get_session().query(Sample).filter(Sample.name == 'SampleC').one()
    all_sample_ids = {s.id for s in db.get_session().query(Sample).all()}

    present_set = SampleSet(all_sample_ids)
    mutations_summarizer = MutationFeaturesFromIndexComparator(
        connection=loaded_database_genomic_data_store.connection,
        ignore_annotations=True)

    # Test single category of all
    sample_categories = [present_set]
    comparison_df = mutations_summarizer.features_comparison(selected_samples=present_set,
                                                             sample_categories=sample_categories,
                                                             category_prefixes=['All'],
                                                             unit='count')
    comparison_df = comparison_df.sort_index()
    assert comparison_df.index.name == 'Mutation'
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type', 'Total',
            'All_count', 'All_Unknown count', 'All_Present and Unknown count',
            'All_total'] == comparison_df.columns.tolist()
    assert {9} == set(comparison_df['Total'].tolist())
    assert {9} == set(comparison_df['All_total'].tolist())
    assert 2 == comparison_df.loc['reference:619:G:C', 'All_count']
    assert 0 == comparison_df.loc['reference:619:G:C', 'All_Unknown count']
    assert 2 == comparison_df.loc['reference:619:G:C', 'All_Present and Unknown count']
    assert 'SNP' == comparison_df.loc['reference:619:G:C', 'Type']
    assert 1 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'All_count']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'All_Unknown count']
    assert 1 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'All_Present and Unknown count']
    assert 'INDEL' == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'Type']
    assert 2 == comparison_df.loc['reference:4693:C:CGA', 'All_count']
    assert 0 == comparison_df.loc['reference:4693:C:CGA', 'All_Unknown count']
    assert 2 == comparison_df.loc['reference:4693:C:CGA', 'All_Present and Unknown count']
    assert 'INDEL' == comparison_df.loc['reference:4693:C:CGA', 'Type']

    # Test two categories, one of A and one of BC
    sample_categories = [SampleSet([sampleA.id]), SampleSet([sampleB.id, sampleC.id])]
    comparison_df = mutations_summarizer.features_comparison(selected_samples=present_set,
                                                             sample_categories=sample_categories,
                                                             category_prefixes=['A', 'BC'],
                                                             unit='count')
    comparison_df = comparison_df.sort_index()
    assert comparison_df.index.name == 'Mutation'
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type', 'Total',
            'A_count', 'A_Unknown count', 'A_Present and Unknown count',
            'BC_count', 'BC_Unknown count', 'BC_Present and Unknown count',
            'A_total', 'BC_total'] == comparison_df.columns.tolist()
    assert {9} == set(comparison_df['Total'].tolist())
    assert {1} == set(comparison_df['A_total'].tolist())
    assert {2} == set(comparison_df['BC_total'].tolist())
    assert 0 == comparison_df.loc['reference:619:G:C', 'A_count']
    assert 0 == comparison_df.loc['reference:619:G:C', 'A_Unknown count']
    assert 0 == comparison_df.loc['reference:619:G:C', 'A_Present and Unknown count']
    assert 2 == comparison_df.loc['reference:619:G:C', 'BC_count']
    assert 0 == comparison_df.loc['reference:619:G:C', 'BC_Unknown count']
    assert 2 == comparison_df.loc['reference:619:G:C', 'BC_Present and Unknown count']
    assert 1 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'A_count']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'A_Unknown count']
    assert 1 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'A_Present and Unknown count']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'BC_count']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'BC_Unknown count']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'BC_Present and Unknown count']
    assert 0 == comparison_df.loc['reference:4693:C:CGA', 'A_count']
    assert 0 == comparison_df.loc['reference:4693:C:CGA', 'A_Unknown count']
    assert 0 == comparison_df.loc['reference:4693:C:CGA', 'A_Present and Unknown count']
    assert 2 == comparison_df.loc['reference:4693:C:CGA', 'BC_count']
    assert 0 == comparison_df.loc['reference:4693:C:CGA', 'BC_Unknown count']
    assert 2 == comparison_df.loc['reference:4693:C:CGA', 'BC_Present and Unknown count']

    # Test two categories, one of AB and one of C
    sample_categories = [SampleSet([sampleA.id, sampleB.id]), SampleSet([sampleC.id])]
    comparison_df = mutations_summarizer.features_comparison(selected_samples=present_set,
                                                             sample_categories=sample_categories,
                                                             category_prefixes=['AB', 'C'],
                                                             unit='count')
    comparison_df = comparison_df.sort_index()
    assert comparison_df.index.name == 'Mutation'
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type', 'Total',
            'AB_count', 'AB_Unknown count', 'AB_Present and Unknown count',
            'C_count', 'C_Unknown count', 'C_Present and Unknown count',
            'AB_total', 'C_total'] == comparison_df.columns.tolist()
    assert {9} == set(comparison_df['Total'].tolist())
    assert {2} == set(comparison_df['AB_total'].tolist())
    assert {1} == set(comparison_df['C_total'].tolist())
    assert 1 == comparison_df.loc['reference:619:G:C', 'AB_count']
    assert 0 == comparison_df.loc['reference:619:G:C', 'AB_Unknown count']
    assert 1 == comparison_df.loc['reference:619:G:C', 'AB_Present and Unknown count']
    assert 1 == comparison_df.loc['reference:619:G:C', 'C_count']
    assert 0 == comparison_df.loc['reference:619:G:C', 'C_Unknown count']
    assert 1 == comparison_df.loc['reference:619:G:C', 'C_Present and Unknown count']
    assert 1 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'AB_count']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'AB_Unknown count']
    assert 1 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'AB_Present and Unknown count']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'C_count']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'C_Unknown count']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'C_Present and Unknown count']
    assert 1 == comparison_df.loc['reference:4693:C:CGA', 'AB_count']
    assert 0 == comparison_df.loc['reference:4693:C:CGA', 'AB_Unknown count']
    assert 1 == comparison_df.loc['reference:4693:C:CGA', 'AB_Present and Unknown count']
    assert 1 == comparison_df.loc['reference:4693:C:CGA', 'C_count']
    assert 0 == comparison_df.loc['reference:4693:C:CGA', 'C_Unknown count']
    assert 1 == comparison_df.loc['reference:4693:C:CGA', 'C_Present and Unknown count']
    assert 1 == comparison_df.loc['reference:528:C:CAG', 'AB_count']
    assert 0 == comparison_df.loc['reference:528:C:CAG', 'AB_Unknown count']
    assert 1 == comparison_df.loc['reference:528:C:CAG', 'AB_Present and Unknown count']
    assert 0 == comparison_df.loc['reference:528:C:CAG', 'C_count']
    assert 0 == comparison_df.loc['reference:528:C:CAG', 'C_Unknown count']
    assert 0 == comparison_df.loc['reference:528:C:CAG', 'C_Present and Unknown count']

    # Test three categories: A, B, and C, and total out of only these 3
    sample_categories = [SampleSet([sampleA.id]), SampleSet([sampleB.id]), SampleSet([sampleC.id])]
    selected_samples = SampleSet([sampleA.id, sampleB.id, sampleC.id])
    comparison_df = mutations_summarizer.features_comparison(selected_samples=selected_samples,
                                                             sample_categories=sample_categories,
                                                             category_prefixes=['A', 'B', 'C'],
                                                             unit='count')
    comparison_df = comparison_df.sort_index()
    assert comparison_df.index.name == 'Mutation'
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type', 'Total',
            'A_count', 'A_Unknown count', 'A_Present and Unknown count',
            'B_count', 'B_Unknown count', 'B_Present and Unknown count',
            'C_count', 'C_Unknown count', 'C_Present and Unknown count',
            'A_total', 'B_total', 'C_total'] == comparison_df.columns.tolist()
    assert {3} == set(comparison_df['Total'].tolist())
    assert {1} == set(comparison_df['A_total'].tolist())
    assert {1} == set(comparison_df['B_total'].tolist())
    assert {1} == set(comparison_df['C_total'].tolist())
    assert 0 == comparison_df.loc['reference:619:G:C', 'A_count']
    assert 0 == comparison_df.loc['reference:619:G:C', 'A_Unknown count']
    assert 0 == comparison_df.loc['reference:619:G:C', 'A_Present and Unknown count']
    assert 1 == comparison_df.loc['reference:619:G:C', 'B_count']
    assert 0 == comparison_df.loc['reference:619:G:C', 'B_Unknown count']
    assert 1 == comparison_df.loc['reference:619:G:C', 'B_Present and Unknown count']
    assert 1 == comparison_df.loc['reference:619:G:C', 'C_count']
    assert 0 == comparison_df.loc['reference:619:G:C', 'C_Unknown count']
    assert 1 == comparison_df.loc['reference:619:G:C', 'C_Present and Unknown count']
    assert 1 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'A_count']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'A_Unknown count']
    assert 1 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'A_Present and Unknown count']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'B_count']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'B_Unknown count']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'B_Present and Unknown count']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'C_count']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'C_Unknown count']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'C_Present and Unknown count']
    assert 0 == comparison_df.loc['reference:4693:C:CGA', 'A_count']
    assert 0 == comparison_df.loc['reference:4693:C:CGA', 'A_Unknown count']
    assert 0 == comparison_df.loc['reference:4693:C:CGA', 'A_Present and Unknown count']
    assert 1 == comparison_df.loc['reference:4693:C:CGA', 'B_count']
    assert 0 == comparison_df.loc['reference:4693:C:CGA', 'B_Unknown count']
    assert 1 == comparison_df.loc['reference:4693:C:CGA', 'B_Present and Unknown count']
    assert 1 == comparison_df.loc['reference:4693:C:CGA', 'C_count']
    assert 0 == comparison_df.loc['reference:4693:C:CGA', 'C_Unknown count']
    assert 1 == comparison_df.loc['reference:4693:C:CGA', 'C_Present and Unknown count']
    assert 0 == comparison_df.loc['reference:528:C:CAG', 'A_count']
    assert 0 == comparison_df.loc['reference:528:C:CAG', 'A_Unknown count']
    assert 0 == comparison_df.loc['reference:528:C:CAG', 'A_Present and Unknown count']
    assert 1 == comparison_df.loc['reference:528:C:CAG', 'B_count']
    assert 0 == comparison_df.loc['reference:528:C:CAG', 'B_Unknown count']
    assert 1 == comparison_df.loc['reference:528:C:CAG', 'B_Present and Unknown count']
    assert 0 == comparison_df.loc['reference:528:C:CAG', 'C_count']
    assert 0 == comparison_df.loc['reference:528:C:CAG', 'C_Unknown count']
    assert 0 == comparison_df.loc['reference:528:C:CAG', 'C_Present and Unknown count']

    # Test two categories: A, and BC, and percent
    sample_categories = [SampleSet([sampleA.id]), SampleSet([sampleB.id, sampleC.id])]
    comparison_df = mutations_summarizer.features_comparison(selected_samples=present_set,
                                                             sample_categories=sample_categories,
                                                             category_prefixes=['A', 'BC'],
                                                             unit='percent')
    comparison_df = comparison_df.sort_index()
    comparison_df['A_percent'] = comparison_df['A_percent'].astype(int)  # Convert to int for easier comparison
    comparison_df['BC_percent'] = comparison_df['BC_percent'].astype(int)  # Convert to int for easier comparison
    assert comparison_df.index.name == 'Mutation'
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type', 'Total',
            'A_percent', 'A_Unknown percent', 'A_Present and Unknown percent',
            'BC_percent', 'BC_Unknown percent', 'BC_Present and Unknown percent',
            'A_total', 'BC_total'] == comparison_df.columns.tolist()
    assert {9} == set(comparison_df['Total'].tolist())
    assert {1} == set(comparison_df['A_total'].tolist())
    assert {2} == set(comparison_df['BC_total'].tolist())
    assert 0 == comparison_df.loc['reference:619:G:C', 'A_percent']
    assert 0 == comparison_df.loc['reference:619:G:C', 'A_Unknown percent']
    assert 0 == comparison_df.loc['reference:619:G:C', 'A_Present and Unknown percent']
    assert 100 == comparison_df.loc['reference:619:G:C', 'BC_percent']
    assert 0 == comparison_df.loc['reference:619:G:C', 'BC_Unknown percent']
    assert 100 == comparison_df.loc['reference:619:G:C', 'BC_Present and Unknown percent']
    assert 100 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'A_percent']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'A_Unknown percent']
    assert 100 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'A_Present and Unknown percent']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'BC_percent']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'BC_Unknown percent']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'BC_Present and Unknown percent']
    assert 0 == comparison_df.loc['reference:4693:C:CGA', 'A_percent']
    assert 0 == comparison_df.loc['reference:4693:C:CGA', 'A_Unknown percent']
    assert 0 == comparison_df.loc['reference:4693:C:CGA', 'A_Present and Unknown percent']
    assert 100 == comparison_df.loc['reference:4693:C:CGA', 'BC_percent']
    assert 0 == comparison_df.loc['reference:4693:C:CGA', 'BC_Unknown percent']
    assert 100 == comparison_df.loc['reference:4693:C:CGA', 'BC_Present and Unknown percent']

    # Test default category_names and compare_kind
    sample_categories = [SampleSet([sampleA.id]), SampleSet([sampleB.id, sampleC.id])]
    comparison_df = mutations_summarizer.features_comparison(selected_samples=present_set,
                                                             sample_categories=sample_categories)
    comparison_df = comparison_df.sort_index()
    comparison_df['Category1_percent'] = comparison_df['Category1_percent'].astype(
        int)  # Convert to int for easier comparison
    comparison_df['Category2_percent'] = comparison_df['Category2_percent'].astype(
        int)  # Convert to int for easier comparison
    assert comparison_df.index.name == 'Mutation'
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type', 'Total',
            'Category1_percent', 'Category1_Unknown percent', 'Category1_Present and Unknown percent',
            'Category2_percent', 'Category2_Unknown percent', 'Category2_Present and Unknown percent',
            'Category1_total', 'Category2_total'] == comparison_df.columns.tolist()
    assert {9} == set(comparison_df['Total'].tolist())
    assert {1} == set(comparison_df['Category1_total'].tolist())
    assert {2} == set(comparison_df['Category2_total'].tolist())
    assert 0 == comparison_df.loc['reference:619:G:C', 'Category1_percent']
    assert 0 == comparison_df.loc['reference:619:G:C', 'Category1_Unknown percent']
    assert 0 == comparison_df.loc['reference:619:G:C', 'Category1_Present and Unknown percent']
    assert 100 == comparison_df.loc['reference:619:G:C', 'Category2_percent']
    assert 0 == comparison_df.loc['reference:619:G:C', 'Category2_Unknown percent']
    assert 100 == comparison_df.loc['reference:619:G:C', 'Category2_Present and Unknown percent']
    assert 100 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'Category1_percent']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'Category1_Unknown percent']
    assert 100 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'Category1_Present and Unknown percent']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'Category2_percent']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'Category2_Unknown percent']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'Category2_Present and Unknown percent']
    assert 0 == comparison_df.loc['reference:4693:C:CGA', 'Category1_percent']
    assert 0 == comparison_df.loc['reference:4693:C:CGA', 'Category1_Unknown percent']
    assert 0 == comparison_df.loc['reference:4693:C:CGA', 'Category1_Present and Unknown percent']
    assert 100 == comparison_df.loc['reference:4693:C:CGA', 'Category2_percent']
    assert 0 == comparison_df.loc['reference:4693:C:CGA', 'Category2_Unknown percent']
    assert 100 == comparison_df.loc['reference:4693:C:CGA', 'Category2_Present and Unknown percent']


def test_features_comparison_annotations(loaded_database_genomic_data_store_annotations_include_unknown: GenomicsDataIndex):
    db = loaded_database_genomic_data_store_annotations_include_unknown.connection.database

    sample_sh14_001 = db.get_session().query(Sample).filter(Sample.name == 'SH14-001').one()
    sample_sh14_014 = db.get_session().query(Sample).filter(Sample.name == 'SH14-014').one()
    sample_sh10_014 = db.get_session().query(Sample).filter(Sample.name == 'SH10-014').one()
    three_samples = {sample_sh14_001.id, sample_sh14_014.id, sample_sh10_014.id}

    present_set = SampleSet(three_samples)
    mutations_summarizer = MutationFeaturesFromIndexComparator(
        connection=loaded_database_genomic_data_store_annotations_include_unknown.connection,
        ignore_annotations=False)

    # Test single category of all
    sample_categories = [present_set]
    comparison_df = mutations_summarizer.features_comparison(selected_samples=present_set,
                                                             sample_categories=sample_categories,
                                                             category_prefixes=['All'],
                                                             unit='count')
    comparison_df = comparison_df.sort_index()
    assert comparison_df.index.name == 'Mutation'
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type', 'Total',
            'All_count', 'All_Unknown count', 'All_Present and Unknown count',
            'All_total',
            'Annotation', 'Annotation_Impact',
            'Gene_Name', 'Gene_ID', 'Feature_Type', 'Transcript_BioType',
            'HGVS.c', 'HGVS.p', 'ID_HGVS.c', 'ID_HGVS.p', 'ID_HGVS_GN.c',
            'ID_HGVS_GN.p'] == list(comparison_df.columns)
    assert 177 == len(comparison_df)
    assert {3} == set(comparison_df['Total'].tolist())
    assert {3} == set(comparison_df['All_total'].tolist())
    assert 3 == comparison_df.loc['NC_011083:140658:C:A', 'All_count']
    assert 0 == comparison_df.loc['NC_011083:140658:C:A', 'All_Unknown count']
    assert 3 == comparison_df.loc['NC_011083:140658:C:A', 'All_Present and Unknown count']
    assert 'hgvs_gn:NC_011083:murF:p.Ala166Glu' == comparison_df.loc[
        'NC_011083:140658:C:A', 'ID_HGVS_GN.p']
    assert 1 == comparison_df.loc['NC_011083:4555461:T:TC', 'All_count']
    assert 2 == comparison_df.loc['NC_011083:4555461:T:TC', 'All_Unknown count']
    assert 3 == comparison_df.loc['NC_011083:4555461:T:TC', 'All_Present and Unknown count']
    assert 'hgvs_gn:NC_011083:n.4555461_4555462insC' == comparison_df.loc[
        'NC_011083:4555461:T:TC', 'ID_HGVS_GN.c']

    # Test 2 categories: one of SH10-014 and one of SH14-001, SH14-014
    sample_categories = [SampleSet([sample_sh10_014.id]), SampleSet([sample_sh14_001.id, sample_sh14_014.id])]
    comparison_df = mutations_summarizer.features_comparison(selected_samples=present_set,
                                                             sample_categories=sample_categories,
                                                             category_prefixes=['10', '14'],
                                                             unit='count')
    comparison_df = comparison_df.sort_index()
    comparison_df = comparison_df.fillna('<NA>')
    assert comparison_df.index.name == 'Mutation'
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type', 'Total',
            '10_count', '10_Unknown count', '10_Present and Unknown count',
            '14_count', '14_Unknown count', '14_Present and Unknown count',
            '10_total', '14_total',
            'Annotation', 'Annotation_Impact',
            'Gene_Name', 'Gene_ID', 'Feature_Type', 'Transcript_BioType',
            'HGVS.c', 'HGVS.p', 'ID_HGVS.c', 'ID_HGVS.p', 'ID_HGVS_GN.c',
            'ID_HGVS_GN.p'] == list(comparison_df.columns)
    assert 177 == len(comparison_df)
    assert {3} == set(comparison_df['Total'].tolist())
    assert {1} == set(comparison_df['10_total'].tolist())
    assert {2} == set(comparison_df['14_total'].tolist())
    assert 1 == comparison_df.loc['NC_011083:140658:C:A', '10_count']
    assert 0 == comparison_df.loc['NC_011083:140658:C:A', '10_Unknown count']
    assert 1 == comparison_df.loc['NC_011083:140658:C:A', '10_Present and Unknown count']
    assert 2 == comparison_df.loc['NC_011083:140658:C:A', '14_count']
    assert 0 == comparison_df.loc['NC_011083:140658:C:A', '14_Unknown count']
    assert 2 == comparison_df.loc['NC_011083:140658:C:A', '14_Present and Unknown count']
    assert 'hgvs_gn:NC_011083:murF:p.Ala166Glu' == comparison_df.loc[
        'NC_011083:140658:C:A', 'ID_HGVS_GN.p']
    assert 1 == comparison_df.loc['NC_011083:4555461:T:TC', '10_count']
    assert 0 == comparison_df.loc['NC_011083:4555461:T:TC', '10_Unknown count']
    assert 1 == comparison_df.loc['NC_011083:4555461:T:TC', '10_Present and Unknown count']
    assert 0 == comparison_df.loc['NC_011083:4555461:T:TC', '14_count']
    assert 2 == comparison_df.loc['NC_011083:4555461:T:TC', '14_Unknown count']
    assert 2 == comparison_df.loc['NC_011083:4555461:T:TC', '14_Present and Unknown count']
    assert 'hgvs_gn:NC_011083:n.4555461_4555462insC' == comparison_df.loc[
        'NC_011083:4555461:T:TC', 'ID_HGVS_GN.c']
    assert 0 == comparison_df.loc['NC_011083:4482211:C:A', '10_count']
    assert 0 == comparison_df.loc['NC_011083:4482211:C:A', '10_Unknown count']
    assert 0 == comparison_df.loc['NC_011083:4482211:C:A', '10_Present and Unknown count']
    assert 1 == comparison_df.loc['NC_011083:4482211:C:A', '14_count']
    assert 0 == comparison_df.loc['NC_011083:4482211:C:A', '14_Unknown count']
    assert 1 == comparison_df.loc['NC_011083:4482211:C:A', '14_Present and Unknown count']
    assert 'hgvs_gn:NC_011083:siiE:p.Arg1263Ser' == comparison_df.loc[
        'NC_011083:4482211:C:A', 'ID_HGVS_GN.p']
    assert 0 == comparison_df.loc['NC_011083:630556:G:A', '10_count']
    assert 0 == comparison_df.loc['NC_011083:630556:G:A', '10_Unknown count']
    assert 0 == comparison_df.loc['NC_011083:630556:G:A', '10_Present and Unknown count']
    assert 2 == comparison_df.loc['NC_011083:630556:G:A', '14_count']
    assert 0 == comparison_df.loc['NC_011083:630556:G:A', '14_Unknown count']
    assert 2 == comparison_df.loc['NC_011083:630556:G:A', '14_Present and Unknown count']
    assert 'hgvs_gn:NC_011083:SEHA_RS03545:p.Trp295*' == comparison_df.loc[
        'NC_011083:630556:G:A', 'ID_HGVS_GN.p']
    assert 0 == comparison_df.loc['NC_011083:3869320:C:A', '10_count']
    assert 0 == comparison_df.loc['NC_011083:3869320:C:A', '10_Unknown count']
    assert 0 == comparison_df.loc['NC_011083:3869320:C:A', '10_Present and Unknown count']
    assert 1 == comparison_df.loc['NC_011083:3869320:C:A', '14_count']
    assert 1 == comparison_df.loc['NC_011083:3869320:C:A', '14_Unknown count']
    assert 2 == comparison_df.loc['NC_011083:3869320:C:A', '14_Present and Unknown count']
    assert 'hgvs_gn:NC_011083:yiaK:p.Gly197Gly' == comparison_df.loc[
        'NC_011083:3869320:C:A', 'ID_HGVS_GN.p']
    assert 1 == comparison_df.loc['NC_011083:3535698:GCC:CAT', '10_count']
    assert 0 == comparison_df.loc['NC_011083:3535698:GCC:CAT', '10_Unknown count']
    assert 1 == comparison_df.loc['NC_011083:3535698:GCC:CAT', '10_Present and Unknown count']
    assert 1 == comparison_df.loc['NC_011083:3535698:GCC:CAT', '14_count']
    assert 1 == comparison_df.loc['NC_011083:3535698:GCC:CAT', '14_Unknown count']
    assert 2 == comparison_df.loc['NC_011083:3535698:GCC:CAT', '14_Present and Unknown count']
    assert 'hgvs_gn:NC_011083:oadA:p.Gly182Met' == comparison_df.loc[
        'NC_011083:3535698:GCC:CAT', 'ID_HGVS_GN.p']
    assert 0 == comparison_df.loc['NC_011083:1676762:CA:C', '10_count']
    assert 1 == comparison_df.loc['NC_011083:1676762:CA:C', '10_Unknown count']
    assert 1 == comparison_df.loc['NC_011083:1676762:CA:C', '10_Present and Unknown count']
    assert 1 == comparison_df.loc['NC_011083:1676762:CA:C', '14_count']
    assert 1 == comparison_df.loc['NC_011083:1676762:CA:C', '14_Unknown count']
    assert 2 == comparison_df.loc['NC_011083:1676762:CA:C', '14_Present and Unknown count']
    assert '<NA>' == comparison_df.loc[
        'NC_011083:1676762:CA:C', 'ID_HGVS_GN.p']
    # All unknown (should not exist in table)
    assert 'NC_011083:1:A:C' not in comparison_df

    # Test 2 categories defaults: one of SH10-014 and one of SH14-001, SH14-014
    sample_categories = [SampleSet([sample_sh10_014.id]), SampleSet([sample_sh14_001.id, sample_sh14_014.id])]
    comparison_df = mutations_summarizer.features_comparison(selected_samples=present_set,
                                                             sample_categories=sample_categories)
    comparison_df = comparison_df.sort_index()
    comparison_df = comparison_df.fillna('<NA>')
    comparison_df['Category1_percent'] = comparison_df['Category1_percent'].astype(
        int)  # Convert to int for easier comparison
    comparison_df['Category2_percent'] = comparison_df['Category2_percent'].astype(
        int)  # Convert to int for easier comparison
    assert comparison_df.index.name == 'Mutation'
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type', 'Total',
            'Category1_percent', 'Category1_Unknown percent', 'Category1_Present and Unknown percent',
            'Category2_percent', 'Category2_Unknown percent', 'Category2_Present and Unknown percent',
            'Category1_total', 'Category2_total',
            'Annotation', 'Annotation_Impact',
            'Gene_Name', 'Gene_ID', 'Feature_Type', 'Transcript_BioType',
            'HGVS.c', 'HGVS.p', 'ID_HGVS.c', 'ID_HGVS.p', 'ID_HGVS_GN.c',
            'ID_HGVS_GN.p'] == list(comparison_df.columns)
    assert 177 == len(comparison_df)
    assert {3} == set(comparison_df['Total'].tolist())
    assert {1} == set(comparison_df['Category1_total'].tolist())
    assert {2} == set(comparison_df['Category2_total'].tolist())
    assert 100 == comparison_df.loc['NC_011083:140658:C:A', 'Category1_percent']
    assert 0 == comparison_df.loc['NC_011083:140658:C:A', 'Category1_Unknown percent']
    assert 100 == comparison_df.loc['NC_011083:140658:C:A', 'Category1_Present and Unknown percent']
    assert 100 == comparison_df.loc['NC_011083:140658:C:A', 'Category2_percent']
    assert 0 == comparison_df.loc['NC_011083:140658:C:A', 'Category2_Unknown percent']
    assert 100 == comparison_df.loc['NC_011083:140658:C:A', 'Category2_Present and Unknown percent']
    assert 'hgvs_gn:NC_011083:murF:p.Ala166Glu' == comparison_df.loc[
        'NC_011083:140658:C:A', 'ID_HGVS_GN.p']
    assert 100 == comparison_df.loc['NC_011083:4555461:T:TC', 'Category1_percent']
    assert 0 == comparison_df.loc['NC_011083:4555461:T:TC', 'Category1_Unknown percent']
    assert 100 == comparison_df.loc['NC_011083:4555461:T:TC', 'Category1_Present and Unknown percent']
    assert 0 == comparison_df.loc['NC_011083:4555461:T:TC', 'Category2_percent']
    assert 100 == comparison_df.loc['NC_011083:4555461:T:TC', 'Category2_Unknown percent']
    assert 100 == comparison_df.loc['NC_011083:4555461:T:TC', 'Category2_Present and Unknown percent']
    assert 'hgvs_gn:NC_011083:n.4555461_4555462insC' == comparison_df.loc[
        'NC_011083:4555461:T:TC', 'ID_HGVS_GN.c']
    assert 0 == comparison_df.loc['NC_011083:4482211:C:A', 'Category1_percent']
    assert 0 == comparison_df.loc['NC_011083:4482211:C:A', 'Category1_Unknown percent']
    assert 0 == comparison_df.loc['NC_011083:4482211:C:A', 'Category1_Present and Unknown percent']
    assert 50 == comparison_df.loc['NC_011083:4482211:C:A', 'Category2_percent']
    assert 0 == comparison_df.loc['NC_011083:4482211:C:A', 'Category2_Unknown percent']
    assert 50 == comparison_df.loc['NC_011083:4482211:C:A', 'Category2_Present and Unknown percent']
    assert 'hgvs_gn:NC_011083:siiE:p.Arg1263Ser' == comparison_df.loc[
        'NC_011083:4482211:C:A', 'ID_HGVS_GN.p']
    assert 0 == comparison_df.loc['NC_011083:630556:G:A', 'Category1_percent']
    assert 100 == comparison_df.loc['NC_011083:630556:G:A', 'Category2_percent']
    assert 0 == comparison_df.loc['NC_011083:630556:G:A', 'Category1_percent']
    assert 0 == comparison_df.loc['NC_011083:630556:G:A', 'Category1_Unknown percent']
    assert 0 == comparison_df.loc['NC_011083:630556:G:A', 'Category1_Present and Unknown percent']
    assert 100 == comparison_df.loc['NC_011083:630556:G:A', 'Category2_percent']
    assert 0 == comparison_df.loc['NC_011083:630556:G:A', 'Category2_Unknown percent']
    assert 100 == comparison_df.loc['NC_011083:630556:G:A', 'Category2_Present and Unknown percent']
    assert 'hgvs_gn:NC_011083:SEHA_RS03545:p.Trp295*' == comparison_df.loc[
        'NC_011083:630556:G:A', 'ID_HGVS_GN.p']
    assert 0 == comparison_df.loc['NC_011083:1676762:CA:C', 'Category1_percent']
    assert 100 == comparison_df.loc['NC_011083:1676762:CA:C', 'Category1_Unknown percent']
    assert 100 == comparison_df.loc['NC_011083:1676762:CA:C', 'Category1_Present and Unknown percent']
    assert 50 == comparison_df.loc['NC_011083:1676762:CA:C', 'Category2_percent']
    assert 50 == comparison_df.loc['NC_011083:1676762:CA:C', 'Category2_Unknown percent']
    assert 100 == comparison_df.loc['NC_011083:1676762:CA:C', 'Category2_Present and Unknown percent']
    assert '<NA>' == comparison_df.loc[
        'NC_011083:1676762:CA:C', 'ID_HGVS_GN.p']

    # Test 2 categories: one of SH10-014 and one of SH14-001, SH14-014, threshold below
    sample_categories = [SampleSet([sample_sh10_014.id]), SampleSet([sample_sh14_001.id, sample_sh14_014.id])]
    comparison_df = mutations_summarizer.features_comparison(selected_samples=present_set,
                                                             sample_categories=sample_categories,
                                                             category_prefixes=['10', '14'],
                                                             category_samples_threshold=1,
                                                             unit='count')
    comparison_df = comparison_df.sort_index()
    assert comparison_df.index.name == 'Mutation'
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type', 'Total',
            '10_count', '10_Unknown count', '10_Present and Unknown count',
            '14_count', '14_Unknown count', '14_Present and Unknown count',
            '10_total', '14_total',
            'Annotation', 'Annotation_Impact',
            'Gene_Name', 'Gene_ID', 'Feature_Type', 'Transcript_BioType',
            'HGVS.c', 'HGVS.p', 'ID_HGVS.c', 'ID_HGVS.p', 'ID_HGVS_GN.c',
            'ID_HGVS_GN.p'] == list(comparison_df.columns)
    assert 177 == len(comparison_df)
    assert {3} == set(comparison_df['Total'].tolist())
    assert {1} == set(comparison_df['10_total'].tolist())
    assert {2} == set(comparison_df['14_total'].tolist())
    assert 1 == comparison_df.loc['NC_011083:140658:C:A', '10_count']
    assert 2 == comparison_df.loc['NC_011083:140658:C:A', '14_count']
    assert 'hgvs_gn:NC_011083:murF:p.Ala166Glu' == comparison_df.loc[
        'NC_011083:140658:C:A', 'ID_HGVS_GN.p']
    assert 1 == comparison_df.loc['NC_011083:4555461:T:TC', '10_count']
    assert 0 == comparison_df.loc['NC_011083:4555461:T:TC', '14_count']
    assert 'hgvs_gn:NC_011083:n.4555461_4555462insC' == comparison_df.loc[
        'NC_011083:4555461:T:TC', 'ID_HGVS_GN.c']
    assert 0 == comparison_df.loc['NC_011083:4482211:C:A', '10_count']
    assert 1 == comparison_df.loc['NC_011083:4482211:C:A', '14_count']
    assert 'hgvs_gn:NC_011083:siiE:p.Arg1263Ser' == comparison_df.loc[
        'NC_011083:4482211:C:A', 'ID_HGVS_GN.p']
    assert 0 == comparison_df.loc['NC_011083:630556:G:A', '10_count']
    assert 2 == comparison_df.loc['NC_011083:630556:G:A', '14_count']
    assert 'hgvs_gn:NC_011083:SEHA_RS03545:p.Trp295*' == comparison_df.loc[
        'NC_011083:630556:G:A', 'ID_HGVS_GN.p']

    # Test 2 categories: one of SH10-014 and one of SH14-001, SH14-014, threshold above
    sample_categories = [SampleSet([sample_sh10_014.id]), SampleSet([sample_sh14_001.id, sample_sh14_014.id])]
    comparison_df = mutations_summarizer.features_comparison(selected_samples=present_set,
                                                             sample_categories=sample_categories,
                                                             category_prefixes=['10', '14'],
                                                             category_samples_threshold=2,
                                                             unit='count')
    comparison_df = comparison_df.sort_index()
    assert comparison_df.index.name == 'Mutation'
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type', 'Total',
            '14_count', '14_Unknown count', '14_Present and Unknown count',
            '14_total',
            'Annotation', 'Annotation_Impact',
            'Gene_Name', 'Gene_ID', 'Feature_Type', 'Transcript_BioType',
            'HGVS.c', 'HGVS.p', 'ID_HGVS.c', 'ID_HGVS.p', 'ID_HGVS_GN.c',
            'ID_HGVS_GN.p'] == list(comparison_df.columns)
    assert 177 == len(comparison_df)
    assert {3} == set(comparison_df['Total'].tolist())
    assert {2} == set(comparison_df['14_total'].tolist())
    assert 2 == comparison_df.loc['NC_011083:140658:C:A', '14_count']
    assert 'hgvs_gn:NC_011083:murF:p.Ala166Glu' == comparison_df.loc[
        'NC_011083:140658:C:A', 'ID_HGVS_GN.p']
    assert 0 == comparison_df.loc['NC_011083:4555461:T:TC', '14_count']
    assert 'hgvs_gn:NC_011083:n.4555461_4555462insC' == comparison_df.loc[
        'NC_011083:4555461:T:TC', 'ID_HGVS_GN.c']
    assert 1 == comparison_df.loc['NC_011083:4482211:C:A', '14_count']
    assert 'hgvs_gn:NC_011083:siiE:p.Arg1263Ser' == comparison_df.loc[
        'NC_011083:4482211:C:A', 'ID_HGVS_GN.p']
    assert 2 == comparison_df.loc['NC_011083:630556:G:A', '14_count']
    assert 'hgvs_gn:NC_011083:SEHA_RS03545:p.Trp295*' == comparison_df.loc[
        'NC_011083:630556:G:A', 'ID_HGVS_GN.p']
