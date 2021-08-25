import warnings

import pandas as pd

from genomics_data_index.api.query.GenomicsDataIndex import GenomicsDataIndex
from genomics_data_index.api.query.features.MutationFeaturesFromIndexComparator import \
    MutationFeaturesFromIndexComparator
from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.db import Sample
from genomics_data_index.test.integration import snippy_all_dataframes


def test_summary_all(loaded_database_genomic_data_store: GenomicsDataIndex):
    db = loaded_database_genomic_data_store.connection.database
    all_sample_ids = {s.id for s in db.get_session().query(Sample).all()}

    dfA = pd.read_csv(snippy_all_dataframes['SampleA'], sep='\t')
    dfB = pd.read_csv(snippy_all_dataframes['SampleB'], sep='\t')
    dfC = pd.read_csv(snippy_all_dataframes['SampleC'], sep='\t')
    expected_df = pd.concat([dfA, dfB, dfC])
    expected_df = expected_df.groupby('Mutation').agg({
        'Sequence': 'first',
        'Position': 'first',
        'Deletion': 'first',
        'Insertion': 'first',
        'Type': 'first',
        'Mutation': 'count',
    }).rename(columns={'Mutation': 'Count'}).sort_index()
    expected_df['Total'] = 9
    expected_df['Percent'] = 100 * (expected_df['Count'] / expected_df['Total'])

    present_set = SampleSet(all_sample_ids)
    mutations_summarizer = MutationFeaturesFromIndexComparator(connection=loaded_database_genomic_data_store.connection,
                                                               ignore_annotations=True)

    mutations_df = mutations_summarizer.summary(present_set)
    mutations_df['Percent'] = mutations_df['Percent'].astype(int)  # Convert to int for easier comparison
    mutations_df = mutations_df.sort_index()

    assert len(expected_df) == len(mutations_df)
    assert list(expected_df.columns) == list(mutations_df.columns)
    assert list(expected_df.index) == list(mutations_df.index)
    assert list(expected_df['Deletion']) == list(mutations_df['Deletion'])
    assert list(expected_df['Count']) == list(mutations_df['Count'])
    assert list(expected_df['Total']) == list(mutations_df['Total'])
    assert list(expected_df['Type']) == list(mutations_df['Type'])
    assert 22 == mutations_df.loc['reference:619:G:C', 'Percent']
    assert 11 == mutations_df.loc['reference:461:AAAT:G', 'Percent']

    # Test with unknown
    mutations_summarizer = MutationFeaturesFromIndexComparator(connection=loaded_database_genomic_data_store.connection,
                                                               ignore_annotations=True, include_unknown=True)
    mutations_df = mutations_summarizer.summary(present_set)
    mutations_df['Percent'] = mutations_df['Percent'].astype(int)  # Convert to int for easier comparison
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

    # Test only include unknown
    mutations_summarizer = MutationFeaturesFromIndexComparator(connection=loaded_database_genomic_data_store.connection,
                                                               ignore_annotations=True, include_unknown=True,
                                                               include_present=False)
    mutations_df = mutations_summarizer.summary(present_set)
    mutations_df['Percent'] = mutations_df['Percent'].astype(int)  # Convert to int for easier comparison
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
                                                               ignore_annotations=True, id_type='spdi')
    mutations_df = mutations_summarizer.summary(present_set)
    mutations_df['Percent'] = mutations_df['Percent'].astype(int)  # Convert to int for easier comparison
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

    # Test with different id type include unknown
    mutations_summarizer = MutationFeaturesFromIndexComparator(connection=loaded_database_genomic_data_store.connection,
                                                               ignore_annotations=True, id_type='spdi',
                                                               include_unknown=True)
    mutations_df = mutations_summarizer.summary(present_set)
    mutations_df['Percent'] = mutations_df['Percent'].astype(int)  # Convert to int for easier comparison
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
                                                               ignore_annotations=True)

    dfA = pd.read_csv(snippy_all_dataframes['SampleA'], sep='\t')
    dfB = pd.read_csv(snippy_all_dataframes['SampleB'], sep='\t')
    dfC = pd.read_csv(snippy_all_dataframes['SampleC'], sep='\t')

    # Unique to A
    present_set = SampleSet({sampleA.id})
    other_set = SampleSet(all_sample_ids - {sampleA.id})
    mutations_df = mutations_summarizer.unique_summary(present_set, other_set=other_set).sort_index()

    expected_df = dfA
    expected_df = expected_df.groupby('Mutation').agg({
        'Sequence': 'first',
        'Position': 'first',
        'Deletion': 'first',
        'Insertion': 'first',
        'Type': 'first',
        'Mutation': 'count',
    }).rename(columns={'Mutation': 'Count'}).sort_index()
    expected_df['Total'] = 1
    expected_df['Percent'] = 100 * (expected_df['Count'] / expected_df['Total'])

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

    dfAC = pd.concat([dfA, dfC])
    expected_df = dfB[~dfB['Mutation'].isin(list(dfAC['Mutation']))]
    expected_df = expected_df.groupby('Mutation').agg({
        'Sequence': 'first',
        'Position': 'first',
        'Deletion': 'first',
        'Insertion': 'first',
        'Type': 'first',
        'Mutation': 'count',
    }).rename(columns={'Mutation': 'Count'}).sort_index()
    expected_df['Total'] = 1
    expected_df['Percent'] = 100 * (expected_df['Count'] / expected_df['Total'])

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

    dfBC = pd.concat([dfB, dfC])
    expected_df = dfBC[~dfBC['Mutation'].isin(list(dfA['Mutation']))]
    expected_df = expected_df.groupby('Mutation').agg({
        'Sequence': 'first',
        'Position': 'first',
        'Deletion': 'first',
        'Insertion': 'first',
        'Type': 'first',
        'Mutation': 'count',
    }).rename(columns={'Mutation': 'Count'}).sort_index()
    expected_df['Total'] = 2
    expected_df['Percent'] = 100 * (expected_df['Count'] / expected_df['Total'])

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
        ignore_annotations=False)

    sample_sh14_001 = db.get_session().query(Sample).filter(Sample.name == 'SH14-001').one()
    sample_sh14_014 = db.get_session().query(Sample).filter(Sample.name == 'SH14-014').one()
    sample_sh10_014 = db.get_session().query(Sample).filter(Sample.name == 'SH10-014').one()
    three_samples = {sample_sh14_001.id, sample_sh14_014.id, sample_sh10_014.id}

    present_set = SampleSet(three_samples)
    mutations_df = mutations_summarizer.summary(present_set)

    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type',
            'Count', 'Total', 'Percent', 'Annotation', 'Annotation_Impact',
            'Gene_Name', 'Gene_ID', 'Feature_Type', 'Transcript_BioType',
            'HGVS.c', 'HGVS.p', 'ID_HGVS.c', 'ID_HGVS.p', 'ID_HGVS_GN.c', 'ID_HGVS_GN.p'] == list(mutations_df.columns)
    assert 177 == len(mutations_df)
    mutations_df['Percent'] = mutations_df['Percent'].astype(int)  # easier to compare percents in assert

    # missense variant (3/3)
    assert ['NC_011083', 140658, 'C', 'A', 'SNP', 3, 3, 100,
            'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
            'c.497C>A', 'p.Ala166Glu',
            'hgvs:NC_011083:SEHA_RS01180:c.497C>A', 'hgvs:NC_011083:SEHA_RS01180:p.Ala166Glu',
            'hgvs_gn:NC_011083:murF:c.497C>A', 'hgvs_gn:NC_011083:murF:p.Ala166Glu'] == list(
        mutations_df.loc['NC_011083:140658:C:A'])

    # Intergenic variant (1/3)
    assert ['NC_011083', 4555461, 'T', 'TC', 'INDEL', 1, 3, 33,
            'intergenic_region', 'MODIFIER', 'SEHA_RS22510-SEHA_RS26685', 'SEHA_RS22510-SEHA_RS26685',
            'intergenic_region', 'NA',
            'n.4555461_4555462insC', 'NA',
            'hgvs:NC_011083:n.4555461_4555462insC', 'NA',
            'hgvs_gn:NC_011083:n.4555461_4555462insC', 'NA'] == list(
        mutations_df.loc['NC_011083:4555461:T:TC'].fillna('NA'))

    # MNP variant (1/3)
    assert ['NC_011083', 3535698, 'GCC', 'CAT', 'MNP', 2, 3, 66,
            'missense_variant', 'MODERATE', 'oadA', 'SEHA_RS17780',
            'transcript', 'protein_coding',
            'c.544_546delGGCinsATG', 'p.Gly182Met',
            'hgvs:NC_011083:SEHA_RS17780:c.544_546delGGCinsATG', 'hgvs:NC_011083:SEHA_RS17780:p.Gly182Met',
            'hgvs_gn:NC_011083:oadA:c.544_546delGGCinsATG', 'hgvs_gn:NC_011083:oadA:p.Gly182Met'] == list(
        mutations_df.loc['NC_011083:3535698:GCC:CAT'])


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
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type',
            'Total', 'All_count', 'All_total'] == comparison_df.columns.tolist()
    assert {9} == set(comparison_df['Total'].tolist())
    assert {9} == set(comparison_df['All_total'].tolist())
    assert 2 == comparison_df.loc['reference:619:G:C', 'All_count']
    assert 'SNP' == comparison_df.loc['reference:619:G:C', 'Type']
    assert 1 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'All_count']
    assert 'INDEL' == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'Type']
    assert 2 == comparison_df.loc['reference:4693:C:CGA', 'All_count']
    assert 'INDEL' == comparison_df.loc['reference:4693:C:CGA', 'Type']

    # Test two categories, one of A and one of BC
    sample_categories = [SampleSet([sampleA.id]), SampleSet([sampleB.id, sampleC.id])]
    comparison_df = mutations_summarizer.features_comparison(selected_samples=present_set,
                                                             sample_categories=sample_categories,
                                                             category_prefixes=['A', 'BC'],
                                                             unit='count')
    comparison_df = comparison_df.sort_index()
    assert comparison_df.index.name == 'Mutation'
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type',
            'Total', 'A_count', 'BC_count', 'A_total', 'BC_total'] == comparison_df.columns.tolist()
    assert {9} == set(comparison_df['Total'].tolist())
    assert {1} == set(comparison_df['A_total'].tolist())
    assert {2} == set(comparison_df['BC_total'].tolist())
    assert 0 == comparison_df.loc['reference:619:G:C', 'A_count']
    assert 2 == comparison_df.loc['reference:619:G:C', 'BC_count']
    assert 1 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'A_count']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'BC_count']
    assert 0 == comparison_df.loc['reference:4693:C:CGA', 'A_count']
    assert 2 == comparison_df.loc['reference:4693:C:CGA', 'BC_count']

    # Test two categories, one of AB and one of C
    sample_categories = [SampleSet([sampleA.id, sampleB.id]), SampleSet([sampleC.id])]
    comparison_df = mutations_summarizer.features_comparison(selected_samples=present_set,
                                                             sample_categories=sample_categories,
                                                             category_prefixes=['AB', 'C'],
                                                             unit='count')
    comparison_df = comparison_df.sort_index()
    assert comparison_df.index.name == 'Mutation'
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type',
            'Total', 'AB_count', 'C_count', 'AB_total', 'C_total'] == comparison_df.columns.tolist()
    assert {9} == set(comparison_df['Total'].tolist())
    assert {2} == set(comparison_df['AB_total'].tolist())
    assert {1} == set(comparison_df['C_total'].tolist())
    assert 1 == comparison_df.loc['reference:619:G:C', 'AB_count']
    assert 1 == comparison_df.loc['reference:619:G:C', 'C_count']
    assert 1 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'AB_count']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'C_count']
    assert 1 == comparison_df.loc['reference:4693:C:CGA', 'AB_count']
    assert 1 == comparison_df.loc['reference:4693:C:CGA', 'C_count']
    assert 1 == comparison_df.loc['reference:528:C:CAG', 'AB_count']
    assert 0 == comparison_df.loc['reference:528:C:CAG', 'C_count']

    # Test three categories: A, B, and C, and total out of only these 3
    sample_categories = [SampleSet([sampleA.id]), SampleSet([sampleB.id]), SampleSet([sampleC.id])]
    selected_samples = SampleSet([sampleA.id, sampleB.id, sampleC.id])
    comparison_df = mutations_summarizer.features_comparison(selected_samples=selected_samples,
                                                             sample_categories=sample_categories,
                                                             category_prefixes=['A', 'B', 'C'],
                                                             unit='count')
    comparison_df = comparison_df.sort_index()
    assert comparison_df.index.name == 'Mutation'
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type',
            'Total',
            'A_count', 'B_count', 'C_count',
            'A_total', 'B_total', 'C_total'] == comparison_df.columns.tolist()
    assert {3} == set(comparison_df['Total'].tolist())
    assert {1} == set(comparison_df['A_total'].tolist())
    assert {1} == set(comparison_df['B_total'].tolist())
    assert {1} == set(comparison_df['C_total'].tolist())
    assert 0 == comparison_df.loc['reference:619:G:C', 'A_count']
    assert 1 == comparison_df.loc['reference:619:G:C', 'B_count']
    assert 1 == comparison_df.loc['reference:619:G:C', 'C_count']
    assert 1 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'A_count']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'B_count']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'C_count']
    assert 0 == comparison_df.loc['reference:4693:C:CGA', 'A_count']
    assert 1 == comparison_df.loc['reference:4693:C:CGA', 'B_count']
    assert 1 == comparison_df.loc['reference:4693:C:CGA', 'C_count']
    assert 0 == comparison_df.loc['reference:528:C:CAG', 'A_count']
    assert 1 == comparison_df.loc['reference:528:C:CAG', 'B_count']
    assert 0 == comparison_df.loc['reference:528:C:CAG', 'C_count']

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
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type',
            'Total', 'A_percent', 'BC_percent', 'A_total', 'BC_total'] == comparison_df.columns.tolist()
    assert {9} == set(comparison_df['Total'].tolist())
    assert {1} == set(comparison_df['A_total'].tolist())
    assert {2} == set(comparison_df['BC_total'].tolist())
    assert 0 == comparison_df.loc['reference:619:G:C', 'A_percent']
    assert 100 == comparison_df.loc['reference:619:G:C', 'BC_percent']
    assert 100 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'A_percent']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'BC_percent']
    assert 0 == comparison_df.loc['reference:4693:C:CGA', 'A_percent']
    assert 100 == comparison_df.loc['reference:4693:C:CGA', 'BC_percent']

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
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type',
            'Total',
            'Category1_percent', 'Category2_percent',
            'Category1_total', 'Category2_total'] == comparison_df.columns.tolist()
    assert {9} == set(comparison_df['Total'].tolist())
    assert {1} == set(comparison_df['Category1_total'].tolist())
    assert {2} == set(comparison_df['Category2_total'].tolist())
    assert 0 == comparison_df.loc['reference:619:G:C', 'Category1_percent']
    assert 100 == comparison_df.loc['reference:619:G:C', 'Category2_percent']
    assert 100 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'Category1_percent']
    assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'Category2_percent']
    assert 0 == comparison_df.loc['reference:4693:C:CGA', 'Category1_percent']
    assert 100 == comparison_df.loc['reference:4693:C:CGA', 'Category2_percent']


def test_features_comparison_annotations(loaded_database_genomic_data_store_annotations: GenomicsDataIndex):
    db = loaded_database_genomic_data_store_annotations.connection.database

    sample_sh14_001 = db.get_session().query(Sample).filter(Sample.name == 'SH14-001').one()
    sample_sh14_014 = db.get_session().query(Sample).filter(Sample.name == 'SH14-014').one()
    sample_sh10_014 = db.get_session().query(Sample).filter(Sample.name == 'SH10-014').one()
    three_samples = {sample_sh14_001.id, sample_sh14_014.id, sample_sh10_014.id}

    present_set = SampleSet(three_samples)
    mutations_summarizer = MutationFeaturesFromIndexComparator(
        connection=loaded_database_genomic_data_store_annotations.connection,
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
            'All_count',
            'All_total',
            'Annotation', 'Annotation_Impact',
            'Gene_Name', 'Gene_ID', 'Feature_Type', 'Transcript_BioType',
            'HGVS.c', 'HGVS.p', 'ID_HGVS.c', 'ID_HGVS.p', 'ID_HGVS_GN.c',
            'ID_HGVS_GN.p'] == list(comparison_df.columns)
    assert 177 == len(comparison_df)
    assert {3} == set(comparison_df['Total'].tolist())
    assert {3} == set(comparison_df['All_total'].tolist())
    assert 3 == comparison_df.loc['NC_011083:140658:C:A', 'All_count']
    assert 'hgvs_gn:NC_011083:murF:p.Ala166Glu' == comparison_df.loc[
        'NC_011083:140658:C:A', 'ID_HGVS_GN.p']
    assert 1 == comparison_df.loc['NC_011083:4555461:T:TC', 'All_count']
    assert 'hgvs_gn:NC_011083:n.4555461_4555462insC' == comparison_df.loc[
        'NC_011083:4555461:T:TC', 'ID_HGVS_GN.c']

    # Test 2 categories: one of SH10-014 and one of SH14-001, SH14-014
    sample_categories = [SampleSet([sample_sh10_014.id]), SampleSet([sample_sh14_001.id, sample_sh14_014.id])]
    comparison_df = mutations_summarizer.features_comparison(selected_samples=present_set,
                                                             sample_categories=sample_categories,
                                                             category_prefixes=['10', '14'],
                                                             unit='count')
    comparison_df = comparison_df.sort_index()
    assert comparison_df.index.name == 'Mutation'
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type', 'Total',
            '10_count', '14_count',
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

    # Test 2 categories defaults: one of SH10-014 and one of SH14-001, SH14-014
    sample_categories = [SampleSet([sample_sh10_014.id]), SampleSet([sample_sh14_001.id, sample_sh14_014.id])]
    comparison_df = mutations_summarizer.features_comparison(selected_samples=present_set,
                                                             sample_categories=sample_categories)
    comparison_df = comparison_df.sort_index()
    comparison_df['Category1_percent'] = comparison_df['Category1_percent'].astype(
        int)  # Convert to int for easier comparison
    comparison_df['Category2_percent'] = comparison_df['Category2_percent'].astype(
        int)  # Convert to int for easier comparison
    assert comparison_df.index.name == 'Mutation'
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type', 'Total',
            'Category1_percent', 'Category2_percent',
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
    assert 100 == comparison_df.loc['NC_011083:140658:C:A', 'Category2_percent']
    assert 'hgvs_gn:NC_011083:murF:p.Ala166Glu' == comparison_df.loc[
        'NC_011083:140658:C:A', 'ID_HGVS_GN.p']
    assert 100 == comparison_df.loc['NC_011083:4555461:T:TC', 'Category1_percent']
    assert 0 == comparison_df.loc['NC_011083:4555461:T:TC', 'Category2_percent']
    assert 'hgvs_gn:NC_011083:n.4555461_4555462insC' == comparison_df.loc[
        'NC_011083:4555461:T:TC', 'ID_HGVS_GN.c']
    assert 0 == comparison_df.loc['NC_011083:4482211:C:A', 'Category1_percent']
    assert 50 == comparison_df.loc['NC_011083:4482211:C:A', 'Category2_percent']
    assert 'hgvs_gn:NC_011083:siiE:p.Arg1263Ser' == comparison_df.loc[
        'NC_011083:4482211:C:A', 'ID_HGVS_GN.p']
    assert 0 == comparison_df.loc['NC_011083:630556:G:A', 'Category1_percent']
    assert 100 == comparison_df.loc['NC_011083:630556:G:A', 'Category2_percent']
    assert 'hgvs_gn:NC_011083:SEHA_RS03545:p.Trp295*' == comparison_df.loc[
        'NC_011083:630556:G:A', 'ID_HGVS_GN.p']

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
            '10_count', '14_count',
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
            '14_count',
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
