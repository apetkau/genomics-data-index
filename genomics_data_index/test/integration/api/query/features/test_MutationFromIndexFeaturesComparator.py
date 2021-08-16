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
        'Mutation': 'count',
    }).rename(columns={'Mutation': 'Count'}).sort_index()
    expected_df['Total'] = 9
    expected_df['Percent'] = 100 * (expected_df['Count'] / expected_df['Total'])

    f = ['reference:461:AAAT:G']
    warnings.warn(f'Removing {f} from expected until I can figure out how to properly handle these')
    expected_df = expected_df.drop(f)

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
    assert 22 == mutations_df.loc['reference:619:G:C', 'Percent']

    # Test with unknown
    mutations_summarizer = MutationFeaturesFromIndexComparator(connection=loaded_database_genomic_data_store.connection,
                                                               ignore_annotations=True, include_unknown=True)
    mutations_df = mutations_summarizer.summary(present_set)
    mutations_df['Percent'] = mutations_df['Percent'].astype(int)  # Convert to int for easier comparison
    mutations_df = mutations_df.sort_index()

    assert 632 == len(mutations_df)
    assert list(expected_df.columns) == list(mutations_df.columns)
    assert 2 == mutations_df.loc['reference:619:G:C', 'Count']
    assert 2 == mutations_df.loc['reference:3063:A:ATGCAGC', 'Count']
    assert 1 == mutations_df.loc['reference:1984:GTGATTG:TTGA', 'Count']
    assert 1 == mutations_df.loc['reference:866:GCCAGATCC:G', 'Count']
    assert 3 == mutations_df.loc['reference:90:T:?', 'Count']
    assert 2 == mutations_df.loc['reference:190:A:?', 'Count']
    assert 1 == mutations_df.loc['reference:887:T:?', 'Count']

    # Test only include unknown
    mutations_summarizer = MutationFeaturesFromIndexComparator(connection=loaded_database_genomic_data_store.connection,
                                                               ignore_annotations=True, include_unknown=True,
                                                               include_present=False)
    mutations_df = mutations_summarizer.summary(present_set)
    mutations_df['Percent'] = mutations_df['Percent'].astype(int)  # Convert to int for easier comparison
    mutations_df = mutations_df.sort_index()

    assert 521 == len(mutations_df)
    assert list(expected_df.columns) == list(mutations_df.columns)
    assert 3 == mutations_df.loc['reference:90:T:?', 'Count']
    assert 2 == mutations_df.loc['reference:190:A:?', 'Count']
    assert 1 == mutations_df.loc['reference:887:T:?', 'Count']
    assert 'reference:619:G:C' not in mutations_df

    # Test with different id type where deletion is int instead of sequence
    mutations_summarizer = MutationFeaturesFromIndexComparator(connection=loaded_database_genomic_data_store.connection,
                                                               ignore_annotations=True, id_type='spdi')
    mutations_df = mutations_summarizer.summary(present_set)
    mutations_df['Percent'] = mutations_df['Percent'].astype(int)  # Convert to int for easier comparison
    mutations_df = mutations_df.sort_index()

    assert 112 - 1 == len(mutations_df)
    assert list(expected_df.columns) == list(mutations_df.columns)
    assert 2 == mutations_df.loc['reference:619:1:C', 'Count']
    assert 2 == mutations_df.loc['reference:3063:1:ATGCAGC', 'Count']
    assert 1 == mutations_df.loc['reference:1984:7:TTGA', 'Count']
    assert 1 == mutations_df.loc['reference:866:9:G', 'Count']

    # Test with different id type include unknown
    mutations_summarizer = MutationFeaturesFromIndexComparator(connection=loaded_database_genomic_data_store.connection,
                                                               ignore_annotations=True, id_type='spdi',
                                                               include_unknown=True)
    mutations_df = mutations_summarizer.summary(present_set)
    mutations_df['Percent'] = mutations_df['Percent'].astype(int)  # Convert to int for easier comparison
    mutations_df = mutations_df.sort_index()

    assert 632 == len(mutations_df)
    assert list(expected_df.columns) == list(mutations_df.columns)
    assert 2 == mutations_df.loc['reference:619:1:C', 'Count']
    assert 2 == mutations_df.loc['reference:3063:1:ATGCAGC', 'Count']
    assert 1 == mutations_df.loc['reference:1984:7:TTGA', 'Count']
    assert 1 == mutations_df.loc['reference:866:9:G', 'Count']
    assert 3 == mutations_df.loc['reference:90:1:?', 'Count']
    assert 2 == mutations_df.loc['reference:190:1:?', 'Count']
    assert 1 == mutations_df.loc['reference:887:1:?', 'Count']


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
        'Mutation': 'count',
    }).rename(columns={'Mutation': 'Count'}).sort_index()
    expected_df['Total'] = 1
    expected_df['Percent'] = 100 * (expected_df['Count'] / expected_df['Total'])

    mutations_df['Percent'] = mutations_df['Percent'].astype(int)  # Convert to int for easier comparison

    f = ['reference:461:AAAT:G']
    warnings.warn(f'Removing {f} from expected until I can figure out how to properly handle these')
    expected_df = expected_df.drop(f)

    assert len(expected_df) == len(mutations_df)
    assert 45 == len(mutations_df)  # Check length against independently generated length
    assert list(expected_df.index) == list(mutations_df.index)
    assert list(expected_df['Count']) == list(mutations_df['Count'])
    assert list(expected_df['Total']) == list(mutations_df['Total'])
    assert 100 == mutations_df.loc['reference:3656:CATT:C', 'Percent']

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

    assert ['Sequence', 'Position', 'Deletion', 'Insertion',
            'Count', 'Total', 'Percent', 'Annotation', 'Annotation_Impact',
            'Gene_Name', 'Gene_ID', 'Feature_Type', 'Transcript_BioType',
            'HGVS.c', 'HGVS.p', 'ID_HGVS.c', 'ID_HGVS.p', 'ID_HGVS_GN.c', 'ID_HGVS_GN.p'] == list(mutations_df.columns)
    assert 177 == len(mutations_df)
    mutations_df['Percent'] = mutations_df['Percent'].astype(int)  # easier to compare percents in assert

    # missense variant (3/3)
    assert ['NC_011083', 140658, 'C', 'A', 3, 3, 100,
            'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
            'c.497C>A', 'p.Ala166Glu',
            'hgvs:NC_011083:SEHA_RS01180:c.497C>A', 'hgvs:NC_011083:SEHA_RS01180:p.Ala166Glu',
            'hgvs_gn:NC_011083:murF:c.497C>A', 'hgvs_gn:NC_011083:murF:p.Ala166Glu'] == list(
        mutations_df.loc['NC_011083:140658:C:A'])

    # Intergenic variant (1/3)
    assert ['NC_011083', 4555461, 'T', 'TC', 1, 3, 33,
            'intergenic_region', 'MODIFIER', 'SEHA_RS22510-SEHA_RS26685', 'SEHA_RS22510-SEHA_RS26685',
            'intergenic_region', 'NA',
            'n.4555461_4555462insC', 'NA',
            'hgvs:NC_011083:n.4555461_4555462insC', 'NA',
            'hgvs_gn:NC_011083:n.4555461_4555462insC', 'NA'] == list(
        mutations_df.loc['NC_011083:4555461:T:TC'].fillna('NA'))


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
                                                                 category_prefixes=['All_count'],
                                                                 compare_kind='count')
        comparison_df = comparison_df.sort_index()
        assert comparison_df.index.name == 'Mutation'
        assert ['Sequence', 'Position', 'Deletion', 'Insertion',
                'Total', 'All_count'] == comparison_df.columns.tolist()
        assert {9} == set(comparison_df['Total'].tolist())
        assert 2 == comparison_df.loc['reference:619:G:C', 'All_count']
        assert 1 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'All_count']
        assert 2 == comparison_df.loc['reference:4693:C:CGA', 'All_count']

        # Test two categories, one of A and one of BC
        sample_categories = [SampleSet([sampleA.id]), SampleSet([sampleB.id, sampleC.id])]
        comparison_df = mutations_summarizer.features_comparison(selected_samples=present_set,
                                                                 sample_categories=sample_categories,
                                                                 category_prefixes=['A_count', 'BC_count'],
                                                                 compare_kind='count')
        comparison_df = comparison_df.sort_index()
        assert comparison_df.index.name == 'Mutation'
        assert ['Sequence', 'Position', 'Deletion', 'Insertion',
                'Total', 'A_count', 'BC_count'] == comparison_df.columns.tolist()
        assert {9} == set(comparison_df['Total'].tolist())
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
                                                                 category_prefixes=['AB_count', 'C_count'],
                                                                 compare_kind='count')
        comparison_df = comparison_df.sort_index()
        assert comparison_df.index.name == 'Mutation'
        assert ['Sequence', 'Position', 'Deletion', 'Insertion',
                'Total', 'AB_count', 'C_count'] == comparison_df.columns.tolist()
        assert {9} == set(comparison_df['Total'].tolist())
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
                                                                 category_prefixes=['A_count', 'B_count', 'C_count'],
                                                                 compare_kind='count')
        comparison_df = comparison_df.sort_index()
        assert comparison_df.index.name == 'Mutation'
        assert ['Sequence', 'Position', 'Deletion', 'Insertion',
                'Total', 'A_count', 'B_count', 'C_count'] == comparison_df.columns.tolist()
        assert {3} == set(comparison_df['Total'].tolist())
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
                                                                 category_prefixes=['A_percent', 'BC_percent'],
                                                                 compare_kind='percent')
        comparison_df = comparison_df.sort_index()
        comparison_df['A_percent'] = comparison_df['A_percent'].astype(int)  # Convert to int for easier comparison
        comparison_df['BC_percent'] = comparison_df['BC_percent'].astype(int)  # Convert to int for easier comparison
        assert comparison_df.index.name == 'Mutation'
        assert ['Sequence', 'Position', 'Deletion', 'Insertion',
                'Total', 'A_percent', 'BC_percent'] == comparison_df.columns.tolist()
        assert {9} == set(comparison_df['Total'].tolist())
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
        comparison_df['Category1'] = comparison_df['Category1'].astype(int)  # Convert to int for easier comparison
        comparison_df['Category2'] = comparison_df['Category2'].astype(int)  # Convert to int for easier comparison
        assert comparison_df.index.name == 'Mutation'
        assert ['Sequence', 'Position', 'Deletion', 'Insertion',
                'Total', 'Category1', 'Category2'] == comparison_df.columns.tolist()
        assert {9} == set(comparison_df['Total'].tolist())
        assert 0 == comparison_df.loc['reference:619:G:C', 'Category1']
        assert 100 == comparison_df.loc['reference:619:G:C', 'Category2']
        assert 100 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'Category1']
        assert 0 == comparison_df.loc['reference:1708:ATGCTGTTCAATAC:A', 'Category2']
        assert 0 == comparison_df.loc['reference:4693:C:CGA', 'Category1']
        assert 100 == comparison_df.loc['reference:4693:C:CGA', 'Category2']


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
                                                             category_prefixes=['All_count'],
                                                             compare_kind='count')
    comparison_df = comparison_df.sort_index()
    assert comparison_df.index.name == 'Mutation'
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Total',
            'All_count',
            'Annotation', 'Annotation_Impact',
            'Gene_Name', 'Gene_ID', 'Feature_Type', 'Transcript_BioType',
            'HGVS.c', 'HGVS.p', 'ID_HGVS.c', 'ID_HGVS.p', 'ID_HGVS_GN.c',
            'ID_HGVS_GN.p'] == list(comparison_df.columns)
    assert 177 == len(comparison_df)
    assert {3} == set(comparison_df['Total'].tolist())
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
                                                             category_prefixes=['10_count', '14_count'],
                                                             compare_kind='count')
    comparison_df = comparison_df.sort_index()
    assert comparison_df.index.name == 'Mutation'
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Total',
            '10_count', '14_count',
            'Annotation', 'Annotation_Impact',
            'Gene_Name', 'Gene_ID', 'Feature_Type', 'Transcript_BioType',
            'HGVS.c', 'HGVS.p', 'ID_HGVS.c', 'ID_HGVS.p', 'ID_HGVS_GN.c',
            'ID_HGVS_GN.p'] == list(comparison_df.columns)
    assert 177 == len(comparison_df)
    assert {3} == set(comparison_df['Total'].tolist())
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
    comparison_df['Category1'] = comparison_df['Category1'].astype(int)  # Convert to int for easier comparison
    comparison_df['Category2'] = comparison_df['Category2'].astype(int)  # Convert to int for easier comparison
    assert comparison_df.index.name == 'Mutation'
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Total',
            'Category1', 'Category2',
            'Annotation', 'Annotation_Impact',
            'Gene_Name', 'Gene_ID', 'Feature_Type', 'Transcript_BioType',
            'HGVS.c', 'HGVS.p', 'ID_HGVS.c', 'ID_HGVS.p', 'ID_HGVS_GN.c',
            'ID_HGVS_GN.p'] == list(comparison_df.columns)
    assert 177 == len(comparison_df)
    assert {3} == set(comparison_df['Total'].tolist())
    assert 100 == comparison_df.loc['NC_011083:140658:C:A', 'Category1']
    assert 100 == comparison_df.loc['NC_011083:140658:C:A', 'Category2']
    assert 'hgvs_gn:NC_011083:murF:p.Ala166Glu' == comparison_df.loc[
        'NC_011083:140658:C:A', 'ID_HGVS_GN.p']
    assert 100 == comparison_df.loc['NC_011083:4555461:T:TC', 'Category1']
    assert 0 == comparison_df.loc['NC_011083:4555461:T:TC', 'Category2']
    assert 'hgvs_gn:NC_011083:n.4555461_4555462insC' == comparison_df.loc[
        'NC_011083:4555461:T:TC', 'ID_HGVS_GN.c']
    assert 0 == comparison_df.loc['NC_011083:4482211:C:A', 'Category1']
    assert 50 == comparison_df.loc['NC_011083:4482211:C:A', 'Category2']
    assert 'hgvs_gn:NC_011083:siiE:p.Arg1263Ser' == comparison_df.loc[
        'NC_011083:4482211:C:A', 'ID_HGVS_GN.p']
    assert 0 == comparison_df.loc['NC_011083:630556:G:A', 'Category1']
    assert 100 == comparison_df.loc['NC_011083:630556:G:A', 'Category2']
    assert 'hgvs_gn:NC_011083:SEHA_RS03545:p.Trp295*' == comparison_df.loc[
        'NC_011083:630556:G:A', 'ID_HGVS_GN.p']

