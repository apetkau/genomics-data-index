import shutil
import tempfile
from pathlib import Path
from tempfile import TemporaryDirectory

import pytest

from genomics_data_index.api.query.GenomicsDataIndex import GenomicsDataIndex
from genomics_data_index.configuration.Project import Project
from genomics_data_index.storage.io.mutation.NucleotideSampleDataPackage import NucleotideSampleDataPackage
from genomics_data_index.storage.io.processor.SerialSampleFilesProcessor import SerialSampleFilesProcessor
from genomics_data_index.storage.service import EntityExistsError
from genomics_data_index.test.integration import sample_dirs, reference_file


def test_get_reference_tree(loaded_database_genomic_data_store_with_tree: GenomicsDataIndex):
    tree = loaded_database_genomic_data_store_with_tree.reference_tree('genome')
    assert tree is not None

    with pytest.raises(EntityExistsError) as execinfo:
        loaded_database_genomic_data_store_with_tree.reference_tree('invalid_reference')

    assert 'No reference genome with name=[invalid_reference]' in str(execinfo.value)


def test_mlst_schemes(loaded_database_genomic_data_store: GenomicsDataIndex):
    mlst_schemes = loaded_database_genomic_data_store.mlst_schemes()
    assert 3 == len(mlst_schemes)
    assert isinstance(mlst_schemes, list)
    assert {'lmonocytogenes', 'ecoli', 'campylobacter'} == set(mlst_schemes)


def test_summaries_loaded_data(loaded_database_genomic_data_store: GenomicsDataIndex):
    gds = loaded_database_genomic_data_store

    # Samples
    assert 9 == gds.count_samples()
    assert {'2014C-3598', '2014C-3599', '2014D-0067', '2014D-0068',
            'CFSAN002349', 'CFSAN023463',
            'SampleA', 'SampleB', 'SampleC'} == set(gds.sample_names())

    # References
    assert 1 == gds.count_references()
    assert ['genome'] == gds.reference_names()

    assert 112 == gds.count_mutations('genome')
    assert 112 + 440 == gds.count_mutations('genome', include_unknown=True)

    # Mutations ignore unknown spdi
    ms = gds.mutations_summary('genome', id_type='spdi', ignore_annotations=True,
                               include_unknown_samples=False)
    assert 112 == len(ms)
    assert 'Mutation' == ms.index.name
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type',
            'Count', 'Unknown Count', 'Present and Unknown Count', 'Total',
            'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(ms.columns)

    ## Convert percent to int to make it easier to compare in assert statements
    ms['Percent'] = ms['Percent'].astype(int)
    ms['Unknown Count'] = ms['Unknown Count'].fillna(-1).astype(int)
    ms['Present and Unknown Count'] = ms['Present and Unknown Count'].fillna(-1).astype(int)
    ms['Unknown Percent'] = ms['Unknown Percent'].fillna(-1).astype(int)
    ms['Present and Unknown Percent'] = ms['Present and Unknown Percent'].fillna(-1).astype(int)

    assert ['reference', 839, 1, 'G', 'SNP', 2, -1, -1, 3, 66, -1, -1] == ms.loc['reference:839:1:G'].values.tolist()
    assert ['reference', 866, 9, 'G', 'INDEL', 1, -1, -1, 3, 33, -1, -1] == ms.loc['reference:866:9:G'].values.tolist()
    assert ['reference', 1048, 1, 'G', 'SNP', 1, -1, -1, 3, 33, -1, -1] == ms.loc['reference:1048:1:G'].values.tolist()
    assert ['reference', 3897, 5, 'G', 'INDEL', 2, -1, -1, 3, 66, -1, -1] == ms.loc[
        'reference:3897:5:G'].values.tolist()

    # Mutations include unknown spdi
    ms = gds.mutations_summary('genome', id_type='spdi', ignore_annotations=True, include_unknown_features=True,
                               include_unknown_samples=False)
    assert 112 + 440 == len(ms)
    assert 'Mutation' == ms.index.name
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type',
            'Count', 'Unknown Count', 'Present and Unknown Count', 'Total',
            'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(ms.columns)

    ## Convert percent to int to make it easier to compare in assert statements
    ms['Percent'] = ms['Percent'].astype(int)
    ms['Unknown Count'] = ms['Unknown Count'].fillna(-1).astype(int)
    ms['Present and Unknown Count'] = ms['Present and Unknown Count'].fillna(-1).astype(int)
    ms['Unknown Percent'] = ms['Unknown Percent'].fillna(-1).astype(int)
    ms['Present and Unknown Percent'] = ms['Present and Unknown Percent'].fillna(-1).astype(int)

    assert ['reference', 839, 1, 'G', 'SNP', 2, -1, -1, 3, 66, -1, -1] == ms.loc['reference:839:1:G'].values.tolist()
    assert ['reference', 866, 9, 'G', 'INDEL', 1, -1, -1, 3, 33, -1, -1] == ms.loc['reference:866:9:G'].values.tolist()
    assert ['reference', 1048, 1, 'G', 'SNP', 1, -1, -1, 3, 33, -1, -1] == ms.loc['reference:1048:1:G'].values.tolist()
    assert ['reference', 3897, 5, 'G', 'INDEL', 2, -1, -1, 3, 66, -1, -1] == ms.loc[
        'reference:3897:5:G'].values.tolist()
    assert ['reference', 89, 1, '?', 'UNKNOWN_MISSING', 3, -1, -1, 3, 100, -1, -1] == ms.loc[
        'reference:89:1:?'].values.tolist()
    assert ['reference', 5100, 1, '?', 'UNKNOWN_MISSING', 2, -1, -1, 3, 66, -1, -1] == ms.loc[
        'reference:5100:1:?'].values.tolist()
    assert ['reference', 649, 1, '?', 'UNKNOWN_MISSING', 1, -1, -1, 3, 33, -1, -1] == ms.loc[
        'reference:649:1:?'].values.tolist()

    # Mutations spdi ref
    ms = gds.mutations_summary('genome', id_type='spdi_ref', ignore_annotations=True,
                               include_unknown_samples=False)
    assert 112 == len(ms)
    assert 'Mutation' == ms.index.name
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type',
            'Count', 'Unknown Count', 'Present and Unknown Count', 'Total',
            'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(ms.columns)

    ## Convert percent to int to make it easier to compare in assert statements
    ms['Percent'] = ms['Percent'].astype(int)
    ms['Unknown Count'] = ms['Unknown Count'].fillna(-1).astype(int)
    ms['Present and Unknown Count'] = ms['Present and Unknown Count'].fillna(-1).astype(int)
    ms['Unknown Percent'] = ms['Unknown Percent'].fillna(-1).astype(int)
    ms['Present and Unknown Percent'] = ms['Present and Unknown Percent'].fillna(-1).astype(int)

    assert ['reference', 839, 'C', 'G', 'SNP', 2, -1, -1, 3, 66, -1, -1] == ms.loc['reference:839:C:G'].values.tolist()
    assert ['reference', 866, 'GCCAGATCC', 'G', 'INDEL', 1, -1, -1, 3, 33, -1, -1] == ms.loc[
        'reference:866:GCCAGATCC:G'].values.tolist()
    assert ['reference', 1048, 'C', 'G', 'SNP', 1, -1, -1, 3, 33, -1, -1] == ms.loc[
        'reference:1048:C:G'].values.tolist()
    assert ['reference', 3897, 'GCGCA', 'G', 'INDEL', 2, -1, -1, 3, 66, -1, -1] == ms.loc[
        'reference:3897:GCGCA:G'].values.tolist()

    # Mutations spdi ref include unknowns
    ms = gds.mutations_summary('genome', id_type='spdi_ref', ignore_annotations=True, include_unknown_features=True,
                               include_unknown_samples=False)
    assert 112 + 440 == len(ms)
    assert 'Mutation' == ms.index.name
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type',
            'Count', 'Unknown Count', 'Present and Unknown Count', 'Total',
            'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(ms.columns)

    ## Convert percent to int to make it easier to compare in assert statements
    ms['Percent'] = ms['Percent'].astype(int)
    ms['Unknown Count'] = ms['Unknown Count'].fillna(-1).astype(int)
    ms['Present and Unknown Count'] = ms['Present and Unknown Count'].fillna(-1).astype(int)
    ms['Unknown Percent'] = ms['Unknown Percent'].fillna(-1).astype(int)
    ms['Present and Unknown Percent'] = ms['Present and Unknown Percent'].fillna(-1).astype(int)

    assert ['reference', 839, 'C', 'G', 'SNP', 2, -1, -1, 3, 66, -1, -1] == ms.loc['reference:839:C:G'].values.tolist()
    assert ['reference', 866, 'GCCAGATCC', 'G', 'INDEL', 1, -1, -1, 3, 33, -1, -1] == ms.loc[
        'reference:866:GCCAGATCC:G'].values.tolist()
    assert ['reference', 1048, 'C', 'G', 'SNP', 1, -1, -1, 3, 33, -1, -1] == ms.loc[
        'reference:1048:C:G'].values.tolist()
    assert ['reference', 3897, 'GCGCA', 'G', 'INDEL', 2, -1, -1, 3, 66, -1, -1] == ms.loc[
        'reference:3897:GCGCA:G'].values.tolist()
    assert ['reference', 89, 'A', '?', 'UNKNOWN_MISSING', 3, -1, -1, 3, 100, -1, -1] == ms.loc[
        'reference:89:A:?'].values.tolist()
    assert ['reference', 5100, 'T', '?', 'UNKNOWN_MISSING', 2, -1, -1, 3, 66, -1, -1] == ms.loc[
        'reference:5100:T:?'].values.tolist()
    assert ['reference', 649, 'T', '?', 'UNKNOWN_MISSING', 1, -1, -1, 3, 33, -1, -1] == ms.loc[
        'reference:649:T:?'].values.tolist()

    # Mutations include annotations (which should all be empty)
    ms = gds.mutations_summary('genome', id_type='spdi', ignore_annotations=False,
                               include_unknown_samples=False)
    assert 112 == len(ms)
    assert 'Mutation' == ms.index.name
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type',
            'Count', 'Unknown Count', 'Present and Unknown Count', 'Total',
            'Percent', 'Unknown Percent', 'Present and Unknown Percent',
            'Annotation', 'Annotation_Impact',
            'Gene_Name', 'Gene_ID', 'Feature_Type', 'Transcript_BioType',
            'HGVS.c', 'HGVS.p', 'ID_HGVS.c', 'ID_HGVS.p', 'ID_HGVS_GN.c', 'ID_HGVS_GN.p'] == list(ms.columns)

    ## Convert percent to int to make it easier to compare in assert statements
    ms['Percent'] = ms['Percent'].astype(int)
    ms['Unknown Count'] = ms['Unknown Count'].fillna(-1).astype(int)
    ms['Present and Unknown Count'] = ms['Present and Unknown Count'].fillna(-1).astype(int)
    ms['Unknown Percent'] = ms['Unknown Percent'].fillna(-1).astype(int)
    ms['Present and Unknown Percent'] = ms['Present and Unknown Percent'].fillna(-1).astype(int)

    assert ['reference', 839, 1, 'G', 'SNP', 2, -1, -1, 3, 66, -1, -1] + ['NA'] * 12 == ms.loc[
        'reference:839:1:G'].fillna(
        'NA').values.tolist()

    # Test case of directly calling features_summary
    ms = gds.features_summary(kind='mutations', scope='genome', id_type='spdi', ignore_annotations=True,
                              include_unknown_samples=False)
    assert 112 == len(ms)
    assert 'Mutation' == ms.index.name
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type',
            'Count', 'Unknown Count', 'Present and Unknown Count', 'Total',
            'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(ms.columns)

    ## Convert percent to int to make it easier to compare in assert statements
    ms['Percent'] = ms['Percent'].astype(int)
    ms['Unknown Count'] = ms['Unknown Count'].fillna(-1).astype(int)
    ms['Present and Unknown Count'] = ms['Present and Unknown Count'].fillna(-1).astype(int)
    ms['Unknown Percent'] = ms['Unknown Percent'].fillna(-1).astype(int)
    ms['Present and Unknown Percent'] = ms['Present and Unknown Percent'].fillna(-1).astype(int)

    assert ['reference', 839, 1, 'G', 'SNP', 2, -1, -1, 3, 66, -1, -1] == ms.loc['reference:839:1:G'].values.tolist()
    assert ['reference', 866, 9, 'G', 'INDEL', 1, -1, -1, 3, 33, -1, -1] == ms.loc['reference:866:9:G'].values.tolist()
    assert ['reference', 1048, 1, 'G', 'SNP', 1, -1, -1, 3, 33, -1, -1] == ms.loc['reference:1048:1:G'].values.tolist()
    assert ['reference', 3897, 5, 'G', 'INDEL', 2, -1, -1, 3, 66, -1, -1] == ms.loc[
        'reference:3897:5:G'].values.tolist()

    # Test case of only including unknowns
    ms = gds.mutations_summary('genome', id_type='spdi', include_present_features=False, include_unknown_features=True,
                               include_unknown_samples=False)
    assert 440 == len(ms)
    assert 'reference:649:1:?' in set(ms.index.tolist())
    assert 'reference:839:1:G' not in (ms.index.tolist())

    # Test case of no present or unknowns
    ms = gds.mutations_summary('genome', id_type='spdi', include_present_features=False, include_unknown_features=False,
                               include_unknown_samples=False)
    assert 0 == len(ms)


def test_summaries_mlst_data(loaded_database_genomic_data_store: GenomicsDataIndex):
    gds = loaded_database_genomic_data_store

    # MLST summaries for lmonocytogenes
    summary_df = gds.features_summary(kind='mlst', scope='lmonocytogenes')
    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    summary_df['Unknown Count'] = summary_df['Unknown Count'].fillna(-1).astype(int)
    summary_df['Present and Unknown Count'] = summary_df['Present and Unknown Count'].fillna(-1).astype(int)
    summary_df['Unknown Percent'] = summary_df['Unknown Percent'].fillna(-1).astype(int)
    summary_df['Present and Unknown Percent'] = summary_df['Present and Unknown Percent'].fillna(-1).astype(int)
    assert 10 == len(summary_df)
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Unknown Count', 'Present and Unknown Count',
            'Total', 'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'abcZ', '1', 5, 0, 5, 5, 100, 0, 100] == summary_df.loc[
        'mlst:lmonocytogenes:abcZ:1'].tolist()
    assert ['lmonocytogenes', 'bglA', '51', 3, 0, 3, 5, 60, 0, 60] == summary_df.loc[
        'mlst:lmonocytogenes:bglA:51'].tolist()
    assert ['lmonocytogenes', 'lhkA', '4', 1, 0, 1, 5, 20, 0, 20] == summary_df.loc[
        'mlst:lmonocytogenes:lhkA:4'].tolist()
    assert ['lmonocytogenes', 'lhkA', '5', 4, 0, 4, 5, 80, 0, 80] == summary_df.loc[
        'mlst:lmonocytogenes:lhkA:5'].tolist()
    assert ['lmonocytogenes', 'ldh', '5', 4, 1, 5, 5, 80, 20, 100] == summary_df.loc[
        'mlst:lmonocytogenes:ldh:5'].tolist()

    # MLST summaries for lmonocytogenes include unknown
    summary_df = gds.features_summary(kind='mlst', scope='lmonocytogenes', include_unknown_features=True)
    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    summary_df['Unknown Count'] = summary_df['Unknown Count'].fillna(-1).astype(int)
    summary_df['Present and Unknown Count'] = summary_df['Present and Unknown Count'].fillna(-1).astype(int)
    summary_df['Unknown Percent'] = summary_df['Unknown Percent'].fillna(-1).astype(int)
    summary_df['Present and Unknown Percent'] = summary_df['Present and Unknown Percent'].fillna(-1).astype(int)
    assert 11 == len(summary_df)
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Unknown Count', 'Present and Unknown Count',
            'Total', 'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'abcZ', '1', 5, 0, 5, 5, 100, 0, 100] == summary_df.loc[
        'mlst:lmonocytogenes:abcZ:1'].tolist()
    assert ['lmonocytogenes', 'bglA', '51', 3, 0, 3, 5, 60, 0, 60] == summary_df.loc[
        'mlst:lmonocytogenes:bglA:51'].tolist()
    assert ['lmonocytogenes', 'lhkA', '4', 1, 0, 1, 5, 20, 0, 20] == summary_df.loc[
        'mlst:lmonocytogenes:lhkA:4'].tolist()
    assert ['lmonocytogenes', 'lhkA', '5', 4, 0, 4, 5, 80, 0, 80] == summary_df.loc[
        'mlst:lmonocytogenes:lhkA:5'].tolist()
    assert ['lmonocytogenes', 'ldh', '5', 4, 1, 5, 5, 80, 20, 100] == summary_df.loc[
        'mlst:lmonocytogenes:ldh:5'].tolist()
    assert ['lmonocytogenes', 'ldh', '?', 1, -1, -1, 5, 20, -1, -1] == summary_df.loc[
        'mlst:lmonocytogenes:ldh:?'].tolist()

    # MLST summaries for lmonocytogenes with specific locus id
    summary_df = gds.features_summary(kind='mlst', scope='lmonocytogenes', locus='bglA')
    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    summary_df['Unknown Count'] = summary_df['Unknown Count'].fillna(-1).astype(int)
    summary_df['Present and Unknown Count'] = summary_df['Present and Unknown Count'].fillna(-1).astype(int)
    summary_df['Unknown Percent'] = summary_df['Unknown Percent'].fillna(-1).astype(int)
    summary_df['Present and Unknown Percent'] = summary_df['Present and Unknown Percent'].fillna(-1).astype(int)
    assert 2 == len(summary_df)
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Unknown Count', 'Present and Unknown Count',
            'Total', 'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'bglA', '51', 3, 0, 3, 5, 60, 0, 60] == summary_df.loc[
        'mlst:lmonocytogenes:bglA:51'].tolist()
    assert ['lmonocytogenes', 'bglA', '52', 2, 0, 2, 5, 40, 0, 40] == summary_df.loc[
        'mlst:lmonocytogenes:bglA:52'].tolist()

    # MLST summaries for lmonocytogenes include unknown and not present
    summary_df = gds.features_summary(kind='mlst', scope='lmonocytogenes', include_present_features=False,
                                      include_unknown_features=True)
    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    summary_df['Unknown Count'] = summary_df['Unknown Count'].fillna(-1).astype(int)
    summary_df['Present and Unknown Count'] = summary_df['Present and Unknown Count'].fillna(-1).astype(int)
    summary_df['Unknown Percent'] = summary_df['Unknown Percent'].fillna(-1).astype(int)
    summary_df['Present and Unknown Percent'] = summary_df['Present and Unknown Percent'].fillna(-1).astype(int)
    assert 1 == len(summary_df)
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Unknown Count', 'Present and Unknown Count',
            'Total', 'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(summary_df.columns)
    assert ['lmonocytogenes', 'ldh', '?', 1, -1, -1, 5, 20, -1, -1] == summary_df.loc[
        'mlst:lmonocytogenes:ldh:?'].tolist()

    # MLST summaries for lmonocytogenes not include present or unknown
    summary_df = gds.features_summary(kind='mlst', scope='lmonocytogenes', include_present_features=False,
                                      include_unknown_features=False)
    summary_df['Percent'] = summary_df['Percent'].astype(int)  # Convert to int for easier comparison
    summary_df['Unknown Count'] = summary_df['Unknown Count'].fillna(-1).astype(int)
    summary_df['Present and Unknown Count'] = summary_df['Present and Unknown Count'].fillna(-1).astype(int)
    summary_df['Unknown Percent'] = summary_df['Unknown Percent'].fillna(-1).astype(int)
    summary_df['Present and Unknown Percent'] = summary_df['Present and Unknown Percent'].fillna(-1).astype(int)
    assert 0 == len(summary_df)
    assert 'MLST Feature' == summary_df.index.name
    assert ['Scheme', 'Locus', 'Allele', 'Count', 'Unknown Count', 'Present and Unknown Count',
            'Total', 'Percent', 'Unknown Percent', 'Present and Unknown Percent'] == list(summary_df.columns)

    # Summaries using 'mlst_summery()'
    summary_df = gds.mlst_summary(scheme_name='lmonocytogenes')
    assert 10 == len(summary_df)


def test_summaries_variant_annotations(
        loaded_database_genomic_data_store_annotations_include_unknown: GenomicsDataIndex):
    gds = loaded_database_genomic_data_store_annotations_include_unknown

    # Samples
    assert 3 == gds.count_samples()
    assert {'SH10-014', 'SH14-014', 'SH14-001'} == set(gds.sample_names())

    # References
    assert 1 == gds.count_references()
    assert ['NC_011083'] == gds.reference_names()

    # Mutations
    assert 177 == gds.count_mutations('NC_011083')

    # spdi
    ms = gds.mutations_summary('NC_011083', id_type='spdi', ignore_annotations=False)
    assert 177 == len(ms)
    assert 'Mutation' == ms.index.name
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type',
            'Count', 'Unknown Count', 'Present and Unknown Count', 'Total',
            'Percent', 'Unknown Percent', 'Present and Unknown Percent',
            'Annotation', 'Annotation_Impact',
            'Gene_Name', 'Gene_ID', 'Feature_Type', 'Transcript_BioType',
            'HGVS.c', 'HGVS.p', 'ID_HGVS.c', 'ID_HGVS.p', 'ID_HGVS_GN.c', 'ID_HGVS_GN.p'] == list(ms.columns)

    ## Convert percent to int to make it easier to compare in assert statements
    ms['Percent'] = ms['Percent'].astype(int)
    ms['Unknown Count'] = ms['Unknown Count'].fillna(-1).astype(int)
    ms['Present and Unknown Count'] = ms['Present and Unknown Count'].fillna(-1).astype(int)
    ms['Unknown Percent'] = ms['Unknown Percent'].fillna(-1).astype(int)
    ms['Present and Unknown Percent'] = ms['Present and Unknown Percent'].fillna(-1).astype(int)

    ## missense variant
    assert ['NC_011083', 140658, 1, 'A', 'SNP', 3, 0, 3, 3, 100, 0, 100,
            'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
            'c.497C>A', 'p.Ala166Glu',
            'hgvs:NC_011083:SEHA_RS01180:c.497C>A', 'hgvs:NC_011083:SEHA_RS01180:p.Ala166Glu',
            'hgvs_gn:NC_011083:murF:c.497C>A', 'hgvs_gn:NC_011083:murF:p.Ala166Glu'] == list(
        ms.loc['NC_011083:140658:1:A'])

    ## inframe deletion
    assert ['NC_011083', 4465400, len('GGCCGAA'), 'G', 'INDEL', 3, 0, 3, 3, 100, 0, 100,
            'conservative_inframe_deletion', 'MODERATE', 'tyrB', 'SEHA_RS22180', 'transcript', 'protein_coding',
            'c.157_162delGAAGCC', 'p.Glu53_Ala54del',
            'hgvs:NC_011083:SEHA_RS22180:c.157_162delGAAGCC', 'hgvs:NC_011083:SEHA_RS22180:p.Glu53_Ala54del',
            'hgvs_gn:NC_011083:tyrB:c.157_162delGAAGCC', 'hgvs_gn:NC_011083:tyrB:p.Glu53_Ala54del'] == list(
        ms.loc['NC_011083:4465400:7:G'])

    ## Intergenic variant (with some NA values in fields)
    assert ['NC_011083', 4555461, len('T'), 'TC', 'INDEL', 1, 2, 3, 3, 33, 66, 100,
            'intergenic_region', 'MODIFIER', 'SEHA_RS22510-SEHA_RS26685', 'SEHA_RS22510-SEHA_RS26685',
            'intergenic_region', 'NA',
            'n.4555461_4555462insC', 'NA',
            'hgvs:NC_011083:n.4555461_4555462insC', 'NA',
            'hgvs_gn:NC_011083:n.4555461_4555462insC', 'NA'] == list(ms.loc['NC_011083:4555461:1:TC'].fillna('NA'))

    ## MNP variant (2/3, 1/3)
    assert ['NC_011083', 3535698, len('GCC'), 'CAT', 'MNP', 2, 1, 3, 3, 66, 33, 100,
            'missense_variant', 'MODERATE', 'oadA', 'SEHA_RS17780',
            'transcript', 'protein_coding',
            'c.544_546delGGCinsATG', 'p.Gly182Met',
            'hgvs:NC_011083:SEHA_RS17780:c.544_546delGGCinsATG', 'hgvs:NC_011083:SEHA_RS17780:p.Gly182Met',
            'hgvs_gn:NC_011083:oadA:c.544_546delGGCinsATG', 'hgvs_gn:NC_011083:oadA:p.Gly182Met'] == list(
        ms.loc['NC_011083:3535698:3:CAT'])

    ## Long MNP variant (2/3, 1/3)
    assert ['NC_011083', 3535143, len('AATGCCTGCC'), 'TATCCCGGCG', 'MNP', 2, 1, 3, 3, 66, 33, 100,
            'synonymous_variant', 'LOW', 'oadA', 'SEHA_RS17780', 'transcript', 'protein_coding',
            'c.1092_1101delGGCAGGCATTinsCGCCGGGATA', 'p.368',
            'hgvs:NC_011083:SEHA_RS17780:c.1092_1101delGGCAGGCATTinsCGCCGGGATA',
            'hgvs:NC_011083:SEHA_RS17780:p.368',
            'hgvs_gn:NC_011083:oadA:c.1092_1101delGGCAGGCATTinsCGCCGGGATA',
            'hgvs_gn:NC_011083:oadA:p.368'] == list(
        ms.loc[f'NC_011083:3535143:{len("AATGCCTGCC")}:TATCCCGGCG'])

    ## synonymous variant, no unknown (1/3, 0/3)
    assert ['NC_011083', 508378, len('C'), 'T', 'SNP', 1, 0, 1, 3, 33, 0, 33,
            'synonymous_variant', 'LOW', 'tgt', 'SEHA_RS02965', 'transcript', 'protein_coding',
            'c.423C>T', 'p.Ile141Ile',
            'hgvs:NC_011083:SEHA_RS02965:c.423C>T', 'hgvs:NC_011083:SEHA_RS02965:p.Ile141Ile',
            'hgvs_gn:NC_011083:tgt:c.423C>T', 'hgvs_gn:NC_011083:tgt:p.Ile141Ile'] == list(
        ms.loc['NC_011083:508378:1:T'])

    # deletion where there is an overlap with present and unknown, (present 3/3, unknown 2/3,
    # overlap 1676762 and 1676763)
    assert ['NC_011083', 1676762, len('CA'), 'C', 'INDEL', 1, 2, 3, 3, 33, 66, 100,
            'intergenic_region', 'MODIFIER', 'SEHA_RS26130-SEHA_RS08880', 'SEHA_RS26130-SEHA_RS08880',
            'intergenic_region', '<NA>',
            'n.1676763delA', '<NA>',
            'hgvs:NC_011083:n.1676763delA', '<NA>',
            'hgvs_gn:NC_011083:n.1676763delA', '<NA>'] == list(
        ms.loc['NC_011083:1676762:2:C'].fillna('<NA>'))

    # spdi_ref
    ms = gds.mutations_summary('NC_011083', id_type='spdi_ref')
    assert 177 == len(ms)
    assert 'Mutation' == ms.index.name
    assert ['Sequence', 'Position', 'Deletion', 'Insertion', 'Type',
            'Count', 'Unknown Count', 'Present and Unknown Count', 'Total',
            'Percent', 'Unknown Percent', 'Present and Unknown Percent',
            'Annotation', 'Annotation_Impact',
            'Gene_Name', 'Gene_ID', 'Feature_Type', 'Transcript_BioType',
            'HGVS.c', 'HGVS.p', 'ID_HGVS.c', 'ID_HGVS.p', 'ID_HGVS_GN.c', 'ID_HGVS_GN.p'] == list(ms.columns)

    ## Convert percent to int to make it easier to compare in assert statements
    ms['Percent'] = ms['Percent'].astype(int)
    ms['Unknown Count'] = ms['Unknown Count'].fillna(-1).astype(int)
    ms['Present and Unknown Count'] = ms['Present and Unknown Count'].fillna(-1).astype(int)
    ms['Unknown Percent'] = ms['Unknown Percent'].fillna(-1).astype(int)
    ms['Present and Unknown Percent'] = ms['Present and Unknown Percent'].fillna(-1).astype(int)

    ## missense variant
    assert ['NC_011083', 140658, 'C', 'A', 'SNP', 3, 0, 3, 3, 100, 0, 100,
            'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
            'c.497C>A', 'p.Ala166Glu',
            'hgvs:NC_011083:SEHA_RS01180:c.497C>A', 'hgvs:NC_011083:SEHA_RS01180:p.Ala166Glu',
            'hgvs_gn:NC_011083:murF:c.497C>A', 'hgvs_gn:NC_011083:murF:p.Ala166Glu'] == list(
        ms.loc['NC_011083:140658:C:A'])

    ## inframe deletion
    assert ['NC_011083', 4465400, 'GGCCGAA', 'G', 'INDEL', 3, 0, 3, 3, 100, 0, 100,
            'conservative_inframe_deletion', 'MODERATE', 'tyrB', 'SEHA_RS22180', 'transcript', 'protein_coding',
            'c.157_162delGAAGCC', 'p.Glu53_Ala54del',
            'hgvs:NC_011083:SEHA_RS22180:c.157_162delGAAGCC', 'hgvs:NC_011083:SEHA_RS22180:p.Glu53_Ala54del',
            'hgvs_gn:NC_011083:tyrB:c.157_162delGAAGCC', 'hgvs_gn:NC_011083:tyrB:p.Glu53_Ala54del'] == list(
        ms.loc['NC_011083:4465400:GGCCGAA:G'])

    ## Intergenic variant (with some NA values in fields)
    assert ['NC_011083', 4555461, 'T', 'TC', 'INDEL', 1, 2, 3, 3, 33, 66, 100,
            'intergenic_region', 'MODIFIER', 'SEHA_RS22510-SEHA_RS26685', 'SEHA_RS22510-SEHA_RS26685',
            'intergenic_region', 'NA',
            'n.4555461_4555462insC', 'NA',
            'hgvs:NC_011083:n.4555461_4555462insC', 'NA',
            'hgvs_gn:NC_011083:n.4555461_4555462insC', 'NA'] == list(ms.loc['NC_011083:4555461:T:TC'].fillna('NA'))

    ## MNP variant (2/3, 1/3)
    assert ['NC_011083', 3535698, 'GCC', 'CAT', 'MNP', 2, 1, 3, 3, 66, 33, 100,
            'missense_variant', 'MODERATE', 'oadA', 'SEHA_RS17780',
            'transcript', 'protein_coding',
            'c.544_546delGGCinsATG', 'p.Gly182Met',
            'hgvs:NC_011083:SEHA_RS17780:c.544_546delGGCinsATG', 'hgvs:NC_011083:SEHA_RS17780:p.Gly182Met',
            'hgvs_gn:NC_011083:oadA:c.544_546delGGCinsATG', 'hgvs_gn:NC_011083:oadA:p.Gly182Met'] == list(
        ms.loc['NC_011083:3535698:GCC:CAT'])

    ## Long MNP variant (2/3, 1/3)
    assert ['NC_011083', 3535143, 'AATGCCTGCC', 'TATCCCGGCG', 'MNP', 2, 1, 3, 3, 66, 33, 100,
            'synonymous_variant', 'LOW', 'oadA', 'SEHA_RS17780', 'transcript', 'protein_coding',
            'c.1092_1101delGGCAGGCATTinsCGCCGGGATA', 'p.368',
            'hgvs:NC_011083:SEHA_RS17780:c.1092_1101delGGCAGGCATTinsCGCCGGGATA',
            'hgvs:NC_011083:SEHA_RS17780:p.368',
            'hgvs_gn:NC_011083:oadA:c.1092_1101delGGCAGGCATTinsCGCCGGGATA',
            'hgvs_gn:NC_011083:oadA:p.368'] == list(
        ms.loc['NC_011083:3535143:AATGCCTGCC:TATCCCGGCG'])

    ## synonymous variant, no unknown (1/3, 0/3)
    assert ['NC_011083', 508378, 'C', 'T', 'SNP', 1, 0, 1, 3, 33, 0, 33,
            'synonymous_variant', 'LOW', 'tgt', 'SEHA_RS02965', 'transcript', 'protein_coding',
            'c.423C>T', 'p.Ile141Ile',
            'hgvs:NC_011083:SEHA_RS02965:c.423C>T', 'hgvs:NC_011083:SEHA_RS02965:p.Ile141Ile',
            'hgvs_gn:NC_011083:tgt:c.423C>T', 'hgvs_gn:NC_011083:tgt:p.Ile141Ile'] == list(
        ms.loc['NC_011083:508378:C:T'])

    # deletion where there is an overlap with present and unknown, (present 3/3, unknown 2/3,
    # overlap 1676762 and 1676763)
    assert ['NC_011083', 1676762, 'CA', 'C', 'INDEL', 1, 2, 3, 3, 33, 66, 100,
            'intergenic_region', 'MODIFIER', 'SEHA_RS26130-SEHA_RS08880', 'SEHA_RS26130-SEHA_RS08880',
            'intergenic_region', '<NA>',
            'n.1676763delA', '<NA>',
            'hgvs:NC_011083:n.1676763delA', '<NA>',
            'hgvs_gn:NC_011083:n.1676763delA', '<NA>'] == list(
        ms.loc['NC_011083:1676762:CA:C'].fillna('<NA>'))


def test_connect_to_project_from_dir():
    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)
        project_dir = tmp_file / 'project'
        Project.initialize_project(project_dir)

        ds = GenomicsDataIndex.connect(project_dir=project_dir)
        assert ds is not None
        assert ds.connection.reference_service is not None
        assert ds.connection.filesystem_storage.variation_dir.parent == project_dir / '.gdi-data'


def test_connect_to_project_from_project():
    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)
        project_dir = tmp_file / 'project'
        project = Project.initialize_project(project_dir)

        ds = GenomicsDataIndex.connect(project=project)
        assert ds is not None
        assert ds.connection.reference_service is not None
        assert ds.connection.filesystem_storage.variation_dir.parent == project_dir / '.gdi-data'


def test_connect_and_tree_after_moving_project(loaded_data_store_from_project_dir):
    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)
        project_dir = tmp_file / 'project'
        Project.initialize_project(project_dir)
        ds = GenomicsDataIndex.connect(project_dir=project_dir)
        database_connection = ds.connection

        # Load Nucleotide variation
        database_connection.reference_service.add_reference_genome(reference_file)
        snippy_tmp_dir = Path(tempfile.mkdtemp())
        data_package = NucleotideSampleDataPackage.create_from_snippy(sample_dirs,
                                                                      SerialSampleFilesProcessor(snippy_tmp_dir))
        database_connection.variation_service.insert(data_package, feature_scope_name='genome')

        # I should be able to build tree initially
        tq1 = ds.samples_query().build_tree(kind='mutation', scope='genome')
        assert tq1.tree is not None

        # Move project
        project_dir_2 = tmp_file / 'project2'
        shutil.move(project_dir, project_dir_2)

        ds2 = GenomicsDataIndex.connect(project_dir=project_dir_2)

        # I should still be able to build tree even after project has been moved
        tq2 = ds2.samples_query().build_tree(kind='mutation', scope='genome')
        assert tq2.tree is not None
