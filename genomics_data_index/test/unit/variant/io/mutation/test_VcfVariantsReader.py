import pandas as pd

from genomics_data_index.storage.io.mutation.VcfVariantsReader import VcfVariantsReader
from genomics_data_index.storage.MaskedGenomicRegions import MaskedGenomicRegions
from pybedtools import BedTool


def test_mask_to_features():
    reader = VcfVariantsReader({})
    mask = MaskedGenomicRegions(BedTool([('ref', 10, 15)]))
    expected_number_features = 5

    features_df = reader.mask_to_features(genomic_mask=mask)
    features_df = features_df.sort_values('POS')

    assert ['CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'VARIANT_ID'] == features_df.columns.tolist()
    assert expected_number_features == len(features_df)
    assert ['ref'] * expected_number_features == features_df['CHROM'].tolist()
    assert list(range(10 + 1, 15 + 1)) == features_df['POS'].tolist()
    assert [1] * expected_number_features == features_df['REF'].tolist()
    assert ['?'] * expected_number_features == features_df['ALT'].tolist()
    assert ['UNKNOWN_MISSING'] * expected_number_features == features_df['TYPE'].tolist()
    assert [f'ref:{p}:1:?' for p in range(10 + 1, 15 + 1)] == features_df['VARIANT_ID'].tolist()


def test_mask_to_features_multi_range():
    reader = VcfVariantsReader({})
    mask = MaskedGenomicRegions(BedTool([('ref', 10, 15), ('ref', 20, 22)]))
    expected_number_features = 7

    features_df = reader.mask_to_features(genomic_mask=mask)
    features_df = features_df.sort_values('POS')

    assert ['CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'VARIANT_ID'] == features_df.columns.tolist()
    assert expected_number_features == len(features_df)
    assert ['ref'] * expected_number_features == features_df['CHROM'].tolist()
    assert [11, 12, 13, 14, 15, 21, 22] == features_df['POS'].tolist()
    assert [1] * expected_number_features == features_df['REF'].tolist()
    assert ['?'] * expected_number_features == features_df['ALT'].tolist()
    assert ['UNKNOWN_MISSING'] * expected_number_features == features_df['TYPE'].tolist()
    assert [f'ref:{p}:1:?' for p in [11, 12, 13, 14, 15, 21, 22]] == features_df['VARIANT_ID'].tolist()


def test_mask_to_features_multi_sequence():
    reader = VcfVariantsReader({})
    mask = MaskedGenomicRegions(BedTool([('ref', 10, 15), ('ref2', 12, 17)]))
    expected_nf = 5
    expected_nf2 = 5
    expected_t = 10

    features_df = reader.mask_to_features(genomic_mask=mask)
    features_df = features_df.sort_values(['CHROM', 'POS'])

    assert ['CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'VARIANT_ID'] == features_df.columns.tolist()
    assert expected_t == len(features_df)
    assert ['ref'] * expected_nf + ['ref2'] * expected_nf2 == features_df['CHROM'].tolist()
    assert [11, 12, 13, 14, 15, 13, 14, 15, 16, 17] == features_df['POS'].tolist()
    assert [1] * expected_t == features_df['REF'].tolist()
    assert ['?'] * expected_t == features_df['ALT'].tolist()
    assert ['UNKNOWN_MISSING'] * expected_t == features_df['TYPE'].tolist()
    assert [f'ref:{p}:1:?' for p in [11, 12, 13, 14, 15]] + [
        f'ref2:{p}:1:?' for p in [13, 14, 15, 16, 17]] == features_df['VARIANT_ID'].tolist()


def test_combine_vcf_mask():
    num_annotations = 9

    data_vcf = [
        ['SampleA', 'ref', 10, 'A', 'T', 'SNP', 'file', 'ref:10:A:T'] + [pd.NA] * num_annotations,
    ]
    data_mask = [
        ['SampleA', 'ref', 1, 1, '?', 'UNKNOWN_MISSING', 'file', 'ref:1:1:?']
    ]
    columns = ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
               'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
               'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p']

    reader = VcfVariantsReader({})

    vcf_df = pd.DataFrame(data_vcf, columns=columns)
    mask_df = pd.DataFrame(data_mask, columns=['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID'])

    combined_df = reader.combine_vcf_mask(vcf_frame=vcf_df, mask_frame=mask_df)
    combined_df = combined_df.sort_values('POS').fillna('NA')

    assert combined_df.columns.tolist() == columns
    assert 2 == len(combined_df)
    assert ['ref', 'ref'] == combined_df['CHROM'].tolist()
    assert [1, 10] == combined_df['POS'].tolist()
    assert [1, 'A'] == combined_df['REF'].tolist()
    assert ['?', 'T'] == combined_df['ALT'].tolist()
    assert ['UNKNOWN_MISSING', 'SNP'] == combined_df['TYPE'].tolist()
    assert ['ref:1:1:?', 'ref:10:A:T'] == combined_df['VARIANT_ID'].tolist()
    assert ['NA', 'NA'] == combined_df['ANN.Annotation'].tolist()


def test_combine_vcf_mask_overlap_feature():
    num_annotations = 9

    data_vcf = [
        ['SampleA', 'ref', 10, 'A', 'T', 'SNP', 'file', 'ref:10:A:T'] + [pd.NA] * num_annotations,
    ]
    data_mask = [
        ['SampleA', 'ref', 10, 1, '?', 'UNKNOWN_MISSING', 'file', 'ref:10:1:?']
    ]
    columns = ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
               'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
               'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p']

    reader = VcfVariantsReader({})

    vcf_df = pd.DataFrame(data_vcf, columns=columns)
    mask_df = pd.DataFrame(data_mask, columns=['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID'])

    combined_df = reader.combine_vcf_mask(vcf_frame=vcf_df, mask_frame=mask_df)
    combined_df = combined_df.sort_values('POS').fillna('NA')

    assert combined_df.columns.tolist() == columns
    assert 1 == len(combined_df)
    assert ['ref'] == combined_df['CHROM'].tolist()
    assert [10] == combined_df['POS'].tolist()
    assert [1] == combined_df['REF'].tolist()
    assert ['?'] == combined_df['ALT'].tolist()
    assert ['UNKNOWN_MISSING'] == combined_df['TYPE'].tolist()
    assert ['ref:10:1:?'] == combined_df['VARIANT_ID'].tolist()
    assert ['NA'] == combined_df['ANN.Annotation'].tolist()


def test_combine_vcf_mask_no_mask_features():
    num_annotations = 9

    data_vcf = [
        ['SampleA', 'ref', 10, 'A', 'T', 'SNP', 'file', 'ref:10:A:T'] + [pd.NA] * num_annotations,
    ]
    data_mask = [
    ]
    columns = ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
               'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
               'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p']

    reader = VcfVariantsReader({})

    vcf_df = pd.DataFrame(data_vcf, columns=columns)
    mask_df = pd.DataFrame(data_mask, columns=['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID'])

    combined_df = reader.combine_vcf_mask(vcf_frame=vcf_df, mask_frame=mask_df)
    combined_df = combined_df.sort_values('POS').fillna('NA')

    assert combined_df.columns.tolist() == columns
    assert 1 == len(combined_df)
    assert ['ref'] == combined_df['CHROM'].tolist()
    assert [10] == combined_df['POS'].tolist()
    assert ['A'] == combined_df['REF'].tolist()
    assert ['T'] == combined_df['ALT'].tolist()
    assert ['SNP'] == combined_df['TYPE'].tolist()
    assert ['ref:10:A:T'] == combined_df['VARIANT_ID'].tolist()
    assert ['NA'] == combined_df['ANN.Annotation'].tolist()


def test_combine_vcf_mask_no_vcf_feature():
    data_vcf = [
    ]
    data_mask = [
        ['SampleA', 'ref', 10, 1, '?', 'UNKNOWN_MISSING', 'file', 'ref:10:1:?']
    ]
    columns = ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
               'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
               'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p']

    reader = VcfVariantsReader({})

    vcf_df = pd.DataFrame(data_vcf, columns=columns)
    mask_df = pd.DataFrame(data_mask, columns=['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID'])

    combined_df = reader.combine_vcf_mask(vcf_frame=vcf_df, mask_frame=mask_df)
    combined_df = combined_df.sort_values('POS').fillna('NA')

    assert combined_df.columns.tolist() == columns
    assert 1 == len(combined_df)
    assert ['ref'] == combined_df['CHROM'].tolist()
    assert [10] == combined_df['POS'].tolist()
    assert [1] == combined_df['REF'].tolist()
    assert ['?'] == combined_df['ALT'].tolist()
    assert ['UNKNOWN_MISSING'] == combined_df['TYPE'].tolist()
    assert ['ref:10:1:?'] == combined_df['VARIANT_ID'].tolist()
    assert ['NA'] == combined_df['ANN.Annotation'].tolist()


def test_combine_vcf_mask_same_position_different_sequence():
    num_annotations = 9

    data_vcf = [
        ['SampleA', 'ref', 10, 'A', 'T', 'SNP', 'file', 'ref:10:A:T'] + [pd.NA] * num_annotations,
    ]
    data_mask = [
        ['SampleA', 'ref2', 10, 1, '?', 'UNKNOWN_MISSING', 'file', 'ref2:10:1:?']
    ]
    columns = ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
               'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
               'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p']

    reader = VcfVariantsReader({})

    vcf_df = pd.DataFrame(data_vcf, columns=columns)
    mask_df = pd.DataFrame(data_mask, columns=['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID'])

    combined_df = reader.combine_vcf_mask(vcf_frame=vcf_df, mask_frame=mask_df)
    combined_df = combined_df.sort_values(['CHROM', 'POS']).fillna('NA')

    assert combined_df.columns.tolist() == columns
    assert 2 == len(combined_df)
    assert ['ref', 'ref2'] == combined_df['CHROM'].tolist()
    assert [10, 10] == combined_df['POS'].tolist()
    assert ['A', 1] == combined_df['REF'].tolist()
    assert ['T', '?'] == combined_df['ALT'].tolist()
    assert ['SNP', 'UNKNOWN_MISSING'] == combined_df['TYPE'].tolist()
    assert ['ref:10:A:T', 'ref2:10:1:?'] == combined_df['VARIANT_ID'].tolist()
    assert ['NA', 'NA'] == combined_df['ANN.Annotation'].tolist()


def test_combine_vcf_mask_multiple_features():
    num_annotations = 9

    data_vcf = [
        ['SampleA', 'ref', 10, 'A', 'T', 'SNP', 'file', 'ref:10:A:T'] + [pd.NA] * num_annotations,
        ['SampleA', 'ref', 20, 'ATT', 'T', 'INDEL', 'file', 'ref:20:ATT:T'] + [pd.NA] * num_annotations,
    ]
    data_mask = [
        ['SampleA', 'ref', 10, 1, '?', 'UNKNOWN_MISSING', 'file', 'ref:10:1:?'],
        ['SampleA', 'ref', 30, 1, '?', 'UNKNOWN_MISSING', 'file', 'ref:30:1:?']
    ]
    columns = ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
               'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
               'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p']

    reader = VcfVariantsReader({})

    vcf_df = pd.DataFrame(data_vcf, columns=columns)
    mask_df = pd.DataFrame(data_mask, columns=['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID'])

    combined_df = reader.combine_vcf_mask(vcf_frame=vcf_df, mask_frame=mask_df)
    combined_df = combined_df.sort_values('POS').fillna('NA')

    assert combined_df.columns.tolist() == columns
    assert 3 == len(combined_df)
    assert ['ref', 'ref', 'ref'] == combined_df['CHROM'].tolist()
    assert [10, 20, 30] == combined_df['POS'].tolist()
    assert [1, 'ATT', 1] == combined_df['REF'].tolist()
    assert ['?', 'T', '?'] == combined_df['ALT'].tolist()
    assert ['UNKNOWN_MISSING', 'INDEL', 'UNKNOWN_MISSING'] == combined_df['TYPE'].tolist()
    assert ['ref:10:1:?', 'ref:20:ATT:T', 'ref:30:1:?'] == combined_df['VARIANT_ID'].tolist()
    assert ['NA', 'NA', 'NA'] == combined_df['ANN.Annotation'].tolist()
