import pandas as pd

from genomics_data_index.storage.io.mutation.NucleotideSampleData import NucleotideSampleData


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

    sample_data = NucleotideSampleData('SampleA', vcf_file=None, vcf_file_index=None,
                                       mask_bed_file=None, preprocessed=True)

    vcf_df = pd.DataFrame(data_vcf, columns=columns)
    mask_df = pd.DataFrame(data_mask, columns=['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID'])

    combined_df = sample_data.combine_vcf_mask(vcf_frame=vcf_df, mask_frame=mask_df)
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

    sample_data = NucleotideSampleData('SampleA', vcf_file=None, vcf_file_index=None,
                                       mask_bed_file=None, preprocessed=True)

    vcf_df = pd.DataFrame(data_vcf, columns=columns)
    mask_df = pd.DataFrame(data_mask, columns=['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID'])

    combined_df = sample_data.combine_vcf_mask(vcf_frame=vcf_df, mask_frame=mask_df)
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

    sample_data = NucleotideSampleData('SampleA', vcf_file=None, vcf_file_index=None,
                                       mask_bed_file=None, preprocessed=True)

    vcf_df = pd.DataFrame(data_vcf, columns=columns)
    mask_df = pd.DataFrame(data_mask, columns=['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID'])

    combined_df = sample_data.combine_vcf_mask(vcf_frame=vcf_df, mask_frame=mask_df)
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

    sample_data = NucleotideSampleData('SampleA', vcf_file=None, vcf_file_index=None,
                                       mask_bed_file=None, preprocessed=True)

    vcf_df = pd.DataFrame(data_vcf, columns=columns)
    mask_df = pd.DataFrame(data_mask, columns=['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID'])

    combined_df = sample_data.combine_vcf_mask(vcf_frame=vcf_df, mask_frame=mask_df)
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

    sample_data = NucleotideSampleData('SampleA', vcf_file=None, vcf_file_index=None,
                                       mask_bed_file=None, preprocessed=True)

    vcf_df = pd.DataFrame(data_vcf, columns=columns)
    mask_df = pd.DataFrame(data_mask, columns=['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID'])

    combined_df = sample_data.combine_vcf_mask(vcf_frame=vcf_df, mask_frame=mask_df)
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

    sample_data = NucleotideSampleData('SampleA', vcf_file=None, vcf_file_index=None,
                                       mask_bed_file=None, preprocessed=True)

    vcf_df = pd.DataFrame(data_vcf, columns=columns)
    mask_df = pd.DataFrame(data_mask, columns=['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID'])

    combined_df = sample_data.combine_vcf_mask(vcf_frame=vcf_df, mask_frame=mask_df)
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
