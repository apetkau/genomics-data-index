from pathlib import Path


from genomics_data_index.storage.io.mutation.VcfVariantsReader import VcfVariantsReader
from genomics_data_index.storage.MaskedGenomicRegions import MaskedGenomicRegions
from pybedtools import BedTool


def test_mask_to_features():
    reader = VcfVariantsReader({})
    vcf_file = Path('/tmp/file.vcf')
    mask = MaskedGenomicRegions(BedTool([('ref', 10, 15)]))
    expected_number_features = 5

    features_df = reader.mask_to_features(genomic_mask=mask, vcf_file=vcf_file, sample='SampleA')
    features_df = features_df.sort_values('POS')

    assert ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID'] == features_df.columns.tolist()
    assert expected_number_features == len(features_df)
    assert ['SampleA'] * expected_number_features == features_df['SAMPLE'].tolist()
    assert ['ref'] * expected_number_features == features_df['CHROM'].tolist()
    assert list(range(10+1, 15+1)) == features_df['POS'].tolist()
    assert [1] * expected_number_features == features_df['REF'].tolist()
    assert ['?'] * expected_number_features == features_df['ALT'].tolist()
    assert ['UNKNOWN'] * expected_number_features == features_df['TYPE'].tolist()
    assert ['file.vcf'] * expected_number_features == features_df['FILE'].tolist()
    assert [f'ref:{p}:1:?' for p in range(10+1, 15+1)] == features_df['VARIANT_ID'].tolist()


def test_mask_to_features_multi_range():
    reader = VcfVariantsReader({})
    vcf_file = Path('/tmp/file.vcf')
    mask = MaskedGenomicRegions(BedTool([('ref', 10, 15), ('ref', 20, 22)]))
    expected_number_features = 7

    features_df = reader.mask_to_features(genomic_mask=mask, vcf_file=vcf_file, sample='SampleA')
    features_df = features_df.sort_values('POS')

    assert ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID'] == features_df.columns.tolist()
    assert expected_number_features == len(features_df)
    assert ['SampleA'] * expected_number_features == features_df['SAMPLE'].tolist()
    assert ['ref'] * expected_number_features == features_df['CHROM'].tolist()
    assert [11, 12, 13, 14, 15, 21, 22] == features_df['POS'].tolist()
    assert [1] * expected_number_features == features_df['REF'].tolist()
    assert ['?'] * expected_number_features == features_df['ALT'].tolist()
    assert ['UNKNOWN'] * expected_number_features == features_df['TYPE'].tolist()
    assert ['file.vcf'] * expected_number_features == features_df['FILE'].tolist()
    assert [f'ref:{p}:1:?' for p in [11, 12, 13, 14, 15, 21, 22]] == features_df['VARIANT_ID'].tolist()


def test_mask_to_features_multi_sequence():
    reader = VcfVariantsReader({})
    vcf_file = Path('/tmp/file.vcf')
    mask = MaskedGenomicRegions(BedTool([('ref', 10, 15), ('ref2', 12, 17)]))
    expected_nf = 5
    expected_nf2 = 5
    expected_t = 10

    features_df = reader.mask_to_features(genomic_mask=mask, vcf_file=vcf_file, sample='SampleA')
    features_df = features_df.sort_values(['CHROM', 'POS'])

    assert ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID'] == features_df.columns.tolist()
    assert expected_t == len(features_df)
    assert ['SampleA'] * expected_t == features_df['SAMPLE'].tolist()
    assert ['ref'] * expected_nf + ['ref2'] * expected_nf2 == features_df['CHROM'].tolist()
    assert [11, 12, 13, 14, 15, 13, 14, 15, 16, 17] == features_df['POS'].tolist()
    assert [1] * expected_t == features_df['REF'].tolist()
    assert ['?'] * expected_t == features_df['ALT'].tolist()
    assert ['UNKNOWN'] * expected_t == features_df['TYPE'].tolist()
    assert ['file.vcf'] * expected_t == features_df['FILE'].tolist()
    assert [f'ref:{p}:1:?' for p in [11, 12, 13, 14, 15]] + [
        f'ref2:{p}:1:?' for p in [13, 14, 15, 16, 17]] == features_df['VARIANT_ID'].tolist()
