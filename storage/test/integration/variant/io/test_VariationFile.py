import pytest
import tempfile
from pathlib import Path

from storage.test.integration import data_dir, regular_vcf_dir, variation_dir, reference_file, consensus_dir
from storage.test.integration import extra_snippy_dir
from storage.test.integration.variant.io import read_vcf_df
from storage.variant.MaskedGenomicRegions import MaskedGenomicRegions
from storage.variant.io.mutation.VariationFile import VariationFile
from storage.variant.util import parse_sequence_file


def test_write():
    sample_vcf = data_dir / 'SampleA' / 'snps.vcf.gz'
    with tempfile.TemporaryDirectory() as out_dir:
        out_file = Path(out_dir) / 'out.vcf.gz'

        assert not out_file.exists()
        file, index = VariationFile(sample_vcf).write(out_file)
        assert out_file.exists()
        assert file == out_file
        assert index.exists()
        assert str(file) + '.csi' == str(index)

        df = read_vcf_df(out_file)
        assert 'SNP' == df[df['POS'] == 293]['TYPE'].tolist()[0]
        assert 'INDEL' == df[df['POS'] == 302]['TYPE'].tolist()[0]
        assert 'INDEL' == df[df['POS'] == 324]['TYPE'].tolist()[0]
        assert 'INDEL' == df[df['POS'] == 374]['TYPE'].tolist()[0]
        assert 'OTHER' == df[df['POS'] == 461]['TYPE'].tolist()[0]
        assert 'SNP' == df[df['POS'] == 506]['TYPE'].tolist()[0]


def test_write_2():
    sample_vcf = data_dir / 'SampleC' / 'snps.vcf.gz'
    with tempfile.TemporaryDirectory() as out_dir:
        out_file = Path(out_dir) / 'out.vcf.gz'

        assert not out_file.exists()
        file, index = VariationFile(sample_vcf).write(out_file)
        assert out_file.exists()
        assert file == out_file
        assert index.exists()
        assert str(file) + '.csi' == str(index)

        df = read_vcf_df(out_file)
        assert 'INDEL' == df[df['POS'] == 347]['TYPE'].tolist()[0]
        assert 'SNP' == df[df['POS'] == 619]['TYPE'].tolist()[0]
        assert 'OTHER' == df[df['POS'] == 1984]['TYPE'].tolist()[0]


def test_write_missing_type_tag():
    sample_vcf = regular_vcf_dir / 'SampleA.vcf.gz'
    with tempfile.TemporaryDirectory() as out_dir:
        out_file = Path(out_dir) / 'out.vcf.gz'

        assert not out_file.exists()
        VariationFile(sample_vcf).write(out_file)
        assert out_file.exists()

        df = read_vcf_df(out_file)
        assert 'SNP' == df[df['POS'] == 293]['TYPE'].tolist()[0]
        assert 'INDEL' == df[df['POS'] == 302]['TYPE'].tolist()[0]
        assert 'INDEL' == df[df['POS'] == 324]['TYPE'].tolist()[0]
        assert 'INDEL' == df[df['POS'] == 374]['TYPE'].tolist()[0]
        assert 'OTHER' == df[df['POS'] == 461]['TYPE'].tolist()[0]
        assert 'SNP' == df[df['POS'] == 506]['TYPE'].tolist()[0]


def test_write_2_missing_type_tag():
    sample_vcf = regular_vcf_dir / 'SampleC.vcf.gz'
    with tempfile.TemporaryDirectory() as out_dir:
        out_file = Path(out_dir) / 'out.vcf.gz'

        assert not out_file.exists()
        VariationFile(sample_vcf).write(out_file)
        assert out_file.exists()

        df = read_vcf_df(out_file)
        assert 'INDEL' == df[df['POS'] == 347]['TYPE'].tolist()[0]
        assert 'SNP' == df[df['POS'] == 619]['TYPE'].tolist()[0]
        assert 'OTHER' == df[df['POS'] == 1984]['TYPE'].tolist()[0]


def test_write_bcf():
    sample_vcf = data_dir / 'SampleA' / 'snps.vcf.gz'
    with tempfile.TemporaryDirectory() as out_dir:
        out_file = Path(out_dir) / 'out.bcf'

        assert not out_file.exists()
        VariationFile(sample_vcf).write(out_file)
        assert out_file.exists()


def test_consensus_no_mask():
    sample_bcf = variation_dir / 'SampleA.bcf'

    expected_consensus_file = consensus_dir / 'SampleA-consensus-nomask.fasta.gz'
    name, expected_consensus_records = parse_sequence_file(expected_consensus_file)
    expected_consensus_record = expected_consensus_records[0]

    seq_records = VariationFile(sample_bcf).consensus(reference_file=reference_file)

    assert 1 == len(seq_records)
    actual_seq_record = seq_records[0]
    assert 5180 == len(actual_seq_record)
    assert expected_consensus_record.id == actual_seq_record.id
    assert expected_consensus_record.seq == actual_seq_record.seq


def test_consensus_empty_mask():
    sample_bcf = variation_dir / 'SampleA.bcf'
    empty_mask = MaskedGenomicRegions.empty_mask()

    expected_consensus_file = consensus_dir / 'SampleA-consensus-nomask.fasta.gz'
    name, expected_consensus_records = parse_sequence_file(expected_consensus_file)
    expected_consensus_record = expected_consensus_records[0]

    with tempfile.TemporaryDirectory() as out_dir:
        mask_file = Path(out_dir) / 'mask.bed.gz'
        empty_mask.write(mask_file)

        seq_records = VariationFile(sample_bcf).consensus(reference_file=reference_file,
                                                          mask_file=mask_file)
        assert 1 == len(seq_records)
        actual_seq_record = seq_records[0]
        assert 5180 == len(actual_seq_record)
        assert expected_consensus_record.id == actual_seq_record.id
        assert expected_consensus_record.seq == actual_seq_record.seq


def test_consensus_mask():
    sample_bcf = variation_dir / 'SampleA.bcf'
    sample_mask_file = variation_dir / 'SampleA.bed.gz'

    expected_consensus_file = consensus_dir / 'SampleA-consensus-withmask.fasta.gz'
    name, expected_consensus_records = parse_sequence_file(expected_consensus_file)
    expected_consensus_record = expected_consensus_records[0]

    seq_records = VariationFile(sample_bcf).consensus(reference_file=reference_file,
                                                      mask_file=sample_mask_file)
    assert 1 == len(seq_records)
    actual_seq_record = seq_records[0]
    assert 5180 == len(actual_seq_record)
    assert expected_consensus_record.id == actual_seq_record.id
    assert expected_consensus_record.seq == actual_seq_record.seq


def test_consensus_mask_over_mutation():
    sample_bcf = variation_dir / 'SampleA.bcf'
    sample_mask_file = variation_dir / 'SampleA-mask-over-mutation.bed.gz'

    expected_consensus_file = consensus_dir / 'SampleA-consensus-withmask-over-mutation.fasta.gz'
    name, expected_consensus_records = parse_sequence_file(expected_consensus_file)
    expected_consensus_record = expected_consensus_records[0]

    seq_records = VariationFile(sample_bcf).consensus(reference_file=reference_file,
                                                      mask_file=sample_mask_file)
    assert 1 == len(seq_records)
    actual_seq_record = seq_records[0]
    assert 5180 == len(actual_seq_record)
    assert expected_consensus_record.id == actual_seq_record.id
    assert expected_consensus_record.seq == actual_seq_record.seq


def test_union_all_files():
    # List like this to guarantee a specific order
    variant_files = [
        data_dir / 'SampleA' / 'snps.vcf.gz',
        data_dir / 'SampleB' / 'snps.vcf.gz',
        data_dir / 'SampleC' / 'snps.vcf.gz'
    ]
    union_df = VariationFile.union_all_files(variant_files)

    assert 60 == len(union_df)
    assert ['ID', 'CHROM', 'POS', 'REF', 'ALT', 'COUNT'] == union_df.columns.tolist()
    assert 1 == union_df[(union_df['CHROM'] == 'reference') & (union_df['POS'] == 190)]['COUNT'].values[0]
    assert 1 == union_df[(union_df['CHROM'] == 'reference') & (union_df['POS'] == 5061)]['COUNT'].values[0]
    assert 2 == union_df[(union_df['CHROM'] == 'reference') & (union_df['POS'] == 4975)]['COUNT'].values[0]
    assert 1 == union_df[(union_df['CHROM'] == 'reference') & (union_df['POS'] == 2076)]['COUNT'].values[0]


def test_union_one_file():
    sample_bcf = variation_dir / 'SampleA.bcf'
    union_df = VariationFile.union_all_files([sample_bcf])

    assert 26 == len(union_df)
    assert ['ID', 'CHROM', 'POS', 'REF', 'ALT', 'COUNT'] == union_df.columns.tolist()
    assert 1 == union_df[(union_df['CHROM'] == 'reference') & (union_df['POS'] == 293)]['COUNT'].values[0]
    assert 1 == union_df[(union_df['CHROM'] == 'reference') & (union_df['POS'] == 4929)]['COUNT'].values[0]


def test_union_batch_size_1():
    # List like this to guarantee a specific order
    variant_files = [
        data_dir / 'SampleA' / 'snps.vcf.gz',
        data_dir / 'SampleB' / 'snps.vcf.gz',
        data_dir / 'SampleC' / 'snps.vcf.gz'
    ]
    union_df = VariationFile.union_all_files(variant_files, batch_size=1)

    assert 60 == len(union_df)
    assert ['ID', 'CHROM', 'POS', 'REF', 'ALT', 'COUNT'] == union_df.columns.tolist()
    assert 1 == union_df[(union_df['CHROM'] == 'reference') & (union_df['POS'] == 190)]['COUNT'].values[0]
    assert 1 == union_df[(union_df['CHROM'] == 'reference') & (union_df['POS'] == 5061)]['COUNT'].values[0]
    assert 2 == union_df[(union_df['CHROM'] == 'reference') & (union_df['POS'] == 4975)]['COUNT'].values[0]
    assert 1 == union_df[(union_df['CHROM'] == 'reference') & (union_df['POS'] == 2076)]['COUNT'].values[0]


def test_union_batch_size_2_all_data():
    # List like this to guarantee a specific order
    variant_files = [
        data_dir / 'SampleA' / 'snps.vcf.gz',
        data_dir / 'SampleB' / 'snps.vcf.gz',
        data_dir / 'SampleC' / 'snps.vcf.gz'
    ]
    union_df = VariationFile.union_all_files(variant_files, include_expression=None, batch_size=2)

    assert 112 == len(union_df)
    assert ['ID', 'CHROM', 'POS', 'REF', 'ALT', 'COUNT'] == union_df.columns.tolist()
    assert 1 == union_df[(union_df['CHROM'] == 'reference') & (union_df['POS'] == 190)]['COUNT'].values[0]
    assert 1 == union_df[(union_df['CHROM'] == 'reference') & (union_df['POS'] == 5061)]['COUNT'].values[0]
    assert 2 == union_df[(union_df['CHROM'] == 'reference') & (union_df['POS'] == 4975)]['COUNT'].values[0]
    assert 1 == union_df[(union_df['CHROM'] == 'reference') & (union_df['POS'] == 2076)]['COUNT'].values[0]


def test_union_many_files_batch_size_2_more_data():
    # List like this to guarantee a specific order
    variant_files = [
        data_dir / 'SampleA' / 'snps.vcf.gz',
        data_dir / 'SampleB' / 'snps.vcf.gz',
        data_dir / 'SampleC' / 'snps.vcf.gz',
        extra_snippy_dir / 'SampleB2.snps.fill-tags.vcf.gz',
        extra_snippy_dir / 'SampleB3.snps.fill-tags.vcf.gz',
        extra_snippy_dir / 'SampleB4.snps.fill-tags.vcf.gz',
        extra_snippy_dir / 'SampleB5.snps.fill-tags.vcf.gz',
        extra_snippy_dir / 'SampleB5-different-allele.fill-tags.vcf.gz',
    ]
    union_df = VariationFile.union_all_files(variant_files, batch_size=2, include_expression=None)

    assert 115 == len(union_df)
    assert ['ID', 'CHROM', 'POS', 'REF', 'ALT', 'COUNT'] == union_df.columns.tolist()
    assert 6 == union_df[union_df['ID'] == 'reference:190:A:G']['COUNT'].values[0]
    assert 6 == union_df[union_df['ID'] == 'reference:5061:G:A']['COUNT'].values[0]

    assert 6 == union_df[union_df['ID'] == 'reference:4975:T:C']['COUNT'].values[0]
    assert 1 == union_df[union_df['ID'] == 'reference:4975:T:CAT']['COUNT'].values[0]

    assert 1 == union_df[union_df['ID'] == 'reference:2076:A:T']['COUNT'].values[0]

    assert 2 == union_df[union_df['ID'] == 'reference:1483:AAAGAGGGGCTGCTGGAGCCG:A']['COUNT'].values[0]
    assert 3 == union_df[union_df['ID'] == 'reference:1640:C:A']['COUNT'].values[0]

    assert 6 == union_df[union_df['ID'] == 'reference:4693:C:CGA']['COUNT'].values[0]
    assert 1 == union_df[union_df['ID'] == 'reference:4693:C:G']['COUNT'].values[0]

    assert 3 == union_df[union_df['ID'] == 'reference:883:CACATG:C']['COUNT'].values[0]

    assert 5 == union_df[union_df['ID'] == 'reference:349:AAGT:A']['COUNT'].values[0]
    assert 1 == union_df[union_df['ID'] == 'reference:349:AAGT:T']['COUNT'].values[0]


def test_union_many_files_batch_size_2_with_empty_vcf():
    # List like this to guarantee a specific order
    variant_files = [
        data_dir / 'SampleA' / 'snps.vcf.gz',
        data_dir / 'SampleB' / 'snps.vcf.gz',
        data_dir / 'SampleC' / 'snps.vcf.gz',
        extra_snippy_dir / 'SampleB2.snps.fill-tags.vcf.gz',
        extra_snippy_dir / 'SampleB3.snps.fill-tags.vcf.gz',
        extra_snippy_dir / 'SampleB4.snps.fill-tags.vcf.gz',
        extra_snippy_dir / 'SampleB5.snps.fill-tags.vcf.gz',
        extra_snippy_dir / 'SampleB5-different-allele.fill-tags.vcf.gz',
        extra_snippy_dir / 'SampleB-empty.snps.fill-tags.vcf.gz'
    ]
    union_df = VariationFile.union_all_files(variant_files, batch_size=2, include_expression=None)
    print(union_df)

    assert 115 == len(union_df)
    assert ['ID', 'CHROM', 'POS', 'REF', 'ALT', 'COUNT'] == union_df.columns.tolist()
    assert 6 == union_df[union_df['ID'] == 'reference:190:A:G']['COUNT'].values[0]
    assert 6 == union_df[union_df['ID'] == 'reference:5061:G:A']['COUNT'].values[0]

    assert 6 == union_df[union_df['ID'] == 'reference:4975:T:C']['COUNT'].values[0]
    assert 1 == union_df[union_df['ID'] == 'reference:4975:T:CAT']['COUNT'].values[0]

    assert 1 == union_df[union_df['ID'] == 'reference:2076:A:T']['COUNT'].values[0]

    assert 2 == union_df[union_df['ID'] == 'reference:1483:AAAGAGGGGCTGCTGGAGCCG:A']['COUNT'].values[0]
    assert 3 == union_df[union_df['ID'] == 'reference:1640:C:A']['COUNT'].values[0]

    assert 6 == union_df[union_df['ID'] == 'reference:4693:C:CGA']['COUNT'].values[0]
    assert 1 == union_df[union_df['ID'] == 'reference:4693:C:G']['COUNT'].values[0]

    assert 3 == union_df[union_df['ID'] == 'reference:883:CACATG:C']['COUNT'].values[0]

    assert 5 == union_df[union_df['ID'] == 'reference:349:AAGT:A']['COUNT'].values[0]
    assert 1 == union_df[union_df['ID'] == 'reference:349:AAGT:T']['COUNT'].values[0]


def test_union_many_files_batch_size_2_single_empty_vcf():
    # List like this to guarantee a specific order
    variant_files = [
        extra_snippy_dir / 'SampleB-empty.snps.fill-tags.vcf.gz'
    ]
    union_df = VariationFile.union_all_files(variant_files, batch_size=2, include_expression=None)
    print(union_df)

    assert 0 == len(union_df)
    assert ['ID', 'CHROM', 'POS', 'REF', 'ALT', 'COUNT'] == union_df.columns.tolist()
