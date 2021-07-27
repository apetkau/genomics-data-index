import tempfile
from pathlib import Path

import pytest

from genomics_data_index.storage.MaskedGenomicRegions import MaskedGenomicRegions
from genomics_data_index.storage.io.mutation.SequenceFile import SequenceFile
from genomics_data_index.storage.io.mutation.VariationFile import VariationFile
from genomics_data_index.storage.io.mutation.VcfSnpEffAnnotationParser import VcfSnpEffAnnotationParser
from genomics_data_index.test.integration import data_dir, regular_vcf_dir, variation_dir, reference_file, consensus_dir
from genomics_data_index.test.integration import extra_snippy_dir
from genomics_data_index.test.integration import reference_file_5000_snpeff, snpeff_vcf_file
from genomics_data_index.test.integration import snpeff_sample_vcfs
from genomics_data_index.test.integration.storage.io import read_vcf_df


@pytest.fixture
def snpeff_parser() -> VcfSnpEffAnnotationParser:
    return VcfSnpEffAnnotationParser()


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
    name, expected_consensus_records = SequenceFile(expected_consensus_file).parse_sequence_file()
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
    name, expected_consensus_records = SequenceFile(expected_consensus_file).parse_sequence_file()
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
    name, expected_consensus_records = SequenceFile(expected_consensus_file).parse_sequence_file()
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
    name, expected_consensus_records = SequenceFile(expected_consensus_file).parse_sequence_file()
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
    union_df = VariationFile.union_all_files(variant_files, include_expression='TYPE="SNP"')

    assert 60 == len(union_df)
    assert ['ID', 'CHROM', 'POS', 'REF', 'ALT', 'COUNT'] == union_df.columns.tolist()
    assert 1 == union_df[(union_df['CHROM'] == 'reference') & (union_df['POS'] == 190)]['COUNT'].values[0]
    assert 1 == union_df[(union_df['CHROM'] == 'reference') & (union_df['POS'] == 5061)]['COUNT'].values[0]
    assert 2 == union_df[(union_df['CHROM'] == 'reference') & (union_df['POS'] == 4975)]['COUNT'].values[0]
    assert 1 == union_df[(union_df['CHROM'] == 'reference') & (union_df['POS'] == 2076)]['COUNT'].values[0]


def test_union_one_file():
    sample_bcf = variation_dir / 'SampleA.bcf'
    union_df = VariationFile.union_all_files([sample_bcf], include_expression='TYPE="SNP"')

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
    union_df = VariationFile.union_all_files(variant_files, include_expression='TYPE="SNP"', batch_size=1)

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
    union_df = VariationFile.union_all_files(variant_files, batch_size=2)

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
    union_df = VariationFile.union_all_files(variant_files, batch_size=2)

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
    union_df = VariationFile.union_all_files(variant_files, batch_size=2)
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


def test_union_many_files_batch_size_odd_cores_3():
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
    union_df = VariationFile.union_all_files(variant_files, ncores=3, batch_size=3)
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


def test_union_many_files_ambiguous():
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
        extra_snippy_dir / 'SampleB5-different-allele-ambiguous.vcf.gz',
    ]
    union_df = VariationFile.union_all_files(variant_files)
    print(union_df)

    assert 119 == len(union_df)
    assert ['ID', 'CHROM', 'POS', 'REF', 'ALT', 'COUNT'] == union_df.columns.tolist()
    assert 6 == union_df[union_df['ID'] == 'reference:190:A:G']['COUNT'].values[0]
    assert 1 == union_df[union_df['ID'] == 'reference:190:A:N']['COUNT'].values[0]

    assert 7 == union_df[union_df['ID'] == 'reference:5061:G:A']['COUNT'].values[0]
    assert 2 == union_df[union_df['ID'] == 'reference:1483:AAAGAGGGGCTGCTGGAGCCG:A']['COUNT'].values[0]

    assert 5 == union_df[union_df['ID'] == 'reference:349:AAGT:A']['COUNT'].values[0]
    assert 1 == union_df[union_df['ID'] == 'reference:349:AAGT:T']['COUNT'].values[0]
    assert 1 == union_df[union_df['ID'] == 'reference:349:ANGT:T']['COUNT'].values[0]

    assert 6 == union_df[union_df['ID'] == 'reference:4693:C:CGA']['COUNT'].values[0]
    assert 1 == union_df[union_df['ID'] == 'reference:4693:C:G']['COUNT'].values[0]
    assert 1 == union_df[union_df['ID'] == 'reference:4693:N:G']['COUNT'].values[0]

    assert 6 == union_df[union_df['ID'] == 'reference:4975:T:C']['COUNT'].values[0]
    assert 1 == union_df[union_df['ID'] == 'reference:4975:T:CAT']['COUNT'].values[0]
    assert 1 == union_df[union_df['ID'] == 'reference:4975:T:CNT']['COUNT'].values[0]


def test_union_many_files_batch_size_2_single_empty_vcf():
    # List like this to guarantee a specific order
    variant_files = [
        extra_snippy_dir / 'SampleB-empty.snps.fill-tags.vcf.gz'
    ]
    union_df = VariationFile.union_all_files(variant_files, batch_size=2)
    print(union_df)

    assert 0 == len(union_df)
    assert ['ID', 'CHROM', 'POS', 'REF', 'ALT', 'COUNT'] == union_df.columns.tolist()


def test_read_features(snpeff_parser):
    vcf_file = data_dir / 'SampleA' / 'snps.vcf.gz'
    df = VariationFile(vcf_file).read_features('SampleA', snpeff_parser=snpeff_parser)

    assert 46 == len(df), 'Data fram has incorrect length'

    assert {'snps.vcf.gz'} == set(df['FILE'].tolist()), 'Incorrect filename'
    assert {'SampleA'} == set(df['SAMPLE'].tolist()), 'Incorrect sample name'

    v = df[df['POS'] == 461]
    assert 'AAAT' == v['REF'].values[0], 'Incorrect reference'
    assert 'G' == v['ALT'].values[0], 'Incorrect alt'

    v = df[df['POS'] == 1048]
    assert 'C' == v['REF'].values[0], 'Incorrect reference'
    assert 'G' == v['ALT'].values[0], 'Incorrect alt'

    v = df[df['POS'] == 1253]
    assert 'T' == v['REF'].values[0], 'Incorrect reference'
    assert 'TAA' == v['ALT'].values[0], 'Incorrect alt'

    v = df[df['POS'] == 3656]
    assert 'CATT' == v['REF'].values[0], 'Incorrect reference'
    assert 'C' == v['ALT'].values[0], 'Incorrect alt'


def test_read_features_snpeff(snpeff_parser):
    sample_10_014 = VariationFile(
        snpeff_sample_vcfs['SH10-014']).read_features('SH10-014', snpeff_parser=snpeff_parser).sort_values('POS')
    sample_14_001 = VariationFile(
        snpeff_sample_vcfs['SH14-001']).read_features('SH14-001', snpeff_parser=snpeff_parser).sort_values('POS')
    sample_14_014 = VariationFile(
        snpeff_sample_vcfs['SH14-014']).read_features('SH14-014', snpeff_parser=snpeff_parser).sort_values('POS')

    assert 139 == len(sample_10_014)
    assert ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
            'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p'] == list(sample_10_014.columns)

    # snv/snp
    sample_10_014_varA = sample_10_014[sample_10_014['POS'] == 140658]
    assert 1 == len(sample_10_014_varA)
    assert ['SH10-014', 'NC_011083', 140658, 'C', 'A', 'snp', 'SH10-014.vcf.gz', 'NC_011083:140658:C:A',
            'A', 'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
            'c.497C>A', 'p.Ala166Glu'] == sample_10_014_varA[
               sample_10_014_varA['ANN.Annotation'] == 'missense_variant'].iloc[0].tolist()

    # del
    sample_10_014_varB = sample_10_014[sample_10_014['POS'] == 1125996]
    assert 1 == len(sample_10_014_varB)
    assert ['SH10-014', 'NC_011083', 1125996, 'CG', 'C', 'del', 'SH10-014.vcf.gz', 'NC_011083:1125996:CG:C',
            'C', 'frameshift_variant', 'HIGH', 'SEHA_RS05995', 'SEHA_RS05995', 'transcript', 'protein_coding',
            'c.418delG', 'p.Glu140fs'] == sample_10_014_varB[
               sample_10_014_varB['ANN.Annotation'] == 'frameshift_variant'].iloc[0].tolist()

    # ins
    sample_10_014_varC = sample_10_014[sample_10_014['POS'] == 1246085]
    assert 1 == len(sample_10_014_varC)
    assert ['SH10-014', 'NC_011083', 1246085, 'C', 'CG', 'ins', 'SH10-014.vcf.gz', 'NC_011083:1246085:C:CG',
            'CG', 'frameshift_variant', 'HIGH', 'mdtG', 'SEHA_RS06605', 'transcript', 'protein_coding',
            'c.722dupC', 'p.Leu242fs'] == sample_10_014_varC[
               sample_10_014_varC['ANN.Annotation'] == 'frameshift_variant'].iloc[0].tolist()

    # complex
    sample_10_014_varD = sample_10_014[sample_10_014['POS'] == 3535121]
    assert 1 == len(sample_10_014_varD)
    assert ['SH10-014', 'NC_011083', 3535121, 'CGCGA', 'TGTGG', 'complex', 'SH10-014.vcf.gz',
            'NC_011083:3535121:CGCGA:TGTGG',
            'TGTGG', 'missense_variant', 'MODERATE', 'oadA', 'SEHA_RS17780', 'transcript', 'protein_coding',
            'c.1119_1123delTCGCGinsCCACA', 'p.ArgAla374HisThr'] == sample_10_014_varD[
               sample_10_014_varD['ANN.Annotation'] == 'missense_variant'].iloc[0].tolist()

    assert 115 == len(sample_14_001)
    assert ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
            'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p'] == list(sample_14_001.columns)
    sample_14_001_var = sample_14_001[sample_14_001['POS'] == 140658]
    assert 1 == len(sample_14_001_var)
    assert ['SH14-001', 'NC_011083', 140658, 'C', 'A', 'snp', 'SH14-001.vcf.gz', 'NC_011083:140658:C:A',
            'A', 'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
            'c.497C>A', 'p.Ala166Glu'] == sample_14_001_var[
               sample_14_001_var['ANN.Annotation'] == 'missense_variant'].iloc[0].tolist()

    assert 107 == len(sample_14_014)
    assert ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
            'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p'] == list(sample_14_014.columns)
    sample_14_014_var = sample_14_014[sample_14_014['POS'] == 298472]
    assert 1 == len(sample_14_014_var)
    assert ['SH14-014', 'NC_011083', 298472, 'A', 'C', 'snp', 'SH14-014.vcf.gz', 'NC_011083:298472:A:C',
            'C', 'intergenic_region', 'MODIFIER', 'SEHA_RS01880-SEHA_RS01885', 'SEHA_RS01880-SEHA_RS01885',
            'intergenic_region', 'n.298472A>C'] == sample_14_014_var[
               sample_14_014_var['ANN.Annotation'] == 'intergenic_region'].drop(
        ['ANN.Transcript_BioType', 'ANN.HGVS.p'], axis='columns').iloc[0].tolist()
    assert {True} == set(sample_14_014_var[sample_14_014_var['ANN.Annotation'] == 'intergenic_region'] \
                             [['ANN.Transcript_BioType', 'ANN.HGVS.p']].iloc[0].isna().tolist())


def test_annotate(snpeff_parser):
    with tempfile.TemporaryDirectory() as out_dir:
        database_dir = Path(out_dir)
        output_vcf_file = database_dir / 'output.vcf.gz'
        variation_file = VariationFile(snpeff_vcf_file)

        snpeff_database = SequenceFile(reference_file_5000_snpeff).create_snpeff_database(database_dir)
        annotated_variation_file = variation_file.annotate(snpeff_database=snpeff_database,
                                                           annotated_vcf=output_vcf_file)

        assert output_vcf_file == annotated_variation_file.file

        # Verify VCF annotation contents
        vcf_annotation_df = annotated_variation_file.read_features('SampleA', snpeff_parser=snpeff_parser).sort_values(
            'POS')
        assert 2 == len(vcf_annotation_df)
        assert ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
                'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
                'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p'] == list(
            vcf_annotation_df.columns)
        assert ['NC_011083.1:195:C:G', 'NC_011083.1:207:C:G'] == vcf_annotation_df['VARIANT_ID'].tolist()
        assert ['SNP', 'SNP'] == vcf_annotation_df['TYPE'].tolist()
        assert ['missense_variant', 'synonymous_variant'] == vcf_annotation_df['ANN.Annotation'].tolist()
        assert ['SEHA_RS00560', 'SEHA_RS00560'] == vcf_annotation_df['ANN.Gene_ID'].tolist()
        assert ['thrL', 'thrL'] == vcf_annotation_df['ANN.Gene_Name'].tolist()
        assert ['c.6C>G', 'c.18C>G'] == vcf_annotation_df['ANN.HGVS.c'].tolist()
        assert ['p.N2K', 'p.T6T'] == vcf_annotation_df['ANN.HGVS.p'].tolist()

        # Original file should still exist and be unannotated
        vcf_no_annotation_df = variation_file.read_features('SampleA', snpeff_parser=snpeff_parser).sort_values('POS')
        assert 2 == len(vcf_no_annotation_df)
        assert ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
                'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
                'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p'] == list(
            vcf_no_annotation_df.columns)
        assert ['NC_011083.1:195:C:G', 'NC_011083.1:207:C:G'] == vcf_no_annotation_df['VARIANT_ID'].tolist()
        assert ['SNP', 'SNP'] == vcf_no_annotation_df['TYPE'].tolist()
        assert all(vcf_no_annotation_df['ANN.Annotation'].isna())
        assert all(vcf_no_annotation_df['ANN.Gene_ID'].isna())
        assert all(vcf_no_annotation_df['ANN.HGVS.c'].isna())
        assert all(vcf_no_annotation_df['ANN.HGVS.p'].isna())
