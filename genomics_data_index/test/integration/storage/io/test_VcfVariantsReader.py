import tempfile
from os import path, listdir
from pathlib import Path
from typing import List, Dict, cast
import pandas as pd

import pytest

from genomics_data_index.storage.io.mutation.NucleotideSampleData import NucleotideSampleData
from genomics_data_index.storage.io.mutation.NucleotideSampleDataPackage import NucleotideSampleDataPackage
from genomics_data_index.storage.io.mutation.VcfVariantsReader import VcfVariantsReader
from genomics_data_index.storage.io.processor.SerialSampleFilesProcessor import SerialSampleFilesProcessor
from genomics_data_index.test.integration import data_dir
from genomics_data_index.test.integration import data_dir_empty
from genomics_data_index.test.integration import snpeff_sample_vcfs


@pytest.fixture
def sample_dirs() -> List[Path]:
    return [data_dir / d for d in listdir(data_dir) if path.isdir(data_dir / d)]


@pytest.fixture
def sample_dirs_empty() -> List[Path]:
    return [data_dir_empty / d for d in listdir(data_dir_empty) if path.isdir(data_dir_empty / d)]


def vcf_and_mask_files(sample_dirs) -> Dict[str, Dict[str, Path]]:
    sample_vcf_map = {}
    sample_genomic_files_mask = {}
    for d in sample_dirs:
        sample_name = path.basename(d)
        vcf_file = Path(d, 'snps.vcf.gz')
        genomic_file_mask = Path(d, 'snps.aligned.fa')

        sample_vcf_map[sample_name] = vcf_file
        sample_genomic_files_mask[sample_name] = genomic_file_mask

    return {
        'vcfs': sample_vcf_map,
        'masks': sample_genomic_files_mask
    }


def variants_reader_internal(sample_dirs) -> VcfVariantsReader:
    tmp_dir = Path(tempfile.mkdtemp())
    vcf_masks = vcf_and_mask_files(sample_dirs)
    file_processor = SerialSampleFilesProcessor(tmp_dir)
    data_package = NucleotideSampleDataPackage.create_from_sequence_masks(sample_vcf_map=vcf_masks['vcfs'],
                                                                          masked_genomic_files_map=vcf_masks[
                                                                              'masks'],
                                                                          sample_files_processor=file_processor)
    processed_files = cast(Dict[str, NucleotideSampleData], data_package.process_all_data())
    return VcfVariantsReader.create(processed_files)


@pytest.fixture
def variants_reader(sample_dirs) -> VcfVariantsReader:
    return variants_reader_internal(sample_dirs)


@pytest.fixture
def variants_reader_empty(sample_dirs_empty) -> VcfVariantsReader:
    return variants_reader_internal(sample_dirs_empty)


@pytest.fixture
def variants_reader_default_no_data() -> VcfVariantsReader:
    return VcfVariantsReader(sample_files_map={})


@pytest.fixture
def variants_reader_empty_masks(sample_dirs) -> VcfVariantsReader:
    sample_vcf_map = {}
    for d in sample_dirs:
        sample_name = path.basename(d)
        vcf_file = Path(d, 'snps.vcf.gz')

        sample_vcf_map[sample_name] = vcf_file

    tmp_dir = Path(tempfile.mkdtemp())
    file_processor = SerialSampleFilesProcessor(tmp_dir)
    data_package = NucleotideSampleDataPackage.create_from_sequence_masks(sample_vcf_map=sample_vcf_map,
                                                                          sample_files_processor=file_processor)
    processed_files = cast(Dict[str, NucleotideSampleData], data_package.process_all_data())
    return VcfVariantsReader.create(processed_files)


def variants_reader_from_snippy_internal(sample_dirs) -> VcfVariantsReader:
    tmp_dir = Path(tempfile.mkdtemp())
    file_processor = SerialSampleFilesProcessor(tmp_dir)
    data_package = NucleotideSampleDataPackage.create_from_snippy(sample_dirs=sample_dirs,
                                                                  sample_files_processor=file_processor)
    processed_files = cast(Dict[str, NucleotideSampleData], data_package.process_all_data())
    return VcfVariantsReader.create(processed_files)


@pytest.fixture
def variants_reader_from_snippy(sample_dirs) -> VcfVariantsReader:
    return variants_reader_from_snippy_internal(sample_dirs)


@pytest.fixture
def variants_reader_snpeff() -> VcfVariantsReader:
    tmp_dir = Path(tempfile.mkdtemp())
    file_processor = SerialSampleFilesProcessor(tmp_dir)
    data_package = NucleotideSampleDataPackage.create_from_sequence_masks(sample_vcf_map=snpeff_sample_vcfs,
                                                                          sample_files_processor=file_processor)
    processed_files = cast(Dict[str, NucleotideSampleData], data_package.process_all_data())
    return VcfVariantsReader.create(processed_files)


def test_get_variants_table(variants_reader):
    df = variants_reader.get_features_table()

    assert 129 == len(df), 'Data has incorrect length'
    assert {'SampleA', 'SampleB', 'SampleC'} == set(df['SAMPLE'].tolist()), 'Incorrect sample names'

    assert ['SNP'] == df.loc[(df['SAMPLE'] == 'SampleA') & (df['POS'] == 293), 'TYPE'].tolist()
    assert ['INDEL'] == df.loc[(df['SAMPLE'] == 'SampleA') & (df['POS'] == 302), 'TYPE'].tolist()
    assert ['INDEL'] == df.loc[(df['SAMPLE'] == 'SampleA') & (df['POS'] == 324), 'TYPE'].tolist()
    assert ['INDEL'] == df.loc[(df['SAMPLE'] == 'SampleA') & (df['POS'] == 374), 'TYPE'].tolist()
    assert ['OTHER'] == df.loc[(df['SAMPLE'] == 'SampleA') & (df['POS'] == 461), 'TYPE'].tolist()
    assert ['OTHER'] == df.loc[(df['SAMPLE'] == 'SampleB') & (df['POS'] == 1325), 'TYPE'].tolist()
    assert ['OTHER'] == df.loc[(df['SAMPLE'] == 'SampleC') & (df['POS'] == 1984), 'TYPE'].tolist()


def test_get_genomic_masks(variants_reader):
    mask = variants_reader.get_genomic_masked_region('SampleA')
    assert 437 == len(mask)
    assert {'reference'} == mask.sequence_names()

    mask = variants_reader.get_genomic_masked_region('SampleB')
    assert 276 == len(mask)
    assert {'reference'} == mask.sequence_names()

    mask = variants_reader.get_genomic_masked_region('SampleC')
    assert 329 == len(mask)
    assert {'reference'} == mask.sequence_names()


def test_get_genomic_masks_empty(variants_reader_empty_masks):
    mask = variants_reader_empty_masks.get_genomic_masked_region('SampleA')
    assert mask.is_empty()

    mask = variants_reader_empty_masks.get_genomic_masked_region('SampleB')
    assert mask.is_empty()

    mask = variants_reader_empty_masks.get_genomic_masked_region('SampleC')
    assert mask.is_empty()


def test_get_samples_list(variants_reader):
    assert {'SampleA', 'SampleB', 'SampleC'} == set(variants_reader.samples_list())


def test_get_variants_table_empty(variants_reader_empty):
    df = variants_reader_empty.get_features_table()

    assert 0 == len(df), 'Data has incorrect length'


def test_get_or_create_feature_file(variants_reader):
    file = variants_reader.get_or_create_feature_file('SampleA')
    assert file.exists()


def test_snippy_read_vcf(variants_reader_from_snippy):
    vcf_file = data_dir / 'SampleA' / 'snps.vcf.gz'

    df = variants_reader_from_snippy.read_vcf(vcf_file, 'SampleA')

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


def test_snippy_get_variants_table(variants_reader_from_snippy):
    df = variants_reader_from_snippy.get_features_table()

    assert 129 == len(df), 'Data has incorrect length'
    assert {'SampleA', 'SampleB', 'SampleC'} == set(df['SAMPLE'].tolist()), 'Incorrect sample names'


def test_snippy_get_genomic_masks(variants_reader_from_snippy):
    mask = variants_reader_from_snippy.get_genomic_masked_region('SampleA')
    assert 437 == len(mask)
    assert {'reference'} == mask.sequence_names()

    mask = variants_reader_from_snippy.get_genomic_masked_region('SampleB')
    assert 276 == len(mask)
    assert {'reference'} == mask.sequence_names()

    mask = variants_reader_from_snippy.get_genomic_masked_region('SampleC')
    assert 329 == len(mask)
    assert {'reference'} == mask.sequence_names()


def test_snippy_get_samples_list(variants_reader_from_snippy):
    assert {'SampleA', 'SampleB', 'SampleC'} == set(variants_reader_from_snippy.samples_list())


def test_snippy_get_samples_list_two_files():
    sample_dirs = [data_dir / 'SampleA', data_dir / 'SampleB']
    reader = variants_reader_from_snippy_internal(sample_dirs)

    assert {'SampleA', 'SampleB'} == set(reader.samples_list())


def test_snippy_read_empty_vcf(sample_dirs_empty):
    reader = variants_reader_from_snippy_internal(sample_dirs_empty)
    df = reader.get_features_table()

    assert 0 == len(df), 'Data has incorrect length'


def test_read_snpeff(variants_reader_default_no_data: VcfVariantsReader):
    vr = variants_reader_default_no_data
    sample_10_014 = vr.read_vcf(file=snpeff_sample_vcfs['SH10-014'], sample_name='SH10-014').sort_values('POS')
    sample_14_001 = vr.read_vcf(file=snpeff_sample_vcfs['SH14-001'], sample_name='SH14-001').sort_values('POS')
    sample_14_014 = vr.read_vcf(file=snpeff_sample_vcfs['SH14-014'], sample_name='SH14-014').sort_values('POS')

    assert 1434 == len(sample_10_014)
    assert ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
            'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p'] == list(sample_10_014.columns)

    # snv/snp
    sample_10_014_varA = sample_10_014[sample_10_014['POS'] == 140658]
    assert 9 == len(sample_10_014_varA)
    assert ['SH10-014', 'NC_011083', 140658, 'C', 'A', 'snp', 'SH10-014.vcf.gz', 'NC_011083:140658:C:A',
            'A', 'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
            'c.497C>A', 'p.Ala166Glu'] == sample_10_014_varA[
        sample_10_014_varA['ANN.Annotation'] == 'missense_variant'].iloc[0].tolist()

    # del
    sample_10_014_varB = sample_10_014[sample_10_014['POS'] == 1125996]
    assert 14 == len(sample_10_014_varB)
    assert ['SH10-014', 'NC_011083', 1125996, 'CG', 'C', 'del', 'SH10-014.vcf.gz', 'NC_011083:1125996:CG:C',
            'C', 'frameshift_variant', 'HIGH', 'SEHA_RS05995', 'SEHA_RS05995', 'transcript', 'protein_coding',
            'c.418delG', 'p.Glu140fs'] == sample_10_014_varB[
        sample_10_014_varB['ANN.Annotation'] == 'frameshift_variant'].iloc[0].tolist()

    # ins
    sample_10_014_varC = sample_10_014[sample_10_014['POS'] == 1246085]
    assert 11 == len(sample_10_014_varC)
    assert ['SH10-014', 'NC_011083', 1246085, 'C', 'CG', 'ins', 'SH10-014.vcf.gz', 'NC_011083:1246085:C:CG',
            'CG', 'frameshift_variant', 'HIGH', 'mdtG', 'SEHA_RS06605', 'transcript', 'protein_coding',
            'c.722dupC', 'p.Leu242fs'] == sample_10_014_varC[
        sample_10_014_varC['ANN.Annotation'] == 'frameshift_variant'].iloc[0].tolist()

    # complex
    sample_10_014_varD = sample_10_014[sample_10_014['POS'] == 3535121]
    assert 10 == len(sample_10_014_varD)
    assert ['SH10-014', 'NC_011083', 3535121, 'CGCGA', 'TGTGG', 'complex', 'SH10-014.vcf.gz', 'NC_011083:3535121:CGCGA:TGTGG',
            'TGTGG', 'missense_variant', 'MODERATE', 'oadA', 'SEHA_RS17780', 'transcript', 'protein_coding',
            'c.1119_1123delTCGCGinsCCACA', 'p.ArgAla374HisThr'] == sample_10_014_varD[
        sample_10_014_varD['ANN.Annotation'] == 'missense_variant'].iloc[0].tolist()

    assert 1218 == len(sample_14_001)
    assert ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
            'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p'] == list(sample_14_001.columns)
    sample_14_001_var = sample_14_001[sample_14_001['POS'] == 140658]
    assert 9 == len(sample_14_001_var)
    assert ['SH14-001', 'NC_011083', 140658, 'C', 'A', 'snp', 'SH14-001.vcf.gz', 'NC_011083:140658:C:A',
            'A', 'missense_variant', 'MODERATE', 'murF', 'SEHA_RS01180', 'transcript', 'protein_coding',
            'c.497C>A', 'p.Ala166Glu'] == sample_14_001_var[
        sample_14_001_var['ANN.Annotation'] == 'missense_variant'].iloc[0].tolist()

    assert 1128 == len(sample_14_014)
    assert ['SAMPLE', 'CHROM', 'POS', 'REF', 'ALT', 'TYPE', 'FILE', 'VARIANT_ID',
            'ANN.Allele', 'ANN.Annotation', 'ANN.Annotation_Impact', 'ANN.Gene_Name', 'ANN.Gene_ID',
            'ANN.Feature_Type', 'ANN.Transcript_BioType', 'ANN.HGVS.c', 'ANN.HGVS.p'] == list(sample_14_014.columns)
    sample_14_014_var = sample_14_014[sample_14_014['POS'] == 298472]
    assert 6 == len(sample_14_014_var)
    assert ['SH14-014', 'NC_011083', 298472, 'A', 'C', 'snp', 'SH14-014.vcf.gz', 'NC_011083:298472:A:C',
            'C', 'intergenic_region', 'MODIFIER', 'SEHA_RS01880-SEHA_RS01885', 'SEHA_RS01880-SEHA_RS01885',
            'intergenic_region', 'n.298472A>C'] == sample_14_014_var[
        sample_14_014_var['ANN.Annotation'] == 'intergenic_region'].drop(
        ['ANN.Transcript_BioType', 'ANN.HGVS.p'], axis='columns').iloc[0].tolist()
    assert {True} == set(sample_14_014_var[sample_14_014_var['ANN.Annotation'] == 'intergenic_region']\
              [['ANN.Transcript_BioType', 'ANN.HGVS.p']].iloc[0].isna().tolist())
