import tempfile
from os import path, listdir
from pathlib import Path
from typing import List, Dict, cast

import pytest

from storage.test.integration.variant import data_dir
from storage.test.integration.variant import data_dir_empty
from storage.variant.io.mutation.VcfVariantsReader import VcfVariantsReader
from storage.variant.io.processor.SerialSampleFilesProcessor import SerialSampleFilesProcessor
from storage.variant.io.mutation.NucleotideSampleDataPackage import NucleotideSampleDataPackage
from storage.variant.io.mutation.NucleotideSampleData import NucleotideSampleData


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
    assert 'SampleA' in str(file)


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
