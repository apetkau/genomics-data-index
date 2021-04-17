from os import path, listdir
from pathlib import Path
from typing import List, Dict
from tempfile import TemporaryDirectory
import tempfile

import pytest

from storage.test.integration.variant import data_dir
from storage.test.integration.variant import data_dir_empty
from storage.variant.io.mutation.VcfVariantsReader import VcfVariantsReader
from storage.variant.io.processor.SerialSampleFilesProcessor import SerialSampleFilesProcessor
from storage.variant.io.processor.MultipleProcessSampleFilesProcessor import MultipleProcessSampleFilesProcessor


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
    return VcfVariantsReader.create_from_sequence_masks(sample_vcf_map=vcf_masks['vcfs'],
                                                        masked_genomic_files_map=vcf_masks['masks'],
                                                        sample_files_processor=SerialSampleFilesProcessor(tmp_dir))


@pytest.fixture
def variants_reader(sample_dirs) -> VcfVariantsReader:
    return variants_reader_internal(sample_dirs)


@pytest.fixture
def variants_reader_empty(sample_dirs_empty) -> VcfVariantsReader:
    return variants_reader_internal(sample_dirs_empty)


@pytest.fixture
def variants_reader_empty_masks(sample_dirs) -> VcfVariantsReader:
    tmp_dir = Path(tempfile.mkdtemp())
    sample_vcf_map = {}
    for d in sample_dirs:
        sample_name = path.basename(d)
        vcf_file = Path(d, 'snps.vcf.gz')

        sample_vcf_map[sample_name] = vcf_file

    return VcfVariantsReader.create_from_sequence_masks(sample_vcf_map=sample_vcf_map,
                                                        sample_files_processor=SerialSampleFilesProcessor(tmp_dir))


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


def test_with_serial_sample_files_processor(sample_dirs):
    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)
        vcf_masks = vcf_and_mask_files(sample_dirs)
        variants_reader = VcfVariantsReader.create_from_sequence_masks(sample_vcf_map=vcf_masks['vcfs'],
                                                        masked_genomic_files_map=vcf_masks['masks'],
                                                        sample_files_processor=SerialSampleFilesProcessor(tmp_file))

        mask = variants_reader.get_genomic_masked_region('SampleA')
        assert 437 == len(mask)
        assert {'reference'} == mask.sequence_names()

        mask = variants_reader.get_genomic_masked_region('SampleB')
        assert 276 == len(mask)
        assert {'reference'} == mask.sequence_names()

        mask = variants_reader.get_genomic_masked_region('SampleC')
        assert 329 == len(mask)
        assert {'reference'} == mask.sequence_names()


def test_with_multiprocess_sample_files_processor(sample_dirs):
    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)
        vcf_masks = vcf_and_mask_files(sample_dirs)
        file_processor = MultipleProcessSampleFilesProcessor(tmp_file, processing_cores=2)
        variants_reader = VcfVariantsReader.create_from_sequence_masks(sample_vcf_map=vcf_masks['vcfs'],
                                                        masked_genomic_files_map=vcf_masks['masks'],
                                                        sample_files_processor=file_processor)

        mask = variants_reader.get_genomic_masked_region('SampleA')
        assert 437 == len(mask)
        assert {'reference'} == mask.sequence_names()

        mask = variants_reader.get_genomic_masked_region('SampleB')
        assert 276 == len(mask)
        assert {'reference'} == mask.sequence_names()

        mask = variants_reader.get_genomic_masked_region('SampleC')
        assert 329 == len(mask)
        assert {'reference'} == mask.sequence_names()


def test_with_null_sample_files_processor(sample_dirs):
    vcf_masks = vcf_and_mask_files(sample_dirs)
    variants_reader = VcfVariantsReader.create_from_sequence_masks(sample_vcf_map=vcf_masks['vcfs'],
                                                    masked_genomic_files_map=vcf_masks['masks'])

    with pytest.raises(Exception) as execinfo:
        variants_reader.get_genomic_masked_region('SampleA')
    assert 'Sample mask file is not preprocessed for sample' in str(execinfo.value)


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
