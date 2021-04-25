import re
from os import path, listdir
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import List, Dict, cast

import pytest

from genomics_data_index.test.integration import data_dir
from genomics_data_index.test.integration import data_dir_empty
from genomics_data_index.storage.io.mutation.NucleotideSampleData import NucleotideSampleData
from genomics_data_index.storage.io.mutation.NucleotideSampleDataPackage import NucleotideSampleDataPackage
from genomics_data_index.storage.io.processor.MultipleProcessSampleFilesProcessor import MultipleProcessSampleFilesProcessor
from genomics_data_index.storage.io.processor.SerialSampleFilesProcessor import SerialSampleFilesProcessor


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


def test_iter_sample_data(sample_dirs):
    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)
        vcf_masks = vcf_and_mask_files(sample_dirs)
        file_processor = SerialSampleFilesProcessor(tmp_file)
        data_package = NucleotideSampleDataPackage.create_from_sequence_masks(sample_vcf_map=vcf_masks['vcfs'],
                                                                              masked_genomic_files_map=vcf_masks[
                                                                                  'masks'],
                                                                              sample_files_processor=file_processor)

        processed_files_dict = {}
        for sample_file in data_package.iter_sample_data():
            processed_files_dict[sample_file.sample_name] = sample_file

        assert 3 == len(processed_files_dict)

        sample_data = processed_files_dict['SampleA']
        assert isinstance(sample_data, NucleotideSampleData)
        sample_data = cast(NucleotideSampleData, sample_data)
        mask = sample_data.get_mask()
        assert 437 == len(mask)
        assert {'reference'} == mask.sequence_names()
        mask_file = sample_data.get_mask_file()
        assert mask_file.parent == tmp_file
        vcf_file, vcf_index = sample_data.get_vcf_file()
        assert vcf_file.exists()
        assert vcf_index.exists()
        assert vcf_file.parent == tmp_file
        assert vcf_index.parent == tmp_file


def test_persisted_sample_data_file_names(sample_dirs):
    def split_all_ext(file: Path) -> str:
        return file.name.split('.')[0]

    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)
        vcf_masks = vcf_and_mask_files(sample_dirs)
        file_processor = SerialSampleFilesProcessor(tmp_file)
        data_package = NucleotideSampleDataPackage.create_from_sequence_masks(sample_vcf_map=vcf_masks['vcfs'],
                                                                              masked_genomic_files_map=vcf_masks[
                                                                                  'masks'],
                                                                              sample_files_processor=file_processor)

        processed_files_dict = data_package.process_all_data()
        sample_data = cast(NucleotideSampleData, processed_files_dict['SampleA'])
        mask_file = sample_data.get_mask_file()
        vcf_file, vcf_index = sample_data.get_vcf_file()

        assert sample_data.sample_name_persistence == split_all_ext(mask_file)
        assert sample_data.sample_name_persistence == split_all_ext(vcf_file)
        assert sample_data.sample_name_persistence == split_all_ext(vcf_index)
        assert 32 == len(sample_data.sample_name_persistence)
        assert bool(re.match(r'^[0-9a-f]*$', sample_data.sample_name_persistence))


def test_with_serial_sample_files_processor(sample_dirs):
    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)
        vcf_masks = vcf_and_mask_files(sample_dirs)
        file_processor = SerialSampleFilesProcessor(tmp_file)
        data_package = NucleotideSampleDataPackage.create_from_sequence_masks(sample_vcf_map=vcf_masks['vcfs'],
                                                                              masked_genomic_files_map=vcf_masks[
                                                                                  'masks'],
                                                                              sample_files_processor=file_processor)

        processed_files_dict = data_package.process_all_data()

        assert 3 == len(processed_files_dict)

        sample_data = processed_files_dict['SampleA']
        assert isinstance(sample_data, NucleotideSampleData)
        sample_data = cast(NucleotideSampleData, sample_data)
        mask = sample_data.get_mask()
        assert 437 == len(mask)
        assert {'reference'} == mask.sequence_names()
        assert sample_data.get_mask_file().parent == tmp_file
        vcf_file, vcf_index = sample_data.get_vcf_file()
        assert vcf_file.exists()
        assert vcf_index.exists()
        assert vcf_file.parent == tmp_file
        assert vcf_index.parent == tmp_file

        sample_data = processed_files_dict['SampleB']
        assert isinstance(sample_data, NucleotideSampleData)
        sample_data = cast(NucleotideSampleData, sample_data)
        mask = sample_data.get_mask()
        assert 276 == len(mask)
        assert {'reference'} == mask.sequence_names()
        assert sample_data.get_mask_file().parent == tmp_file
        vcf_file, vcf_index = sample_data.get_vcf_file()
        assert vcf_file.exists()
        assert vcf_index.exists()
        assert vcf_file.parent == tmp_file
        assert vcf_index.parent == tmp_file

        sample_data = processed_files_dict['SampleC']
        assert isinstance(sample_data, NucleotideSampleData)
        sample_data = cast(NucleotideSampleData, sample_data)
        mask = sample_data.get_mask()
        assert 329 == len(mask)
        assert {'reference'} == mask.sequence_names()
        assert sample_data.get_mask_file().parent == tmp_file
        vcf_file, vcf_index = sample_data.get_vcf_file()
        assert vcf_file.exists()
        assert vcf_index.exists()
        assert vcf_file.parent == tmp_file
        assert vcf_index.parent == tmp_file


def test_with_multiprocess_sample_files_processor(sample_dirs):
    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)
        vcf_masks = vcf_and_mask_files(sample_dirs)
        file_processor = MultipleProcessSampleFilesProcessor(tmp_file, processing_cores=2)
        data_package = NucleotideSampleDataPackage.create_from_sequence_masks(sample_vcf_map=vcf_masks['vcfs'],
                                                                              masked_genomic_files_map=vcf_masks[
                                                                                  'masks'],
                                                                              sample_files_processor=file_processor)

        processed_files_dict = data_package.process_all_data()

        assert 3 == len(processed_files_dict)

        sample_data = processed_files_dict['SampleA']
        assert isinstance(sample_data, NucleotideSampleData)
        sample_data = cast(NucleotideSampleData, sample_data)
        mask = sample_data.get_mask()
        assert 437 == len(mask)
        assert {'reference'} == mask.sequence_names()
        assert sample_data.get_mask_file().parent == tmp_file
        vcf_file, vcf_index = sample_data.get_vcf_file()
        assert vcf_file.exists()
        assert vcf_index.exists()
        assert vcf_file.parent == tmp_file
        assert vcf_index.parent == tmp_file

        sample_data = processed_files_dict['SampleB']
        assert isinstance(sample_data, NucleotideSampleData)
        sample_data = cast(NucleotideSampleData, sample_data)
        mask = sample_data.get_mask()
        assert 276 == len(mask)
        assert {'reference'} == mask.sequence_names()
        assert sample_data.get_mask_file().parent == tmp_file
        vcf_file, vcf_index = sample_data.get_vcf_file()
        assert vcf_file.exists()
        assert vcf_index.exists()
        assert vcf_file.parent == tmp_file
        assert vcf_index.parent == tmp_file

        sample_data = processed_files_dict['SampleC']
        assert isinstance(sample_data, NucleotideSampleData)
        sample_data = cast(NucleotideSampleData, sample_data)
        mask = sample_data.get_mask()
        assert 329 == len(mask)
        assert {'reference'} == mask.sequence_names()
        assert sample_data.get_mask_file().parent == tmp_file
        vcf_file, vcf_index = sample_data.get_vcf_file()
        assert vcf_file.exists()
        assert vcf_index.exists()
        assert vcf_file.parent == tmp_file
        assert vcf_index.parent == tmp_file


def test_with_null_sample_files_processor(sample_dirs):
    vcf_masks = vcf_and_mask_files(sample_dirs)
    data_package = NucleotideSampleDataPackage.create_from_sequence_masks(sample_vcf_map=vcf_masks['vcfs'],
                                                                          masked_genomic_files_map=vcf_masks[
                                                                              'masks'])

    processed_files_dict = data_package.process_all_data()

    assert 3 == len(processed_files_dict)
    sample_data = processed_files_dict['SampleA']
    assert isinstance(sample_data, NucleotideSampleData)
    sample_data = cast(NucleotideSampleData, sample_data)

    with pytest.raises(Exception) as execinfo:
        sample_data.get_mask()
    assert 'Sample mask file is not preprocessed for sample' in str(execinfo.value)

    with pytest.raises(Exception) as execinfo:
        sample_data.get_vcf_file()
    assert 'VCF file for sample' in str(execinfo.value)
