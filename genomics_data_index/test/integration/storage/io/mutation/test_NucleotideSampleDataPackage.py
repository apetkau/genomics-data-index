import re
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import cast

import pytest

from genomics_data_index.storage.io.mutation.NucleotideSampleData import NucleotideSampleData
from genomics_data_index.storage.io.mutation.NucleotideSampleDataPackage import NucleotideSampleDataPackage
from genomics_data_index.storage.io.mutation.variants_processor.MultipleProcessVcfVariantsTableProcessor import \
    MultipleProcessVcfVariantsTableProcessorFactory
from genomics_data_index.storage.io.mutation.variants_processor.SerialVcfVariantsTableProcessor import \
    SerialVcfVariantsTableProcessorFactory
from genomics_data_index.storage.io.processor.MultipleProcessSampleFilesProcessor import \
    MultipleProcessSampleFilesProcessor
from genomics_data_index.storage.io.processor.SerialSampleFilesProcessor import SerialSampleFilesProcessor
from genomics_data_index.test.integration.storage.io.mutation import vcf_and_mask_files


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


def do_test_get_features_reader_and_features_table(data_package: NucleotideSampleDataPackage):
    # If data is not preprocessed, make sure exception is raised
    features_reader = data_package.get_features_reader()
    assert features_reader is not None
    with pytest.raises(Exception) as execinfo:
        features_reader.get_features_table()
    assert 'VCF file for sample' in str(execinfo.value)
    assert 'is not preprocessed' in str(execinfo.value)

    # Process data
    data_package_processed = data_package.process_all_data()

    # Now get features reader from new data package and verify I can read feature tables
    features_reader = data_package_processed.get_features_reader()
    features_df = features_reader.get_features_table()
    assert 1170 == len(features_df)
    assert {'SampleA', 'SampleB', 'SampleC'} == set(features_df['SAMPLE'].tolist())


def test_get_features_reader_and_features_table_serial_variants_processor(sample_dirs):
    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)
        vcf_masks = vcf_and_mask_files(sample_dirs)
        file_processor = SerialSampleFilesProcessor(tmp_file)
        variants_processor_factory = SerialVcfVariantsTableProcessorFactory.instance()
        data_package = NucleotideSampleDataPackage.create_from_sequence_masks(sample_vcf_map=vcf_masks['vcfs'],
                                                                              masked_genomic_files_map=vcf_masks[
                                                                                  'masks'],
                                                                              variants_processor_factory=variants_processor_factory,
                                                                              sample_files_processor=file_processor)
        do_test_get_features_reader_and_features_table(data_package)


def test_get_features_reader_and_features_table_parallel_variants_processor(sample_dirs):
    with TemporaryDirectory() as tmp_file_str:
        tmp_file = Path(tmp_file_str)
        vcf_masks = vcf_and_mask_files(sample_dirs)
        file_processor = SerialSampleFilesProcessor(tmp_file)
        variants_processor_factory = MultipleProcessVcfVariantsTableProcessorFactory(ncores=2)
        data_package = NucleotideSampleDataPackage.create_from_sequence_masks(sample_vcf_map=vcf_masks['vcfs'],
                                                                              masked_genomic_files_map=vcf_masks[
                                                                                  'masks'],
                                                                              variants_processor_factory=variants_processor_factory,
                                                                              sample_files_processor=file_processor)
        do_test_get_features_reader_and_features_table(data_package)


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

        processed_data_package = cast(NucleotideSampleDataPackage, data_package.process_all_data())
        processed_files_dict = processed_data_package.get_sample_data()
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

        processed_data_package = cast(NucleotideSampleDataPackage, data_package.process_all_data())
        processed_files_dict = processed_data_package.get_sample_data()

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

        processed_data_package = cast(NucleotideSampleDataPackage, data_package.process_all_data())
        processed_files_dict = processed_data_package.get_sample_data()

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

    processed_data_package = cast(NucleotideSampleDataPackage, data_package.process_all_data())
    processed_files_dict = processed_data_package.get_sample_data()

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
