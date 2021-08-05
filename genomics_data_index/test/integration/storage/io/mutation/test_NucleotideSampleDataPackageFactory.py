from typing import List
from tempfile import TemporaryDirectory
from pathlib import Path

import pandas as pd

from genomics_data_index.test.integration import data_dir as snippy_dir
from genomics_data_index.test.integration.storage.io.mutation import vcf_and_mask_files
from genomics_data_index.storage.io.mutation.NucleotideSampleDataPackageFactory import \
    NucleotideSnippySampleDataPackageFactory

from genomics_data_index.storage.io.mutation.NucleotideSampleDataPackageFactory import \
    NucleotideInputFilesSampleDataPackageFactory


def create_samples_input_file(tmp_dir: Path, sample_dirs: List[Path]) -> Path:
    file = tmp_dir / 'input-files.tsv'

    vcf_mask_files = vcf_and_mask_files(sample_dirs)
    data = []
    vcfs = vcf_mask_files['vcfs']
    masks = vcf_mask_files['masks']
    for sample in vcfs:
        vcf = vcfs[sample]
        mask = masks[sample]
        data.append([sample, str(vcf.absolute()), str(mask.absolute())])

    input_files_df = pd.DataFrame(data,
                                  columns=['Sample', 'VCF', 'Mask File'])

    input_files_df.to_csv(file, sep='\t', index=False)
    return file


def test_snippy_data_package():
    # Test no index unknown
    with TemporaryDirectory() as preprocess_dir:
        preprocess_dir = Path(preprocess_dir)
        data_package_factory = NucleotideSnippySampleDataPackageFactory(ncores=1, index_unknown=False,
                                                                        preprocess_dir=preprocess_dir,
                                                                        snippy_dir=snippy_dir)
        assert 3 == data_package_factory.number_samples()

        data_package = data_package_factory.create_data_package()
        assert {'SampleA', 'SampleB', 'SampleC'} == data_package.sample_names()

        # Preprocess data
        processed_data_package = data_package.process_all_data()
        assert {'SampleA', 'SampleB', 'SampleC'} == processed_data_package.sample_names()

        features_df = processed_data_package.get_features_reader().get_features_table()
        assert 129 == len(features_df)
        assert {'SampleA', 'SampleB', 'SampleC'} == set(features_df['SAMPLE'].tolist())

        assert 46 == len(features_df[features_df['SAMPLE'] == 'SampleA'])
        assert 50 == len(features_df[features_df['SAMPLE'] == 'SampleB'])
        assert 33 == len(features_df[features_df['SAMPLE'] == 'SampleC'])

    # Test with index unknown
    with TemporaryDirectory() as preprocess_dir:
        preprocess_dir = Path(preprocess_dir)
        data_package_factory = NucleotideSnippySampleDataPackageFactory(ncores=1, index_unknown=True,
                                                                        preprocess_dir=preprocess_dir,
                                                                        snippy_dir=snippy_dir)
        assert 3 == data_package_factory.number_samples()

        data_package = data_package_factory.create_data_package()
        assert {'SampleA', 'SampleB', 'SampleC'} == data_package.sample_names()

        # Preprocess data
        processed_data_package = data_package.process_all_data()
        assert {'SampleA', 'SampleB', 'SampleC'} == processed_data_package.sample_names()

        features_df = processed_data_package.get_features_reader().get_features_table()
        assert 1170 == len(features_df)
        assert {'SampleA', 'SampleB', 'SampleC'} == set(features_df['SAMPLE'].tolist())

        assert 482 == len(features_df[features_df['SAMPLE'] == 'SampleA'])
        assert 326 == len(features_df[features_df['SAMPLE'] == 'SampleB'])
        assert 362 == len(features_df[features_df['SAMPLE'] == 'SampleC'])


def test_vcf_mask_files_data_package(sample_dirs: List[Path]):
    # Test no index unknown
    with TemporaryDirectory() as preprocess_dir:
        preprocess_dir = Path(preprocess_dir)
        input_file = create_samples_input_file(tmp_dir=preprocess_dir, sample_dirs=sample_dirs)
        data_package_factory = NucleotideInputFilesSampleDataPackageFactory(ncores=1, index_unknown=False,
                                                                            preprocess_dir=preprocess_dir,
                                                                            input_files_file=input_file)
        assert 3 == data_package_factory.number_samples()

        data_package = data_package_factory.create_data_package()
        assert {'SampleA', 'SampleB', 'SampleC'} == data_package.sample_names()

        # Preprocess data
        processed_data_package = data_package.process_all_data()
        assert {'SampleA', 'SampleB', 'SampleC'} == processed_data_package.sample_names()

        features_df = processed_data_package.get_features_reader().get_features_table()
        assert 129 == len(features_df)
        assert {'SampleA', 'SampleB', 'SampleC'} == set(features_df['SAMPLE'].tolist())

        assert 46 == len(features_df[features_df['SAMPLE'] == 'SampleA'])
        assert 50 == len(features_df[features_df['SAMPLE'] == 'SampleB'])
        assert 33 == len(features_df[features_df['SAMPLE'] == 'SampleC'])

    # Test with index unknown
    with TemporaryDirectory() as preprocess_dir:
        preprocess_dir = Path(preprocess_dir)
        input_file = create_samples_input_file(tmp_dir=preprocess_dir, sample_dirs=sample_dirs)
        data_package_factory = NucleotideInputFilesSampleDataPackageFactory(ncores=1, index_unknown=True,
                                                                            preprocess_dir=preprocess_dir,
                                                                            input_files_file=input_file)
        assert 3 == data_package_factory.number_samples()

        data_package = data_package_factory.create_data_package()
        assert {'SampleA', 'SampleB', 'SampleC'} == data_package.sample_names()

        # Preprocess data
        processed_data_package = data_package.process_all_data()
        assert {'SampleA', 'SampleB', 'SampleC'} == processed_data_package.sample_names()

        features_df = processed_data_package.get_features_reader().get_features_table()
        assert 1170 == len(features_df)
        assert {'SampleA', 'SampleB', 'SampleC'} == set(features_df['SAMPLE'].tolist())

        assert 482 == len(features_df[features_df['SAMPLE'] == 'SampleA'])
        assert 326 == len(features_df[features_df['SAMPLE'] == 'SampleB'])
        assert 362 == len(features_df[features_df['SAMPLE'] == 'SampleC'])
