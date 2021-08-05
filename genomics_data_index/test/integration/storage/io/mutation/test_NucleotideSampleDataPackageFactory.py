from pathlib import Path
from tempfile import TemporaryDirectory
from typing import List

import pandas as pd

from genomics_data_index.storage.io.mutation.NucleotideSampleDataPackageFactory import \
    NucleotideInputFilesSampleDataPackageFactory
from genomics_data_index.storage.io.mutation.NucleotideSampleDataPackageFactory import \
    NucleotideSnippySampleDataPackageFactory
from genomics_data_index.test.integration import data_dir as snippy_dir
from genomics_data_index.test.integration.storage.io.mutation import vcf_and_mask_files


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


def verify_features_df_for_sample(sample_name: str, features_df: pd.DataFrame) -> None:
    if sample_name == 'SampleA':
        assert 482 == len(features_df)
        assert {'SampleA'} == set(features_df['SAMPLE'].tolist())
    elif sample_name == 'SampleB':
        assert 326 == len(features_df)
        assert {'SampleB'} == set(features_df['SAMPLE'].tolist())
    elif sample_name == 'SampleC':
        assert 362 == len(features_df)
        assert {'SampleC'} == set(features_df['SAMPLE'].tolist())
    else:
        raise Exception(f'sample={sample_name} is not expected')


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


def test_snippy_data_package_iter():
    with TemporaryDirectory() as preprocess_dir:
        preprocess_dir = Path(preprocess_dir)
        data_package_factory = NucleotideSnippySampleDataPackageFactory(ncores=1, index_unknown=True,
                                                                        preprocess_dir=preprocess_dir,
                                                                        snippy_dir=snippy_dir)
        assert 3 == data_package_factory.number_samples()

        number_iterations = 0
        handled_samples = {}
        for data_package in data_package_factory.create_data_package_iter(1):
            assert 1 == len(data_package.sample_names())
            processed_data_package = data_package.process_all_data()
            sample_name = list(processed_data_package.sample_names())[0]
            features_df = processed_data_package.get_features_reader().get_features_table()

            verify_features_df_for_sample(sample_name=sample_name, features_df=features_df)

            handled_samples[sample_name] = True
            number_iterations = number_iterations + 1
        assert 3 == number_iterations
        assert handled_samples['SampleA']
        assert handled_samples['SampleB']
        assert handled_samples['SampleC']


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


def test_vcf_mask_files_data_package_iter(sample_dirs: List[Path]):
    # Test batch size 1
    with TemporaryDirectory() as preprocess_dir:
        preprocess_dir = Path(preprocess_dir)
        input_file = create_samples_input_file(tmp_dir=preprocess_dir, sample_dirs=sample_dirs)
        data_package_factory = NucleotideInputFilesSampleDataPackageFactory(ncores=1, index_unknown=True,
                                                                            preprocess_dir=preprocess_dir,
                                                                            input_files_file=input_file)
        assert 3 == data_package_factory.number_samples()

        number_iterations = 0
        handled_samples = {}
        for data_package in data_package_factory.create_data_package_iter(1):
            assert 1 == len(data_package.sample_names())
            processed_data_package = data_package.process_all_data()
            sample_name = list(processed_data_package.sample_names())[0]
            features_df = processed_data_package.get_features_reader().get_features_table()

            verify_features_df_for_sample(sample_name=sample_name, features_df=features_df)

            handled_samples[sample_name] = True
            number_iterations = number_iterations + 1
        assert 3 == number_iterations
        assert handled_samples['SampleA']
        assert handled_samples['SampleB']
        assert handled_samples['SampleC']

    # Test batch size 2
    with TemporaryDirectory() as preprocess_dir:
        preprocess_dir = Path(preprocess_dir)
        input_file = create_samples_input_file(tmp_dir=preprocess_dir, sample_dirs=sample_dirs)
        data_package_factory = NucleotideInputFilesSampleDataPackageFactory(ncores=1, index_unknown=True,
                                                                            preprocess_dir=preprocess_dir,
                                                                            input_files_file=input_file)
        assert 3 == data_package_factory.number_samples()

        handled_samples = {}
        data_package_iter = data_package_factory.create_data_package_iter(2)

        # Iteration 1
        data_package = next(data_package_iter)

        assert 2 == len(data_package.sample_names())
        processed_data_package = data_package.process_all_data()
        features_df = processed_data_package.get_features_reader().get_features_table()

        sample_names = list(processed_data_package.sample_names())
        for sample_name in sample_names:
            features_df_sample = features_df[features_df['SAMPLE'] == sample_name]
            verify_features_df_for_sample(sample_name=sample_name, features_df=features_df_sample)
            handled_samples[sample_name] = True

        # Iteration 2
        data_package = next(data_package_iter)

        assert 1 == len(data_package.sample_names())
        processed_data_package = data_package.process_all_data()
        features_df = processed_data_package.get_features_reader().get_features_table()

        sample_names = list(processed_data_package.sample_names())
        assert 1 == len(sample_names)
        sample_name = sample_names[0]
        verify_features_df_for_sample(sample_name=sample_name, features_df=features_df)
        handled_samples[sample_name] = True

        assert handled_samples['SampleA']
        assert handled_samples['SampleB']
        assert handled_samples['SampleC']

        # Make sure no data packages left to iterate over
        assert next(data_package_iter, None) is None

    # Test batch size 3
    with TemporaryDirectory() as preprocess_dir:
        preprocess_dir = Path(preprocess_dir)
        input_file = create_samples_input_file(tmp_dir=preprocess_dir, sample_dirs=sample_dirs)
        data_package_factory = NucleotideInputFilesSampleDataPackageFactory(ncores=1, index_unknown=True,
                                                                            preprocess_dir=preprocess_dir,
                                                                            input_files_file=input_file)
        assert 3 == data_package_factory.number_samples()

        handled_samples = {}
        data_package_iter = data_package_factory.create_data_package_iter(3)

        # Iteration 1
        data_package = next(data_package_iter)

        assert 3 == len(data_package.sample_names())
        processed_data_package = data_package.process_all_data()
        features_df = processed_data_package.get_features_reader().get_features_table()

        sample_names = list(processed_data_package.sample_names())
        for sample_name in sample_names:
            features_df_sample = features_df[features_df['SAMPLE'] == sample_name]
            verify_features_df_for_sample(sample_name=sample_name, features_df=features_df_sample)
            handled_samples[sample_name] = True

        assert handled_samples['SampleA']
        assert handled_samples['SampleB']
        assert handled_samples['SampleC']

        # Make sure no data packages left to iterate over
        assert next(data_package_iter, None) is None
