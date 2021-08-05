from tempfile import TemporaryDirectory
from pathlib import Path

from genomics_data_index.test.integration import data_dir as snippy_dir
from genomics_data_index.storage.io.mutation.NucleotideSampleDataPackageFactory import NucleotideSnippySampleDataPackageFactory


def test_snippy_data_package():
    # Test no index unknown
    with TemporaryDirectory() as preprocess_dir:
        preprocess_dir = Path(preprocess_dir)
        data_package_factory = NucleotideSnippySampleDataPackageFactory(ncores=1, index_unknown=False,
                                                                        preprocess_dir=preprocess_dir,
                                                                        snippy_dir=snippy_dir)
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
