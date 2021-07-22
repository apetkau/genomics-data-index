import shutil
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Dict, Any

from genomics_data_index.configuration.connector.FilesystemStorage import FilesystemStorage
from genomics_data_index.storage.io.mutation.NucleotideSampleDataPackage import NucleotideSampleDataPackage
from genomics_data_index.storage.model.db import DatabasePathTranslator, SampleNucleotideVariation, Sample, \
    SampleKmerIndex
from genomics_data_index.storage.service import DatabaseConnection
from genomics_data_index.storage.service.KmerService import KmerService
from genomics_data_index.storage.service.ReferenceService import ReferenceService
from genomics_data_index.storage.service.SampleService import SampleService
from genomics_data_index.storage.service.VariationService import VariationService
from genomics_data_index.test.integration import reference_file
from genomics_data_index.test.integration import sourmash_signatures


def setup_services(root_dir: Path, database_file: Path) -> Dict[str, Any]:
    filesystem_storage = FilesystemStorage(root_dir)
    dpt = DatabasePathTranslator(filesystem_storage.root_dir)
    db_connection = DatabaseConnection(f'sqlite:///{database_file}', dpt)
    reference_service = ReferenceService(db_connection, filesystem_storage.reference_dir)
    sample_service = SampleService(db_connection)
    variation_service = VariationService(database_connection=db_connection,
                                         reference_service=reference_service,
                                         sample_service=sample_service,
                                         variation_dir=filesystem_storage.variation_dir,
                                         index_unknown_missing=False,
                                         sql_select_limit=50)
    kmer_service = KmerService(database_connection=db_connection,
                               sample_service=sample_service,
                               features_dir=filesystem_storage.kmer_dir)

    return {
        'reference_service': reference_service,
        'variation_service': variation_service,
        'database_connection': db_connection,
        'kmer_service': kmer_service
    }


def test_path_translation_nucleotide_data(snippy_nucleotide_data_package: NucleotideSampleDataPackage):
    with TemporaryDirectory() as root_dir_1_str:
        root_dir_1 = Path(root_dir_1_str)
        data_dir_1 = root_dir_1 / 'data'
        database_file = root_dir_1 / 'db.sqlite'
        services = setup_services(data_dir_1, database_file)
        variation_service = services['variation_service']
        reference_service = services['reference_service']
        database_connection = services['database_connection']

        # Insert data
        reference_service.add_reference_genome(reference_file)
        variation_service.insert(feature_scope_name='genome',
                                 data_package=snippy_nucleotide_data_package)

        # Get a single saved object which has attached files on the filesystem
        sampleA_nv = database_connection.get_session().query(SampleNucleotideVariation) \
            .join(SampleNucleotideVariation.sample) \
            .filter(Sample.name == 'SampleA').one()

        # Make sure this file exists and is in the correct location
        nvf: Path = sampleA_nv.nucleotide_variants_file
        bf: Path = sampleA_nv.masked_regions_file
        assert nvf.exists()
        assert nvf.parent.parent == data_dir_1
        assert bf.exists()
        assert bf.parent.parent == data_dir_1

        # Now we move data files into brand new directory
        with TemporaryDirectory() as root_dir_2_str:
            root_dir_2 = Path(root_dir_2_str)
            data_dir_2 = root_dir_2 / 'data'
            shutil.move(data_dir_1, data_dir_2)

            assert not nvf.exists(), 'Data should be moved so original file not exist'
            assert not bf.exists(), 'Data should be moved so original file not exist'

            services_2 = setup_services(data_dir_2, database_file)
            database_connection_2 = services_2['database_connection']

            # Get a single saved object which has attached files on the filesystem
            sampleA_nv_2 = database_connection_2.get_session().query(SampleNucleotideVariation) \
                .join(SampleNucleotideVariation.sample) \
                .filter(Sample.name == 'SampleA').one()
            nvf_2: Path = sampleA_nv_2.nucleotide_variants_file
            bf_2: Path = sampleA_nv_2.masked_regions_file

            assert nvf_2.exists(), 'Path from database should now correspond to moved path'
            assert nvf_2.parent.parent == data_dir_2
            assert bf_2.exists()
            assert bf_2.parent.parent == data_dir_2


def test_path_translation_kmer_data(snippy_nucleotide_data_package: NucleotideSampleDataPackage):
    with TemporaryDirectory() as root_dir_1_str:
        root_dir_1 = Path(root_dir_1_str)
        data_dir_1 = root_dir_1 / 'data'
        database_file = root_dir_1 / 'db.sqlite'

        services = setup_services(data_dir_1, database_file)
        kmer_service = services['kmer_service']
        database_connection = services['database_connection']

        # Insert data
        for sample_name in sourmash_signatures:
            kmer_service.insert_kmer_index(sample_name=sample_name,
                                           kmer_index_path=sourmash_signatures[sample_name])

        # Get a single saved object which has attached files on the filesystem
        sampleA_kmer = database_connection.get_session().query(SampleKmerIndex) \
            .join(SampleKmerIndex.sample) \
            .filter(Sample.name == 'SampleA').one()

        # Make sure this file exists and is in the correct location
        kf: Path = sampleA_kmer.kmer_index_path
        assert kf.exists()
        assert kf.parent.parent == data_dir_1

        # Now we move data files into brand new directory
        with TemporaryDirectory() as root_dir_2_str:
            root_dir_2 = Path(root_dir_2_str)
            data_dir_2 = root_dir_2 / 'data'
            shutil.move(data_dir_1, data_dir_2)

            assert not kf.exists(), 'Data should be moved so original file not exist'

            services_2 = setup_services(data_dir_2, database_file)
            database_connection_2 = services_2['database_connection']

            # Get a single saved object which has attached files on the filesystem
            sampleA_kmer_2 = database_connection_2.get_session().query(SampleKmerIndex) \
                .join(SampleKmerIndex.sample) \
                .filter(Sample.name == 'SampleA').one()
            kf_2: Path = sampleA_kmer_2.kmer_index_path

            assert kf_2.exists(), 'Path from database should now correspond to moved path'
            assert kf_2.parent.parent == data_dir_2
