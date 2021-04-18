from pathlib import Path
from tempfile import TemporaryDirectory
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)

from storage.connector.DataIndexConnection import DataIndexConnection
from storage.test.integration.variant import sample_dirs, reference_file
from storage.variant.io.mutation.NucleotideSampleDataPackage import NucleotideSampleDataPackage
from storage.variant.io.processor.SerialSampleFilesProcessor import SerialSampleFilesProcessor


def test_connection():
    with TemporaryDirectory() as tmp_dir:
        connection = DataIndexConnection.connect(database_connection='sqlite:///:memory:',
                                                 database_dir=Path(tmp_dir))

        assert connection.sample_service is not None
        assert connection.variation_service is not None
        assert connection.reference_service is not None


def test_add_reference_from_connection():
    with TemporaryDirectory() as tmp_dir:
        connection = DataIndexConnection.connect(database_connection='sqlite:///:memory:',
                                                 database_dir=Path(tmp_dir))

        connection.reference_service.add_reference_genome(reference_file)


def test_add_variants_from_connection():
    with TemporaryDirectory() as tmp_dir:
        connection = DataIndexConnection.connect(database_connection='sqlite:///:memory:',
                                                 database_dir=Path(tmp_dir))

        data_package = NucleotideSampleDataPackage.create_from_snippy(sample_dirs,
                                                                      sample_files_processor=SerialSampleFilesProcessor(
                                                                          Path(tmp_dir))
                                                                      )

        connection.reference_service.add_reference_genome(reference_file)
        connection.variation_service.insert(feature_scope_name='genome',
                                            data_package=data_package)
