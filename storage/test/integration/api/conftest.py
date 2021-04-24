import tempfile
import warnings
from pathlib import Path

import pytest

warnings.filterwarnings("ignore", category=DeprecationWarning)

from storage.configuration.Project import Project
from storage.variant.io.mlst.MLSTSampleDataPackage import MLSTSampleDataPackage
from storage.variant.io.mlst.MLSTTSeemannFeaturesReader import MLSTTSeemannFeaturesReader
from storage.api.GenomicDataStore import GenomicDataStore

from storage.test.integration import sample_dirs, reference_file, basic_mlst_file
from storage.configuration.connector.DataIndexConnection import DataIndexConnection
from storage.variant.io.mutation.NucleotideSampleDataPackage import NucleotideSampleDataPackage
from storage.variant.io.processor.SerialSampleFilesProcessor import SerialSampleFilesProcessor


@pytest.fixture
def loaded_database_connection() -> DataIndexConnection:
    tmp_dir = Path(tempfile.mkdtemp())
    database_connection = DataIndexConnection.connect(database_connection='sqlite:///:memory:',
                                                      database_dir=tmp_dir)

    # Load Nucleotide variation
    database_connection.reference_service.add_reference_genome(reference_file)
    snippy_tmp_dir = Path(tempfile.mkdtemp())
    data_package = NucleotideSampleDataPackage.create_from_snippy(sample_dirs,
                                                                  SerialSampleFilesProcessor(snippy_tmp_dir))
    database_connection.variation_service.insert(data_package, feature_scope_name='genome')

    # Load MLST
    mlst_package = MLSTSampleDataPackage(MLSTTSeemannFeaturesReader(mlst_file=basic_mlst_file))
    database_connection.mlst_service.insert(mlst_package)

    # Load MLST results with overlapping samples with Nucleotide variation
    # mlst_package_snippy = MLSTSampleDataPackage(MLSTTSeemannFeaturesReader(mlst_file=mlst_snippy_file))
    # database_connection.mlst_service.insert(mlst_package_snippy)

    return database_connection


@pytest.fixture
def loaded_database_genomic_data_store(loaded_database_connection):
    return GenomicDataStore(connection=loaded_database_connection)


@pytest.fixture
def loaded_database_connection_with_built_tree() -> DataIndexConnection:
    tmp_dir = Path(tempfile.mkdtemp())
    database_connection = DataIndexConnection.connect(database_connection='sqlite:///:memory:',
                                                      database_dir=tmp_dir)

    # Load Nucleotide variation
    database_connection.reference_service.add_reference_genome(reference_file)
    snippy_tmp_dir = Path(tempfile.mkdtemp())
    data_package = NucleotideSampleDataPackage.create_from_snippy(sample_dirs,
                                                                  SerialSampleFilesProcessor(snippy_tmp_dir))
    database_connection.variation_service.insert(data_package, feature_scope_name='genome')
    database_connection.tree_service.rebuild_tree(reference_name='genome', align_type='full')

    # Load MLST
    mlst_package = MLSTSampleDataPackage(MLSTTSeemannFeaturesReader(mlst_file=basic_mlst_file))
    database_connection.mlst_service.insert(mlst_package)

    return database_connection


@pytest.fixture
def loaded_database_only_snippy() -> DataIndexConnection:
    tmp_dir = Path(tempfile.mkdtemp())
    database_connection = DataIndexConnection.connect(database_connection='sqlite:///:memory:',
                                                      database_dir=tmp_dir)

    # Load Nucleotide variation
    database_connection.reference_service.add_reference_genome(reference_file)
    snippy_tmp_dir = Path(tempfile.mkdtemp())
    data_package = NucleotideSampleDataPackage.create_from_snippy(sample_dirs,
                                                                  SerialSampleFilesProcessor(snippy_tmp_dir))
    database_connection.variation_service.insert(data_package, feature_scope_name='genome')

    return database_connection


@pytest.fixture
def loaded_data_store_from_project_dir() -> GenomicDataStore:
    tmp_dir = Path(tempfile.mkdtemp())
    project_dir = tmp_dir / 'project'
    project = Project.create_new_project(project_dir)

    database_connection = project.create_connection()

    # Load Nucleotide variation
    database_connection.reference_service.add_reference_genome(reference_file)
    snippy_tmp_dir = Path(tempfile.mkdtemp())
    data_package = NucleotideSampleDataPackage.create_from_snippy(sample_dirs,
                                                                  SerialSampleFilesProcessor(snippy_tmp_dir))
    database_connection.variation_service.insert(data_package, feature_scope_name='genome')

    # Load MLST
    mlst_package = MLSTSampleDataPackage(MLSTTSeemannFeaturesReader(mlst_file=basic_mlst_file))
    database_connection.mlst_service.insert(mlst_package)
    print(tmp_dir)

    return GenomicDataStore.connect(project_dir=project_dir)
