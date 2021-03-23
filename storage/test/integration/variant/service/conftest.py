import tempfile
import warnings
from pathlib import Path

import pytest

warnings.filterwarnings("ignore", category=DeprecationWarning)

from storage.test.integration.variant import sample_dirs, reference_file
from storage.variant.service import DatabaseConnection
from storage.variant.service.ReferenceService import ReferenceService
from storage.variant.service.VariationService import VariationService
from storage.variant.service.CoreAlignmentService import CoreAlignmentService
from storage.variant.io.SnippyVariantsReader import SnippyVariantsReader
from storage.variant.service.SampleService import SampleService
# from storage.variant.service.TreeService import TreeService
from storage.FilesystemStorage import FilesystemStorage


@pytest.fixture
def database() -> DatabaseConnection:
    return DatabaseConnection('sqlite:///:memory:')


@pytest.fixture
def filesystem_storage() -> FilesystemStorage:
    return FilesystemStorage(Path(tempfile.mkdtemp(prefix='index-test')))


@pytest.fixture
def reference_service(database, filesystem_storage) -> ReferenceService:
    reference_service = ReferenceService(database, filesystem_storage.reference_dir)
    return reference_service


@pytest.fixture
def snippy_variants_reader() -> SnippyVariantsReader:
    return SnippyVariantsReader(sample_dirs)


@pytest.fixture
def reference_service_with_data(reference_service) -> ReferenceService:
    reference_service.add_reference_genome(reference_file)
    return reference_service


@pytest.fixture
def sample_service(database):
    return SampleService(database)


@pytest.fixture
def variation_service(database, reference_service_with_data,
                      snippy_variants_reader, sample_service, filesystem_storage) -> VariationService:
    var_service = VariationService(database_connection=database,
                                   reference_service=reference_service_with_data,
                                   sample_service=sample_service,
                                   variation_dir=filesystem_storage.variation_dir)
    var_service.insert_variants(reference_name='genome',
                                variants_reader=snippy_variants_reader)
    return var_service


@pytest.fixture
def core_alignment_service(database, reference_service_with_data, variation_service,
                           sample_service) -> CoreAlignmentService:
    return CoreAlignmentService(database=database,
                                reference_service=reference_service_with_data,
                                variation_service=variation_service,
                                sample_service=sample_service,
                                )


# @pytest.fixture
# def tree_service(database, reference_service_with_data, core_alignment_service) -> TreeService:
#     return TreeService(database, reference_service_with_data, core_alignment_service)
#
#
# @pytest.fixture
# def tree_service_with_tree_stored(database, reference_service_with_data,
#                                   core_alignment_service, variation_service) -> TreeService:
#     tree_service = TreeService(database, reference_service_with_data, core_alignment_service)
#     tree_service.rebuild_tree('genome',
#                               tree_build_type='iqtree',
#                               extra_params='-m MFP+ASC --seed 42')
#
#     return tree_service
