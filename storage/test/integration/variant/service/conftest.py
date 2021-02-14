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
from storage.variant.service.SampleSequenceService import SampleSequenceService
from storage.variant.service.SampleService import SampleService
from storage.variant.service.TreeService import TreeService


@pytest.fixture
def database() -> DatabaseConnection:
    return DatabaseConnection('sqlite:///:memory:')


@pytest.fixture
def reference_service(database) -> ReferenceService:
    seq_repo_root = Path(tempfile.mkdtemp(prefix='index-test'))
    reference_service = ReferenceService(database, seq_repo_root)
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
                      snippy_variants_reader, sample_service) -> VariationService:
    var_df = snippy_variants_reader.get_variants_table()
    core_masks = snippy_variants_reader.get_core_masks()

    var_service = VariationService(database_connection=database,
                                   reference_service=reference_service_with_data,
                                   sample_service=sample_service)
    var_service.insert_variants(var_df=var_df,
                                reference_name='genome',
                                core_masks=core_masks)
    return var_service


@pytest.fixture
def sample_sequence_service(database, variation_service):
    return SampleSequenceService(database)


@pytest.fixture
def core_alignment_service(database, reference_service_with_data, variation_service,
                           sample_sequence_service) -> CoreAlignmentService:
    return CoreAlignmentService(database=database,
                                reference_service=reference_service_with_data,
                                variation_service=variation_service,
                                sample_sequence_service=sample_sequence_service,
                                )


@pytest.fixture
def tree_service(database, reference_service_with_data, core_alignment_service) -> TreeService:
    return TreeService(database, reference_service_with_data, core_alignment_service)


@pytest.fixture
def tree_service_with_tree_stored(database, reference_service_with_data,
                                  core_alignment_service, variation_service) -> TreeService:
    tree_service = TreeService(database, reference_service_with_data, core_alignment_service)
    tree_service.rebuild_tree('genome',
                              tree_build_type='iqtree',
                              extra_params='-m MFP+ASC --seed 42')

    return tree_service
