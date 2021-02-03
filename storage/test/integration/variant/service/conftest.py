from typing import Dict
import pytest
import tempfile
from pathlib import Path

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

from storage.test.integration.variant import sample_dirs, reference_file
from storage.variant.service import DatabaseConnection
from storage.variant.service.ReferenceService import ReferenceService
from storage.variant.service.VariationService import VariationService
from storage.variant.VariantsReader import SnippyVariantsReader
from storage.variant.model import ReferenceSequence


@pytest.fixture
def database() -> DatabaseConnection:
    return DatabaseConnection('sqlite:///:memory:')


@pytest.fixture
def reference_service(database) -> ReferenceService:
    seq_repo_root = Path(tempfile.mkdtemp(prefix='index-test'))
    reference_service = ReferenceService(database, seq_repo_root)

    reference_service.add_reference_genome(reference_file)

    return reference_service


@pytest.fixture
def snippy_variants_reader() -> SnippyVariantsReader:
    return SnippyVariantsReader(sample_dirs)


@pytest.fixture
def ref_contigs(reference_service) -> Dict[str, ReferenceSequence]:
    reference_service.create_reference_genome(reference_file)
    reference = reference_service.find_reference_genome('genome')

    return {s.sequence_name: s for s in reference.sequences}


@pytest.fixture
def variation_service(database, reference_service,
                      snippy_variants_reader, ref_contigs) -> VariationService:
    var_df = snippy_variants_reader.get_variants_table()
    core_masks = snippy_variants_reader.get_core_masks()

    var_service = VariationService(database)
    var_service.insert_variants(var_df=var_df,
                                ref_contigs=ref_contigs,
                                core_masks=core_masks)

    return var_service
