import tempfile
import warnings
from pathlib import Path

import pytest

warnings.filterwarnings("ignore", category=DeprecationWarning)

from storage.test.integration.variant import sample_dirs, reference_file, regular_vcf_dir, data_dir
from storage.test.integration.variant import mlst_file_single_scheme, basic_mlst_file, mlst_file_unknown
from storage.test.integration.variant import sourmash_signatures
from storage.variant.service import DatabaseConnection
from storage.variant.service.ReferenceService import ReferenceService
from storage.variant.service.VariationService import VariationService
from storage.variant.service.CoreAlignmentService import CoreAlignmentService
from storage.variant.io.mutation.SnippyVariantsReader import SnippyVariantsReader
from storage.variant.io.mutation.VcfVariantsReader import VcfVariantsReader
from storage.variant.io.mlst.MLSTFeaturesReader import MLSTFeaturesReader
from storage.variant.io.mlst.MLSTTSeemannFeaturesReader import MLSTTSeemannFeaturesReader
from storage.variant.service.SampleService import SampleService
from storage.variant.service.TreeService import TreeService
from storage.variant.service.KmerQueryService import KmerQueryService
from storage.variant.service.KmerService import KmerService
from storage.FilesystemStorage import FilesystemStorage
from storage.variant.service.MLSTService import MLSTService


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
def regular_variants_reader() -> VcfVariantsReader:
    vcf_files = {
        'SampleA': Path(regular_vcf_dir, 'SampleA.vcf.gz'),
        'SampleB': Path(regular_vcf_dir, 'SampleB.vcf.gz'),
        'SampleC': Path(regular_vcf_dir, 'SampleC.vcf.gz'),
    }

    mask_files = {
        'SampleA': Path(data_dir, 'SampleA', 'snps.aligned.fa'),
        'SampleB': Path(data_dir, 'SampleB', 'snps.aligned.fa'),
        'SampleC': Path(data_dir, 'SampleC', 'snps.aligned.fa'),
    }

    return VcfVariantsReader(sample_vcf_map=vcf_files, masked_genomic_files_map=mask_files)


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
    var_service.insert(feature_scope_name='genome',
                       features_reader=snippy_variants_reader)
    return var_service


@pytest.fixture
def variation_service_non_snippy_vcfs(database, reference_service_with_data,
                                      regular_variants_reader, sample_service, filesystem_storage) -> VariationService:
    var_service = VariationService(database_connection=database,
                                   reference_service=reference_service_with_data,
                                   sample_service=sample_service,
                                   variation_dir=filesystem_storage.variation_dir)
    var_service.insert(feature_scope_name='genome',
                       features_reader=regular_variants_reader)
    return var_service


@pytest.fixture
def core_alignment_service(database, reference_service_with_data, variation_service,
                           sample_service) -> CoreAlignmentService:
    return CoreAlignmentService(database=database,
                                reference_service=reference_service_with_data,
                                variation_service=variation_service,
                                sample_service=sample_service,
                                )


@pytest.fixture
def core_alignment_service_non_snippy_vcfs(database, reference_service_with_data, variation_service_non_snippy_vcfs,
                                           sample_service) -> CoreAlignmentService:
    return CoreAlignmentService(database=database,
                                reference_service=reference_service_with_data,
                                variation_service=variation_service_non_snippy_vcfs,
                                sample_service=sample_service,
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


@pytest.fixture
def kmer_service_with_data(database, sample_service) -> KmerService:
    kmer_service = KmerService(database_connection=database,
                               sample_service=sample_service)
    for sample_name in sourmash_signatures:
        kmer_service.insert_kmer_index(sample_name=sample_name,
                                       kmer_index_path=sourmash_signatures[sample_name])

    return kmer_service


@pytest.fixture
def kmer_query_service_with_data(sample_service, kmer_service_with_data) -> KmerQueryService:
    return KmerQueryService(sample_service=sample_service)


@pytest.fixture
def mlst_reader_single_scheme() -> MLSTFeaturesReader:
    return MLSTTSeemannFeaturesReader(mlst_file=mlst_file_single_scheme)


@pytest.fixture
def mlst_reader_basic() -> MLSTFeaturesReader:
    return MLSTTSeemannFeaturesReader(mlst_file=basic_mlst_file)


@pytest.fixture
def mlst_reader_unknown() -> MLSTFeaturesReader:
    return MLSTTSeemannFeaturesReader(mlst_file=mlst_file_unknown)


@pytest.fixture
def mlst_service_loaded(mlst_reader_basic, database, sample_service, filesystem_storage) -> MLSTService:
    mlst_service = MLSTService(database_connection=database,
                               sample_service=sample_service,
                               mlst_dir=filesystem_storage.mlst_dir)
    mlst_service.insert(features_reader=mlst_reader_basic)

    return mlst_service


@pytest.fixture
def mlst_service_loaded_unknown(mlst_reader_unknown, database, sample_service, filesystem_storage) -> MLSTService:
    mlst_service = MLSTService(database_connection=database,
                               sample_service=sample_service,
                               mlst_dir=filesystem_storage.mlst_dir)
    mlst_service.insert(features_reader=mlst_reader_unknown)

    return mlst_service
