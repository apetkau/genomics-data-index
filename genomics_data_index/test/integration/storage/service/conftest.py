import tempfile
import warnings
from pathlib import Path

import pytest

warnings.filterwarnings("ignore", category=DeprecationWarning)

from genomics_data_index.test.integration import sample_dirs, reference_file, regular_vcf_dir, data_dir
from genomics_data_index.test.integration import sample_dirs_AB, sample_dirs_C
from genomics_data_index.test.integration import mlst_file_single_scheme, basic_mlst_file, mlst_file_unknown
from genomics_data_index.test.integration import mlst_file_single_scheme2, mlst_file_single_scheme3
from genomics_data_index.test.integration import sourmash_signatures, snpeff_sample_vcfs, reference_file_snpeff
from genomics_data_index.test.integration import snpeff_sample_vcfs_fake_dup
from genomics_data_index.storage.service import DatabaseConnection
from genomics_data_index.storage.service.ReferenceService import ReferenceService
from genomics_data_index.storage.service.VariationService import VariationService
from genomics_data_index.storage.service.CoreAlignmentService import CoreAlignmentService
from genomics_data_index.storage.io.mlst.MLSTFeaturesReader import MLSTFeaturesReader
from genomics_data_index.storage.io.mlst.MLSTTSeemannFeaturesReader import MLSTTSeemannFeaturesReader
from genomics_data_index.storage.service.SampleService import SampleService
from genomics_data_index.storage.service.TreeService import TreeService
from genomics_data_index.storage.service.KmerService import KmerService
from genomics_data_index.configuration.connector.FilesystemStorage import FilesystemStorage
from genomics_data_index.storage.service.MLSTService import MLSTService
from genomics_data_index.storage.io.processor.SerialSampleFilesProcessor import SerialSampleFilesProcessor
from genomics_data_index.storage.io.mutation.NucleotideSampleDataPackage import NucleotideSampleDataPackage
from genomics_data_index.storage.io.mlst.MLSTSampleDataPackage import MLSTSampleDataPackage
from genomics_data_index.storage.model.db.DatabasePathTranslator import DatabasePathTranslator
from genomics_data_index.storage.io.mutation.variants_processor.MultipleProcessVcfVariantsTableProcessor import \
    MultipleProcessVcfVariantsTableProcessorFactory


@pytest.fixture
def filesystem_storage() -> FilesystemStorage:
    return FilesystemStorage(Path(tempfile.mkdtemp(prefix='index-test')))


@pytest.fixture
def database(filesystem_storage) -> DatabaseConnection:
    dpt = DatabasePathTranslator(filesystem_storage.root_dir)
    return DatabaseConnection('sqlite:///:memory:', dpt)


@pytest.fixture
def reference_service(database, filesystem_storage) -> ReferenceService:
    reference_service = ReferenceService(database, filesystem_storage.reference_dir)
    return reference_service


@pytest.fixture
def snippy_nucleotide_data_package() -> NucleotideSampleDataPackage:
    tmp_dir = Path(tempfile.mkdtemp())
    return NucleotideSampleDataPackage.create_from_snippy(sample_dirs,
                                                          sample_files_processor=SerialSampleFilesProcessor(tmp_dir))


@pytest.fixture
def snippy_nucleotide_data_package_no_index_missing() -> NucleotideSampleDataPackage:
    tmp_dir = Path(tempfile.mkdtemp())
    return NucleotideSampleDataPackage.create_from_snippy(sample_dirs,
                                                          sample_files_processor=SerialSampleFilesProcessor(tmp_dir),
                                                          index_unknown_missing=False)


@pytest.fixture
def snippy_nucleotide_data_package_AB() -> NucleotideSampleDataPackage:
    tmp_dir = Path(tempfile.mkdtemp())
    return NucleotideSampleDataPackage.create_from_snippy(sample_dirs_AB,
                                                          sample_files_processor=SerialSampleFilesProcessor(tmp_dir))


@pytest.fixture
def snippy_nucleotide_data_package_C() -> NucleotideSampleDataPackage:
    tmp_dir = Path(tempfile.mkdtemp())
    return NucleotideSampleDataPackage.create_from_snippy(sample_dirs_C,
                                                          sample_files_processor=SerialSampleFilesProcessor(tmp_dir))


@pytest.fixture
def regular_nucleotide_data_package() -> NucleotideSampleDataPackage:
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

    tmp_dir = Path(tempfile.mkdtemp())
    return NucleotideSampleDataPackage.create_from_sequence_masks(sample_vcf_map=vcf_files,
                                                                  masked_genomic_files_map=mask_files,
                                                                  sample_files_processor=SerialSampleFilesProcessor(
                                                                      tmp_dir))


@pytest.fixture
def snpeff_nucleotide_data_package() -> NucleotideSampleDataPackage:
    tmp_dir = Path(tempfile.mkdtemp())
    return NucleotideSampleDataPackage.create_from_sequence_masks(sample_vcf_map=snpeff_sample_vcfs,
                                                                  masked_genomic_files_map=None,
                                                                  sample_files_processor=SerialSampleFilesProcessor(
                                                                      tmp_dir))


@pytest.fixture
def snpeff_nucleotide_data_package_fake_duplicate_gene() -> NucleotideSampleDataPackage:
    tmp_dir = Path(tempfile.mkdtemp())
    return NucleotideSampleDataPackage.create_from_sequence_masks(sample_vcf_map=snpeff_sample_vcfs_fake_dup,
                                                                  masked_genomic_files_map=None,
                                                                  sample_files_processor=SerialSampleFilesProcessor(
                                                                      tmp_dir))


@pytest.fixture
def snpeff_nucleotide_data_package_parallel_variants() -> NucleotideSampleDataPackage:
    tmp_dir = Path(tempfile.mkdtemp())
    variants_processor_factory = MultipleProcessVcfVariantsTableProcessorFactory(ncores=2)
    return NucleotideSampleDataPackage.create_from_sequence_masks(sample_vcf_map=snpeff_sample_vcfs,
                                                                  masked_genomic_files_map=None,
                                                                  sample_files_processor=SerialSampleFilesProcessor(
                                                                      tmp_dir),
                                                                  variants_processor_factory=variants_processor_factory)


@pytest.fixture
def reference_service_with_data(reference_service) -> ReferenceService:
    reference_service.add_reference_genome(reference_file)
    return reference_service


@pytest.fixture
def reference_service_with_snpeff_data(reference_service) -> ReferenceService:
    reference_service.add_reference_genome(reference_file_snpeff)
    return reference_service


@pytest.fixture
def sample_service(database):
    return SampleService(database)


@pytest.fixture
def sample_service_snpeff_annotations(database, reference_service_with_snpeff_data,
                                      snpeff_nucleotide_data_package, sample_service,
                                      filesystem_storage) -> SampleService:
    var_service = VariationService(database_connection=database,
                                   reference_service=reference_service_with_snpeff_data,
                                   sample_service=sample_service,
                                   variation_dir=filesystem_storage.variation_dir,
                                   index_unknown_missing=False,
                                   sql_select_limit=500)
    var_service.insert(feature_scope_name='NC_011083',
                       data_package=snpeff_nucleotide_data_package)
    return sample_service


@pytest.fixture
def sample_service_snpeff_annotations_fake_duplicate_gene(database, reference_service_with_snpeff_data,
                                                          snpeff_nucleotide_data_package_fake_duplicate_gene,
                                                          sample_service,
                                                          filesystem_storage) -> SampleService:
    var_service = VariationService(database_connection=database,
                                   reference_service=reference_service_with_snpeff_data,
                                   sample_service=sample_service,
                                   variation_dir=filesystem_storage.variation_dir,
                                   index_unknown_missing=False,
                                   sql_select_limit=500)
    var_service.insert(feature_scope_name='NC_011083',
                       data_package=snpeff_nucleotide_data_package_fake_duplicate_gene)
    return sample_service


@pytest.fixture
def sample_service_snpeff_annotations_parallel_variants(database, reference_service_with_snpeff_data,
                                                        snpeff_nucleotide_data_package_parallel_variants,
                                                        sample_service,
                                                        filesystem_storage) -> SampleService:
    var_service = VariationService(database_connection=database,
                                   reference_service=reference_service_with_snpeff_data,
                                   sample_service=sample_service,
                                   variation_dir=filesystem_storage.variation_dir,
                                   index_unknown_missing=False,
                                   sql_select_limit=500)
    var_service.insert(feature_scope_name='NC_011083',
                       data_package=snpeff_nucleotide_data_package_parallel_variants)
    return sample_service


@pytest.fixture
def variation_service(database, reference_service_with_data,
                      snippy_nucleotide_data_package, sample_service, filesystem_storage) -> VariationService:
    var_service = VariationService(database_connection=database,
                                   reference_service=reference_service_with_data,
                                   sample_service=sample_service,
                                   variation_dir=filesystem_storage.variation_dir,
                                   index_unknown_missing=False,
                                   sql_select_limit=500)
    var_service.insert(feature_scope_name='genome',
                       data_package=snippy_nucleotide_data_package)
    return var_service


@pytest.fixture
def variation_service_index_unknowns(database, reference_service_with_data,
                                     snippy_nucleotide_data_package, sample_service,
                                     filesystem_storage) -> VariationService:
    var_service = VariationService(database_connection=database,
                                   reference_service=reference_service_with_data,
                                   sample_service=sample_service,
                                   variation_dir=filesystem_storage.variation_dir,
                                   index_unknown_missing=True,
                                   sql_select_limit=500)
    var_service.insert(feature_scope_name='genome',
                       data_package=snippy_nucleotide_data_package)
    return var_service


@pytest.fixture
def variation_service_non_snippy_vcfs(database, reference_service_with_data,
                                      regular_nucleotide_data_package, sample_service,
                                      filesystem_storage) -> VariationService:
    var_service = VariationService(database_connection=database,
                                   reference_service=reference_service_with_data,
                                   sample_service=sample_service,
                                   variation_dir=filesystem_storage.variation_dir,
                                   index_unknown_missing=False,
                                   sql_select_limit=500)
    var_service.insert(feature_scope_name='genome',
                       data_package=regular_nucleotide_data_package)
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
def kmer_service_with_data(database, sample_service, filesystem_storage) -> KmerService:
    kmer_service = KmerService(database_connection=database,
                               sample_service=sample_service,
                               features_dir=filesystem_storage.kmer_dir)
    for sample_name in sourmash_signatures:
        kmer_service.insert_kmer_index(sample_name=sample_name,
                                       kmer_index_path=sourmash_signatures[sample_name])

    return kmer_service


@pytest.fixture
def mlst_reader_single_scheme() -> MLSTFeaturesReader:
    return MLSTTSeemannFeaturesReader(mlst_file=mlst_file_single_scheme)


@pytest.fixture
def mlst_reader_single_scheme2() -> MLSTFeaturesReader:
    return MLSTTSeemannFeaturesReader(mlst_file=mlst_file_single_scheme2)


@pytest.fixture
def mlst_reader_single_scheme3() -> MLSTFeaturesReader:
    return MLSTTSeemannFeaturesReader(mlst_file=mlst_file_single_scheme3)


@pytest.fixture
def mlst_data_package_single_scheme(mlst_reader_single_scheme) -> MLSTSampleDataPackage:
    return MLSTSampleDataPackage(mlst_reader_single_scheme)


@pytest.fixture
def mlst_data_package_single_scheme2(mlst_reader_single_scheme2) -> MLSTSampleDataPackage:
    return MLSTSampleDataPackage(mlst_reader_single_scheme2)


@pytest.fixture
def mlst_data_package_single_scheme3(mlst_reader_single_scheme3) -> MLSTSampleDataPackage:
    return MLSTSampleDataPackage(mlst_reader_single_scheme3)


@pytest.fixture
def mlst_reader_basic() -> MLSTFeaturesReader:
    return MLSTTSeemannFeaturesReader(mlst_file=basic_mlst_file)


@pytest.fixture
def mlst_data_package_basic(mlst_reader_basic) -> MLSTSampleDataPackage:
    return MLSTSampleDataPackage(mlst_reader_basic)


@pytest.fixture
def mlst_reader_unknown() -> MLSTFeaturesReader:
    return MLSTTSeemannFeaturesReader(mlst_file=mlst_file_unknown)


@pytest.fixture
def mlst_data_package_unknown(mlst_reader_unknown) -> MLSTSampleDataPackage:
    return MLSTSampleDataPackage(mlst_reader_unknown)


@pytest.fixture
def mlst_service_loaded(mlst_data_package_basic, database, sample_service, filesystem_storage) -> MLSTService:
    mlst_service = MLSTService(database_connection=database,
                               sample_service=sample_service,
                               mlst_dir=filesystem_storage.mlst_dir)
    mlst_service.insert(data_package=mlst_data_package_basic)

    return mlst_service


@pytest.fixture
def mlst_service_loaded_unknown(mlst_data_package_unknown, database, sample_service, filesystem_storage) -> MLSTService:
    mlst_service = MLSTService(database_connection=database,
                               sample_service=sample_service,
                               mlst_dir=filesystem_storage.mlst_dir)
    mlst_service.insert(data_package=mlst_data_package_unknown)

    return mlst_service
