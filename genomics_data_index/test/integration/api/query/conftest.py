import tempfile
import warnings
from pathlib import Path

import pytest
from ete3 import Tree

warnings.filterwarnings("ignore", category=DeprecationWarning)

from genomics_data_index.configuration.Project import Project
from genomics_data_index.storage.io.mlst.MLSTSampleDataPackage import MLSTSampleDataPackage
from genomics_data_index.storage.io.mlst.MLSTTSeemannFeaturesReader import MLSTTSeemannFeaturesReader
from genomics_data_index.api.query.GenomicsDataIndex import GenomicsDataIndex

from genomics_data_index.test.integration import sample_dirs, reference_file, basic_mlst_file, tree_file, \
    mlst_snippy_file
from genomics_data_index.test.integration import sourmash_signatures, snpeff_sample_vcfs, snpeff_sample_vcfs_fake_dup, \
    reference_file_snpeff
from genomics_data_index.test.integration import snippy_sample2_vcfs_dict, snippy_sample2_mask_sequences_dict
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.storage.io.mutation.NucleotideSampleDataPackage import NucleotideSampleDataPackage
from genomics_data_index.storage.io.processor.SerialSampleFilesProcessor import SerialSampleFilesProcessor
from genomics_data_index.storage.io.mutation.variants_processor.MultipleProcessVcfVariantsTableProcessor import \
    MultipleProcessVcfVariantsTableProcessorFactory


@pytest.fixture
def snippy_data_package_2() -> NucleotideSampleDataPackage:
    tmp_dir = Path(tempfile.mkdtemp())
    data_package = NucleotideSampleDataPackage.create_from_sequence_masks(sample_vcf_map=snippy_sample2_vcfs_dict,
                                                                          masked_genomic_files_map=snippy_sample2_mask_sequences_dict,
                                                                          sample_files_processor=SerialSampleFilesProcessor(
                                                                              tmp_dir))
    return data_package


@pytest.fixture
def loaded_database_connection() -> DataIndexConnection:
    snippy_tmp_dir = Path(tempfile.mkdtemp())
    data_package = NucleotideSampleDataPackage.create_from_snippy(sample_dirs,
                                                                  sample_files_processor=SerialSampleFilesProcessor(
                                                                      snippy_tmp_dir))
    return build_database_connection(data_package=data_package)


@pytest.fixture
def loaded_database_connection_parallel_variants() -> DataIndexConnection:
    snippy_tmp_dir = Path(tempfile.mkdtemp())
    variants_processor_factory = MultipleProcessVcfVariantsTableProcessorFactory(ncores=2)
    data_package = NucleotideSampleDataPackage.create_from_snippy(sample_dirs,
                                                                  sample_files_processor=SerialSampleFilesProcessor(
                                                                      snippy_tmp_dir),
                                                                  variants_processor_factory=variants_processor_factory)
    return build_database_connection(data_package=data_package)


def build_database_connection(data_package: NucleotideSampleDataPackage) -> DataIndexConnection:
    tmp_dir = Path(tempfile.mkdtemp())
    database_connection = DataIndexConnection.connect(database_connection='sqlite:///:memory:',
                                                      database_dir=tmp_dir)

    # Load Nucleotide variation
    database_connection.reference_service.add_reference_genome(reference_file)
    database_connection.variation_service.insert(data_package, feature_scope_name='genome')

    # Load MLST
    mlst_package = MLSTSampleDataPackage(MLSTTSeemannFeaturesReader(mlst_file=basic_mlst_file))
    database_connection.mlst_service.insert(mlst_package)

    # Load kmers
    kmer_service = database_connection.kmer_service
    for sample_name in sourmash_signatures:
        kmer_service.insert_kmer_index(sample_name=sample_name,
                                       kmer_index_path=sourmash_signatures[sample_name])

    # Load MLST results with overlapping samples with Nucleotide variation
    mlst_package_snippy = MLSTSampleDataPackage(MLSTTSeemannFeaturesReader(mlst_file=mlst_snippy_file))
    database_connection.mlst_service.insert(mlst_package_snippy)

    return database_connection


@pytest.fixture
def loaded_database_connection_annotations() -> DataIndexConnection:
    tmp_dir = Path(tempfile.mkdtemp())
    database_connection = DataIndexConnection.connect(database_connection='sqlite:///:memory:',
                                                      database_dir=tmp_dir)

    # Load Nucleotide variation
    database_connection.reference_service.add_reference_genome(reference_file_snpeff)
    vcf_tmp_dir = Path(tempfile.mkdtemp())
    data_package = NucleotideSampleDataPackage.create_from_sequence_masks(sample_vcf_map=snpeff_sample_vcfs,
                                                                          masked_genomic_files_map=None,
                                                                          sample_files_processor=SerialSampleFilesProcessor(
                                                                              vcf_tmp_dir))
    database_connection.variation_service.insert(data_package, feature_scope_name='NC_011083')

    return database_connection


@pytest.fixture
def loaded_database_connection_annotations_duplicate_genes() -> DataIndexConnection:
    tmp_dir = Path(tempfile.mkdtemp())
    database_connection = DataIndexConnection.connect(database_connection='sqlite:///:memory:',
                                                      database_dir=tmp_dir)

    # Load Nucleotide variation
    database_connection.reference_service.add_reference_genome(reference_file_snpeff)
    vcf_tmp_dir = Path(tempfile.mkdtemp())
    data_package = NucleotideSampleDataPackage.create_from_sequence_masks(sample_vcf_map=snpeff_sample_vcfs_fake_dup,
                                                                          masked_genomic_files_map=None,
                                                                          sample_files_processor=SerialSampleFilesProcessor(
                                                                              vcf_tmp_dir))
    database_connection.variation_service.insert(data_package, feature_scope_name='NC_011083')

    return database_connection


@pytest.fixture
def loaded_database_genomic_data_store(loaded_database_connection):
    return GenomicsDataIndex(connection=loaded_database_connection)


@pytest.fixture
def loaded_database_genomic_data_store_annotations(loaded_database_connection_annotations):
    return GenomicsDataIndex(connection=loaded_database_connection_annotations)


@pytest.fixture
def loaded_database_genomic_data_store_with_tree(loaded_database_connection_with_built_tree):
    return GenomicsDataIndex(connection=loaded_database_connection_with_built_tree)


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
def loaded_data_store_from_project_dir() -> GenomicsDataIndex:
    tmp_dir = Path(tempfile.mkdtemp())
    project_dir = tmp_dir / 'project'
    project = Project.initialize_project(project_dir)

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

    return GenomicsDataIndex.connect(project_dir=project_dir)


@pytest.fixture
def prebuilt_tree() -> Tree:
    return Tree(str(tree_file))
