import logging
import multiprocessing
import shutil
import sys
import time
from functools import partial
from os import getcwd, mkdir
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import List, cast

import click
import click_config_file
import coloredlogs
import pandas as pd
from Bio import AlignIO

import genomics_data_index.storage.service.FeatureService as FeatureService
from genomics_data_index import __version__
from genomics_data_index.api.query.GenomicsDataIndex import GenomicsDataIndex
from genomics_data_index.cli import yaml_config_provider
from genomics_data_index.configuration.Project import Project, ProjectConfigurationError
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.pipelines.SnakemakePipelineExecutor import SnakemakePipelineExecutor
from genomics_data_index.storage.index.KmerIndexer import KmerIndexerSourmash, KmerIndexManager
from genomics_data_index.storage.io.SampleDataPackageFactory import SampleDataPackageFactory
from genomics_data_index.storage.io.mlst.MLSTChewbbacaReader import MLSTChewbbacaReader
from genomics_data_index.storage.io.mlst.MLSTSampleDataPackage import MLSTSampleDataPackage
from genomics_data_index.storage.io.mlst.MLSTSistrReader import MLSTSistrReader
from genomics_data_index.storage.io.mlst.MLSTTSeemannFeaturesReader import MLSTTSeemannFeaturesReader
from genomics_data_index.storage.io.mutation.NucleotideSampleDataPackageFactory import \
    NucleotideInputFilesSampleDataPackageFactory
from genomics_data_index.storage.io.mutation.NucleotideSampleDataPackageFactory import \
    NucleotideSnippySampleDataPackageFactory
from genomics_data_index.storage.io.mutation.SequenceFile import SequenceFile
from genomics_data_index.storage.model.QueryFeature import QueryFeature
from genomics_data_index.storage.model.QueryFeatureMLST import QueryFeatureMLST
from genomics_data_index.storage.model.QueryFeatureMutationSPDI import QueryFeatureMutationSPDI
from genomics_data_index.storage.service import EntityExistsError
from genomics_data_index.storage.service.CoreAlignmentService import CoreAlignmentService
from genomics_data_index.storage.service.MLSTService import MLSTService
from genomics_data_index.storage.service.SampleService import SampleService
from genomics_data_index.storage.service.VariationService import VariationService
from genomics_data_index.storage.util import TRACE_LEVEL

logger = logging.getLogger('genomics_data_index')
max_cores = multiprocessing.cpu_count()

LOG_LEVELS = ['TRACE', 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']

click.option = partial(click.option, show_default=True)


def setup_logging(log_level: str) -> None:
    if log_level == 'DEBUG' or log_level == 'TRACE':
        log_format = '%(asctime)s %(levelname)s %(name)s.%(funcName)s,%(lineno)s: %(message)s'

        if log_level == 'TRACE':
            log_level = TRACE_LEVEL
    else:
        log_format = '%(asctime)s %(levelname)s: %(message)s'

    coloredlogs.install(level=log_level, fmt=log_format, logger=logger)


def get_project_exit_on_error(ctx) -> Project:
    project_dir = ctx.obj['project_dir']
    try:
        project = Project(project_dir)
        return project
    except ProjectConfigurationError as e:
        logger.error(f'Could not properly load configuration for project directory [{project_dir}]. '
                     f'Please verify this directory has been configured as a project or set the correct directory '
                     f'with --project-dir')
        logger.debug(e, exc_info=True)
        sys.exit(1)


def get_genomics_index(ctx) -> GenomicsDataIndex:
    project = get_project_exit_on_error(ctx)
    return GenomicsDataIndex.connect(project=project)


@click.group()
@click.pass_context
@click.option('--project-dir', help='A project directory containing the data and connection information.')
@click.option('--ncores', help='Number of cores for any parallel processing', default=max_cores,
              type=click.IntRange(min=1, max=max_cores))
@click.option('--log-level', default='INFO', help='Sets the log level', type=click.Choice(LOG_LEVELS))
@click.version_option(version=__version__)
@click_config_file.configuration_option(provider=yaml_config_provider,
                                        config_file_name='config.yaml',
                                        implicit=False)
def main(ctx, project_dir, ncores, log_level):
    ctx.ensure_object(dict)
    setup_logging(log_level)

    if project_dir is None:
        project_dir = getcwd()

    ctx.obj['project_dir'] = project_dir
    ctx.obj['ncores'] = ncores


@main.command(name='init')
@click.argument('project_dir', type=click.Path(exists=False), nargs=1)
def init(project_dir: str):
    project_dir = Path(project_dir)
    click.echo(f'Initializing empty project in [{project_dir}]')
    Project.initialize_project(project_dir)


@main.group()
@click.pass_context
def load(ctx):
    pass


def load_variants_common(data_index_connection: DataIndexConnection, ncores: int,
                         data_package_factory: SampleDataPackageFactory,
                         reference_file: Path, reference_name: str,
                         input: Path, build_tree: bool,
                         align_type: str, extra_tree_params: str,
                         sample_batch_size: int):
    reference_service = data_index_connection.reference_service
    variation_service = cast(VariationService, data_index_connection.variation_service)
    sample_service = cast(SampleService, data_index_connection.sample_service)
    tree_service = data_index_connection.tree_service

    if reference_name is not None and not reference_service.exists_reference_genome(reference_name):
        logger.error(f'Reference genome [{reference_name}] does not exist in the system. Please try adding '
                     f'a new reference genome by using the --reference-file parameter.')
        sys.exit(1)
    elif reference_name is not None:
        logger.info(f'Using reference genome with name=[{reference_name}]')
    else:
        try:
            logger.info(f'Attempting to load reference genome=[{reference_file}]')
            reference_service.add_reference_genome(reference_file)
        except EntityExistsError as e:
            logger.warning(f'Reference genome [{reference_file}] already exists, will not load')
        reference_name = SequenceFile(reference_file).get_genome_name()

    for data_package in data_package_factory.create_data_package_iter(sample_batch_size):
        samples_exist = sample_service.which_exists(list(data_package.sample_names()))
        if len(samples_exist) > 0:
            max_samples_to_print = 5
            if len(samples_exist) > max_samples_to_print:
                samples_to_print = '[' + ', '.join(samples_exist[0:max_samples_to_print]) + ', ...]'
            else:
                samples_to_print = str(samples_exist)
            logger.error(f'There are {len(samples_exist)} samples which already exist: {samples_to_print}. '
                         f'Will not load any samples.')
        else:
            variation_service.insert(feature_scope_name=reference_name,
                                     data_package=data_package)
        logger.info(f'Loaded variants from [{input}] into database')

        if build_tree:
            tree_service.rebuild_tree(reference_name=reference_name,
                                      align_type=align_type,
                                      num_cores=ncores,
                                      extra_params=extra_tree_params)
            logger.info('Finished building tree of all samples')


@load.command(name='snippy')
@click.pass_context
@click.argument('snippy_dir', type=click.Path(exists=True))
@click.option('--reference-file', help='Reference genome file', required=False, type=click.Path(exists=True))
@click.option('--reference-name', help='Reference genome name', required=False)
@click.option('--index-unknown/--no-index-unknown',
              help='Enable/disable indexing unknown/missing positions. Indexing missing positions can significantly '
                   'slow down the indexing process.',
              required=False, default=True)
@click.option('--sample-batch-size', help='Number of samples to process within a single batch.', default=2000,
              type=click.IntRange(min=1))
@click.option('--build-tree/--no-build-tree', default=False, help='Builds tree of all samples after loading')
@click.option('--align-type', help=f'The type of alignment to generate', default='core',
              type=click.Choice(CoreAlignmentService.ALIGN_TYPES))
@click.option('--extra-tree-params', help='Extra parameters to tree-building software',
              default=None)
def load_snippy(ctx, snippy_dir: str, reference_file: str, reference_name: str,
                index_unknown: bool, sample_batch_size: int, build_tree: bool,
                align_type: str, extra_tree_params: str):
    ncores = ctx.obj['ncores']
    project = get_project_exit_on_error(ctx)
    data_index_connection = project.create_connection()

    if reference_name is None and reference_file is None:
        logger.error(f'Neither --reference-file nor --reference-name are specified. Please define either '
                     f'a new reference genome (--refrence-file [FILE]) or the name of an existing reference genome '
                     f'(--reference-name [NAME]). You can view previously-loaded reference genomes with '
                     f'"gdi --project-dir {project.get_root_dir()} list genomes"')
        sys.exit(1)
    elif reference_file is not None:
        reference_file = Path(reference_file)

    with TemporaryDirectory() as preprocess_dir:
        preprocess_dir = Path(preprocess_dir)
        data_package_factory = NucleotideSnippySampleDataPackageFactory(ncores=ncores, index_unknown=index_unknown,
                                                                        preprocess_dir=preprocess_dir,
                                                                        snippy_dir=Path(snippy_dir))

        load_variants_common(data_index_connection=data_index_connection, ncores=ncores,
                             data_package_factory=data_package_factory,
                             reference_file=reference_file,
                             reference_name=reference_name,
                             input=Path(snippy_dir), build_tree=build_tree, align_type=align_type,
                             extra_tree_params=extra_tree_params,
                             sample_batch_size=sample_batch_size)


@load.command(name='vcf')
@click.pass_context
@click.argument('vcf_fofns', type=click.Path(exists=True))
@click.option('--reference-file', help='Reference genome file', required=False, type=click.Path(exists=True))
@click.option('--reference-name', help='Reference genome name', required=False)
@click.option('--index-unknown/--no-index-unknown',
              help='Enable/disable indexing unknown/missing positions. Indexing missing positions can significantly '
                   'slow down the indexing process.',
              required=False, default=True)
@click.option('--sample-batch-size', help='Number of samples to process within a single batch.', default=2000,
              type=click.IntRange(min=1))
@click.option('--build-tree/--no-build-tree', default=False, help='Builds tree of all samples after loading')
@click.option('--align-type', help=f'The type of alignment to generate', default='core',
              type=click.Choice(CoreAlignmentService.ALIGN_TYPES))
@click.option('--extra-tree-params', help='Extra parameters to tree-building software',
              default=None)
def load_vcf(ctx, vcf_fofns: str, reference_file: str, reference_name: str,
             index_unknown: bool, sample_batch_size: int, build_tree: bool, align_type: str, extra_tree_params: str):
    ncores = ctx.obj['ncores']
    project = get_project_exit_on_error(ctx)
    vcf_fofns = Path(vcf_fofns)

    if reference_name is None and reference_file is None:
        logger.error(f'Neither --reference-file nor --reference-name are specified. Please define either '
                     f'a new reference genome (--refrence-file [FILE]) or the name of an existing reference genome '
                     f'(--reference-name [NAME]). You can view previously-loaded reference genomes with '
                     f'"gdi --project-dir {project.get_root_dir()} list genomes"')
        sys.exit(1)
    elif reference_file is not None:
        reference_file = Path(reference_file)

    data_index_connection = project.create_connection()

    with TemporaryDirectory() as preprocess_dir:
        preprocess_dir = Path(preprocess_dir)
        data_package_factory = NucleotideInputFilesSampleDataPackageFactory(ncores=ncores, index_unknown=index_unknown,
                                                                            preprocess_dir=preprocess_dir,
                                                                            input_files_file=vcf_fofns)

        load_variants_common(data_index_connection=data_index_connection, ncores=ncores,
                             data_package_factory=data_package_factory,
                             reference_file=reference_file,
                             reference_name=reference_name,
                             input=Path(vcf_fofns), build_tree=build_tree, align_type=align_type,
                             extra_tree_params=extra_tree_params,
                             sample_batch_size=sample_batch_size)


@load.command(name='kmer')
@click.pass_context
@click.argument('kmer_fofns', type=click.Path(exists=True))
@click.option('--kmer-size', help='Kmer size for indexing. List multiple for multiple kmer sizes in an index',
              default=[31], multiple=True, type=click.IntRange(min=1, max=201))
def load_kmer(ctx, kmer_fofns, kmer_size):
    data_index_connection = get_project_exit_on_error(ctx).create_connection()
    filesystem_storage = data_index_connection.filesystem_storage
    kmer_service = data_index_connection.kmer_service

    if not isinstance(kmer_size, list):
        kmer_size = list(kmer_size)

    kmer_indexer = KmerIndexerSourmash(
        k=kmer_size,
        scaled=1000,
        abund=False,
        compress=True
    )

    index_manager = KmerIndexManager(filesystem_storage.kmer_dir, kmer_indexer=kmer_indexer)

    files_to_index = []

    files_df = pd.read_csv(kmer_fofns, sep='\t')
    for index, row in files_df.iterrows():
        sample_name = row['Sample']
        if kmer_service.has_kmer_index(sample_name):
            logger.warning(f'Sample [{sample_name}] already has kmer index, will not regenerate')
        else:
            genome_files = row['Files'].split(',')
            files_to_index.append((sample_name, genome_files))

    logger.info(f'Indexing {len(files_to_index)} genomes')

    indexed_genomes = index_manager.index_all_genomes(files_to_index)

    for sample_name in indexed_genomes:
        kmer_service.insert_kmer_index(sample_name=sample_name,
                                       kmer_index_path=indexed_genomes[sample_name])

    print(f'Generated indexes for {len(indexed_genomes)} samples')


@load.command(name='mlst-tseemann')
@click.pass_context
@click.argument('mlst_file', type=click.Path(exists=True), nargs=-1)
@click.option('--scheme-name', help='Override scheme name found in MLST file',
              default=FeatureService.AUTO_SCOPE, type=str)
def load_mlst_tseemann(ctx, mlst_file: List[Path], scheme_name: str):
    data_index_connection = get_project_exit_on_error(ctx).create_connection()
    mlst_service = cast(MLSTService, data_index_connection.mlst_service)
    for file in mlst_file:
        click.echo(f'Loading MLST results from {str(file)}')
        data_package = MLSTSampleDataPackage(MLSTTSeemannFeaturesReader(mlst_file=file))
        mlst_service.insert(data_package=data_package, feature_scope_name=scheme_name)


@load.command(name='mlst-sistr')
@click.pass_context
@click.argument('mlst_file', type=click.Path(exists=True), nargs=-1)
@click.option('--scheme-name', help='Override scheme name found in SISTR MLST file',
              default=FeatureService.AUTO_SCOPE, type=str)
def load_mlst_sistr(ctx, mlst_file: List[Path], scheme_name: str):
    data_index_connection = get_project_exit_on_error(ctx).create_connection()
    mlst_service = cast(MLSTService, data_index_connection.mlst_service)
    for file in mlst_file:
        click.echo(f'Loading cgMLST results from {str(file)}')
        data_package = MLSTSampleDataPackage(MLSTSistrReader(mlst_file=file))
        mlst_service.insert(data_package=data_package, feature_scope_name=scheme_name)


@load.command(name='mlst-chewbbaca')
@click.pass_context
@click.argument('mlst_file', type=click.Path(exists=True), nargs=-1)
@click.option('--scheme-name', help='Set scheme name',
              required=True, type=str)
def load_mlst_sistr(ctx, mlst_file: List[Path], scheme_name: str):
    data_index_connection = get_project_exit_on_error(ctx).create_connection()
    mlst_service = cast(MLSTService, data_index_connection.mlst_service)
    for file in mlst_file:
        click.echo(f'Loading MLST results from [{str(file)}] under scheme [{scheme_name}]')
        data_package = MLSTSampleDataPackage(MLSTChewbbacaReader(mlst_file=file, scheme=scheme_name))
        mlst_service.insert(data_package=data_package, feature_scope_name=FeatureService.AUTO_SCOPE)


@main.group(name='list')
@click.pass_context
def list_data(ctx):
    pass


@list_data.command(name='genomes')
@click.pass_context
def list_genomes(ctx):
    genomics_index = get_genomics_index(ctx)
    reference_genomes = genomics_index.reference_names()
    click.echo('\n'.join(reference_genomes))


@list_data.command(name='samples')
@click.pass_context
def list_samples(ctx):
    genomics_index = get_genomics_index(ctx)
    samples = genomics_index.sample_names()
    click.echo('\n'.join(samples))


def read_genomes_from_file(input_file: Path) -> List[Path]:
    with open(input_file, 'r') as fh:
        genome_paths = [Path(l.strip()) for l in fh.readlines()]
        return genome_paths


@main.command(name='input')
@click.option('--absolute/--no-absolute', help='Convert paths to absolute paths', required=False)
@click.option('--input-genomes-file',
              help='A file listing the genomes to process, one per line. This is an alternative'
                   ' to passing genomes as arguments on the command-line',
              type=click.Path(exists=True),
              required=False)
@click.argument('genomes', type=click.Path(exists=True), nargs=-1)
def input_command(absolute: bool, input_genomes_file: str, genomes: List[str]):
    if input_genomes_file is not None:
        if len(genomes) > 0:
            logger.warning(f'--input-genomes-file=[{input_genomes_file}] is specified so will ignore genomes '
                           f'passed on the command-line.')
        genome_paths = read_genomes_from_file(Path(input_genomes_file))
    elif len(genomes) > 0:
        genome_paths = [Path(f) for f in genomes]
    else:
        logger.error('Must pass a list of genome files or use --input-genomes-file')
        sys.exit(1)

    pipeline_executor = SnakemakePipelineExecutor()
    sample_files = pipeline_executor.create_input_sample_files(input_files=genome_paths)

    # Defaults to writing to stdout
    pipeline_executor.write_input_sample_files(input_sample_files=sample_files, abolute_paths=absolute)


@main.command()
@click.pass_context
@click.option('--reference-file', help='Reference genome', required=True, type=click.Path(exists=True))
@click.option('--load-data/--no-load-data', help='Whether or not to load the processed files into the index or'
                                                 ' just produce the VCFs from assemblies. --no-load-data implies '
                                                 '--no-clean.',
              default=True)
@click.option('--index-unknown/--no-index-unknown',
              help='Enable/disable indexing unknown/missing positions. Indexing missing positions can significantly '
                   'slow down the indexing process.',
              required=False, default=True)
@click.option('--clean/--no-clean', help='Clean up intermediate files when finished.', default=True)
@click.option('--build-tree/--no-build-tree', default=False, help='Builds tree of all samples after loading')
@click.option('--align-type', help=f'The type of alignment to generate', default='core',
              type=click.Choice(CoreAlignmentService.ALIGN_TYPES))
@click.option('--extra-tree-params', help='Extra parameters to tree-building software',
              default=None)
@click.option('--use-conda/--no-use-conda', help="Use (or don't use) conda for dependency management for pipeline.",
              default=False)
@click.option('--include-mlst/--no-include-mlst', help="Enable/disable including basic MLST in results.",
              default=False)
@click.option('--include-kmer/--no-include-kmer', help="Enable/disable including kmer analysis in results.",
              default=False)
@click.option('--ignore-snpeff/--no-ignore-snpeff', help="Enable/disable including snpeff annotations in "
                                                         "mutation results.", default=False)
@click.option('--reads-mincov', help='Minimum coverage when aligning reads to a reference genome',
              default=10, type=click.IntRange(min=1))
@click.option('--reads-minqual', help='Minimum quality score of VCF variants when aligning reads to a reference genome',
              default=100, type=click.IntRange(min=1))
@click.option('--kmer-size', help='Kmer size for indexing. List multiple for multiple kmer sizes in an index',
              default=[31], multiple=True, type=click.IntRange(min=1, max=201))
@click.option('--kmer-scaled', help='The scaled parameter to pass to sourmash. Defines how many kmers to keep in the '
                                    'sketch should be (i.e., a value of 1000 means to keep approx. 1/1000 kmers).',
              default=1000, type=click.IntRange(min=1))
@click.option('--batch-size', help='The maximum number of input files to process before dividing the Snakemake '
                                   'pipeline into batches. The number of jobs scheduled in each batch will be larger '
                                   'than this value.',
              default=2000, type=click.IntRange(min=1))
@click.option('--sample-batch-size', help='The maximum samples to load into an index at once. Increasing this value'
                                          ' may improve runtime at the expense of requiring more memory to construct '
                                          'the index. This only applies if loading data into the index is being done '
                                          'automatically (--load-data).', default=2000,
              type=click.IntRange(min=1))
@click.option('--input-genomes-file',
              help='A file listing the genomes to process, one per line. This is an alternative'
                   ' to passing genomes as arguments on the command-line',
              type=click.Path(exists=True),
              required=False)
@click.option('--input-structured-genomes-file',
              help='A structured file listing the samples and associated files. Used for finer control over sample names'
                   ' and the associated assemblies/reads. You can generate such a file with '
                   '"gdi input *.fasta *.fastq.gz > structured_input.tsv". This is an alternative'
                   ' to passing genomes as arguments on the command-line.',
              type=click.Path(exists=True),
              required=False)
@click.argument('genomes', type=click.Path(exists=True), nargs=-1)
def analysis(ctx, reference_file: str, load_data: bool, index_unknown: bool, clean: bool, build_tree: bool,
             align_type: str,
             extra_tree_params: str, use_conda: bool,
             include_mlst: bool, include_kmer: bool, ignore_snpeff: bool,
             reads_mincov: int, reads_minqual: int,
             kmer_size: List[int], kmer_scaled: int,
             batch_size: int, sample_batch_size: int,
             input_genomes_file: str, input_structured_genomes_file: str, genomes: List[str]):
    project = get_project_exit_on_error(ctx)
    data_index_connection = project.create_connection()
    kmer_service = data_index_connection.kmer_service

    if not load_data:
        logger.debug('--no-load-data is enabled so setting --no-clean')
        clean = False

    if not index_unknown:
        logger.info('--no-index-unknown is enabled so will not load unknown/missing mutation positions in index')

    sample_files = None
    genome_paths = []
    if input_structured_genomes_file is not None:
        logger.debug(f'Using --input-structured-genomes-file=[{input_structured_genomes_file}]')
        sample_files = SnakemakePipelineExecutor().read_input_sample_files(Path(input_structured_genomes_file))
    elif input_genomes_file is not None:
        logger.debug(f'Using --input-genomes-file=[{input_genomes_file}]')
        genome_paths = read_genomes_from_file(Path(input_genomes_file))
    else:
        logger.debug(f'Using {len(genomes)} files passed as command-line arguments')
        genome_paths = [Path(f) for f in genomes]

    timestamp = time.time()
    snakemake_directory = Path(getcwd(), f'snakemake-assemblies.{timestamp}')
    if not snakemake_directory.exists():
        mkdir(snakemake_directory)
    else:
        raise Exception(f'Snakemake working directory [{snakemake_directory}] already exists')

    if isinstance(kmer_size, int):
        kmer_sizes = [kmer_size]
    else:
        kmer_sizes = kmer_size

    pipeline_executor = SnakemakePipelineExecutor(working_directory=snakemake_directory,
                                                  use_conda=use_conda,
                                                  include_mlst=include_mlst,
                                                  include_kmer=include_kmer,
                                                  ignore_snpeff=ignore_snpeff,
                                                  kmer_sizes=kmer_sizes,
                                                  kmer_scaled=kmer_scaled,
                                                  snakemake_input_batch_size=batch_size,
                                                  reads_mincov=reads_mincov,
                                                  reads_minqual=reads_minqual)

    if sample_files is None:
        logger.info(f'Automatically structuring {len(genome_paths)} input files into assemblies/reads')
        sample_files = pipeline_executor.create_input_sample_files(genome_paths)

    logger.info(f'Processing {len(sample_files)} genomes to identify mutations')
    results = pipeline_executor.execute(sample_files=sample_files,
                                        reference_file=Path(reference_file),
                                        ncores=ctx.obj['ncores'])

    processed_files_fofn = results.get_file('gdi-fofn')

    if load_data:
        try:
            logger.info(f'Indexing processed VCF files defined in [{processed_files_fofn}]')
            ctx.invoke(load_vcf, index_unknown=index_unknown, vcf_fofns=str(processed_files_fofn),
                       reference_file=reference_file, build_tree=build_tree, align_type=align_type,
                       extra_tree_params=extra_tree_params, sample_batch_size=sample_batch_size)
        except Exception as e:
            logger.exception(e)
            logger.error(f"Error while indexing. Please verify files in [{snakemake_directory}] are correct.")
            clean = False

        if include_kmer:
            logger.info(f'Inserting kmer sketches for {len(sample_files)} samples into the database')
            files_df = pd.read_csv(processed_files_fofn, sep='\t', index_col=False)
            for idx, row in files_df.iterrows():
                sample_name = row['Sample']
                kmer_sketch = row['Sketch File']
                kmer_service.insert_kmer_index(sample_name=sample_name,
                                               kmer_index_path=Path(kmer_sketch))

        if clean:
            logger.info(f'--clean is enabled so deleting [{snakemake_directory}]')
            shutil.rmtree(snakemake_directory)
    else:
        logger.debug(f'Not indexing processed files defined in [{processed_files_fofn}]')
        reference_name = SequenceFile(Path(reference_file)).get_genome_name()
        click.echo(f'Processed files found in: {processed_files_fofn}')
        click.echo(f'Load with: "gdi --project-dir {project.get_root_dir()} load vcf '
                   f'--reference-name {reference_name} {processed_files_fofn}"')

        if include_mlst:
            mlst_file = results.get_file('mlst')
            click.echo(f'MLST results found in: {mlst_file}')
            click.echo(f'Load with: "gdi --project-dir {project.get_root_dir()} load mlst-tseemann '
                       f'{mlst_file}"')


@main.group()
@click.pass_context
def export(ctx):
    pass


@export.command(name='tree')
@click.pass_context
@click.argument('name', nargs=-1)
@click.option('--ascii/--no-ascii', help='Export as ASCII figure')
def export_tree(ctx, name: List[str], ascii: bool):
    genomics_index = get_genomics_index(ctx)
    if len(name) == 0:
        logger.warning('No reference genome names passed, will not export tree')

    for ref_name in name:
        try:
            reference_tree = genomics_index.reference_tree(ref_name)
            if ascii:
                click.echo(str(reference_tree))
            else:
                click.echo(reference_tree.write())
        except EntityExistsError as e:
            logger.error(str(e))


@main.group()
@click.pass_context
def build(ctx):
    pass


@build.command()
@click.pass_context
@click.option('--output-file', help='Output file', required=True, type=click.Path())
@click.option('--reference-name', help='Reference genome name', required=True, type=str)
@click.option('--align-type', help=f'The type of alignment to generate', default='core',
              type=click.Choice(CoreAlignmentService.ALIGN_TYPES))
@click.option('--sample', help='Sample to include in alignment (can list more than one).',
              multiple=True, type=str)
def alignment(ctx, output_file: Path, reference_name: str, align_type: str, sample: List[str]):
    genomics_index = get_genomics_index(ctx)
    references = genomics_index.reference_names()
    alignment_service = genomics_index.connection.alignment_service

    if reference_name not in references:
        logger.error(f'Reference genome [{reference_name}] does not exist')
        sys.exit(1)

    query = genomics_index.samples_query().isin(list(sample), kind='samples')

    if len(query) != len(sample):
        found_samples = query.toset(names=True)
        logger.error(f'Samples {set(sample) - found_samples} do not exist')
        sys.exit(1)

    alignment_data = alignment_service.construct_alignment(reference_name=reference_name,
                                                           samples=query.tolist(names=True),
                                                           align_type=align_type,
                                                           include_reference=True)

    with open(output_file, 'w') as f:
        AlignIO.write(alignment_data, f, 'fasta')
        click.echo(f'Wrote alignment to [{output_file}]')


supported_tree_build_types = ['iqtree']


@build.command()
@click.pass_context
@click.option('--output-file', help='Output file', required=True, type=click.Path())
@click.option('--reference-name', help='Reference genome name', type=str, required=True)
@click.option('--align-type', help=f'The type of alignment to use for generating the tree', default='core',
              type=click.Choice(CoreAlignmentService.ALIGN_TYPES))
@click.option('--tree-build-type', help=f'The type of tree building software', default='iqtree',
              type=click.Choice(supported_tree_build_types))
@click.option('--sample', help='Sample to include in tree (can list more than one).',
              multiple=True, type=str)
@click.option('--extra-params', help='Extra parameters to tree-building software',
              default=None)
def tree(ctx, output_file: Path, reference_name: str, align_type: str,
         tree_build_type: str, sample: List[str], extra_params: str):
    genomics_index = get_genomics_index(ctx)
    references = genomics_index.reference_names()
    ncores = ctx.obj['ncores']

    if reference_name not in references:
        logger.error(f'Reference genome [{reference_name}] does not exist')
        sys.exit(1)

    query = genomics_index.samples_query().isin(list(sample), kind='samples')

    if len(query) != len(sample):
        found_samples = query.toset(names=True)
        logger.error(f'Samples {set(sample) - found_samples} do not exist')
        sys.exit(1)

    # Eventually I want to add full support for fasttree/other tree builders, so I'm
    # leaving this if/else here
    if align_type == 'full' and tree_build_type == 'fasttree':
        logger.error(f'align_type=[{align_type}] is not supported for tree_build_type=[{tree_build_type}]')
        sys.exit(1)

    tree_query = query.build_tree(kind='mutation',
                                  method=tree_build_type,
                                  align_type=align_type,
                                  scope=reference_name,
                                  include_reference=True,
                                  ncores=ncores,
                                  extra_params=extra_params)

    tree_query.tree.write(outfile=output_file)
    click.echo(f'Wrote tree to [{output_file}]')


@main.group()
@click.pass_context
def rebuild(ctx):
    pass


@rebuild.command(name='tree')
@click.pass_context
@click.argument('reference', type=str, nargs=-1)
@click.option('--align-type', help=f'The type of alignment to use for generating the tree', default='core',
              type=click.Choice(CoreAlignmentService.ALIGN_TYPES))
@click.option('--extra-params', help='Extra parameters to tree-building software',
              default=None)
def rebuild_tree(ctx, reference: List[str], align_type: str, extra_params: str):
    data_index_connection = get_project_exit_on_error(ctx).create_connection()
    tree_service = data_index_connection.tree_service
    reference_service = data_index_connection.reference_service
    ncores = ctx.obj['ncores']

    if len(reference) == 0:
        logger.error('Must define name of reference genome to use. '
                     'To see available genomes try "variants list genomes"')
        sys.exit(1)

    for reference_name in reference:
        if not reference_service.exists_reference_genome(reference_name):
            logger.error(f'Reference genome [{reference_name}] does not exist')
            sys.exit(1)

    for reference_name in reference:
        logger.info(f'Started rebuilding tree for reference genome [{reference_name}]')
        tree_service.rebuild_tree(reference_name=reference_name,
                                  align_type=align_type,
                                  num_cores=ncores,
                                  extra_params=extra_params)
        logger.info(f'Finished rebuilding tree')


@main.group()
@click.pass_context
def query(ctx):
    pass


def query_feature(genomics_index: GenomicsDataIndex, features: List[QueryFeature], summarize: bool) -> None:
    query = genomics_index.samples_query()
    for feature in features:
        query = query.hasa(feature)

    if summarize:
        results_df = query.summary()
    else:
        results_df = query.toframe(include_unknown=True)

    results_df.to_csv(sys.stdout, sep='\t', index=False, float_format='%0.4g', na_rep='-')


@query.command(name='mutation')
@click.pass_context
@click.argument('name', nargs=-1)
@click.option('--summarize/--no-summarize', help='Print summary information on query')
def query_mutation(ctx, name: List[str], summarize: bool):
    genomics_index = get_genomics_index(ctx)
    features = [QueryFeatureMutationSPDI(n) for n in name]
    query_feature(genomics_index=genomics_index, features=features, summarize=summarize)


@query.command(name='mlst')
@click.pass_context
@click.argument('name', nargs=-1)
@click.option('--summarize/--no-summarize', help='Print summary information on query')
def query_mlst(ctx, name: List[str], summarize: bool):
    genomics_index = get_genomics_index(ctx)
    features = [QueryFeatureMLST(n) for n in name]
    query_feature(genomics_index=genomics_index, features=features, summarize=summarize)


@main.group(name='db')
@click.pass_context
def db(ctx):
    pass


UNITS = ['B', 'KB', 'MB', 'GB']


@db.command(name='size')
@click.pass_context
@click.option('--unit', default='B', help='The unit to display data sizes as.', type=click.Choice(UNITS))
def db_size(ctx, unit):
    genomics_index = get_genomics_index(ctx)
    size_df = genomics_index.db_size(unit)
    size_df.to_csv(sys.stdout, sep='\t', index=False, float_format='%0.2f', na_rep='-')
