import logging
import multiprocessing
import shutil
import sys
import time
from functools import partial
from os import path, listdir, getcwd, mkdir
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
from genomics_data_index.cli import yaml_config_provider
from genomics_data_index.configuration.Project import Project, ProjectConfigurationError
from genomics_data_index.configuration.connector.DataIndexConnection import DataIndexConnection
from genomics_data_index.pipelines.SnakemakePipelineExecutor import SnakemakePipelineExecutor
from genomics_data_index.storage.index.KmerIndexer import KmerIndexerSourmash, KmerIndexManager
from genomics_data_index.storage.io.mlst.MLSTChewbbacaReader import MLSTChewbbacaReader
from genomics_data_index.storage.io.mlst.MLSTSampleDataPackage import MLSTSampleDataPackage
from genomics_data_index.storage.io.mlst.MLSTSistrReader import MLSTSistrReader
from genomics_data_index.storage.io.mlst.MLSTTSeemannFeaturesReader import MLSTTSeemannFeaturesReader
from genomics_data_index.storage.io.mutation.NucleotideSampleDataPackage import NucleotideSampleDataPackage
from genomics_data_index.storage.io.mutation.SequenceFile import SequenceFile
from genomics_data_index.storage.io.processor.MultipleProcessSampleFilesProcessor import \
    MultipleProcessSampleFilesProcessor
from genomics_data_index.storage.io.processor.NullSampleFilesProcessor import NullSampleFilesProcessor
from genomics_data_index.storage.model.QueryFeatureMLST import QueryFeatureMLST
from genomics_data_index.storage.model.QueryFeatureMutationSPDI import QueryFeatureMutationSPDI
from genomics_data_index.storage.service import EntityExistsError
from genomics_data_index.storage.service.CoreAlignmentService import CoreAlignmentService
from genomics_data_index.storage.service.MLSTService import MLSTService
from genomics_data_index.storage.service.SampleService import SampleService
from genomics_data_index.storage.service.TreeService import TreeService
from genomics_data_index.storage.service.VariationService import VariationService

logger = logging.getLogger('genomics_data_index')
max_cores = multiprocessing.cpu_count()

LOG_LEVELS = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']

click.option = partial(click.option, show_default=True)


def setup_logging(log_level: str) -> None:
    if log_level == 'DEBUG':
        log_format = '%(asctime)s %(levelname)s %(name)s.%(funcName)s,%(lineno)s: %(message)s'
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
                         data_package: NucleotideSampleDataPackage, reference_file, input, build_tree, align_type,
                         extra_tree_params: str):
    reference_service = data_index_connection.reference_service
    variation_service = cast(VariationService, data_index_connection.variation_service)
    sample_service = cast(SampleService, data_index_connection.sample_service)
    tree_service = data_index_connection.tree_service

    try:
        reference_service.add_reference_genome(reference_file)
    except EntityExistsError as e:
        logger.warning(f'Reference genome [{reference_file}] already exists, will not load')

    samples_exist = sample_service.which_exists(list(data_package.sample_names()))
    if len(samples_exist) > 0:
        logger.error(f'Samples {samples_exist} already exist, will not load any variants')
    else:
        reference_name = SequenceFile(reference_file).get_genome_name()

        variation_service.insert(feature_scope_name=reference_name,
                                 data_package=data_package)
        click.echo(f'Loaded variants from [{input}] into database')

        if build_tree:
            tree_service.rebuild_tree(reference_name=reference_name,
                                      align_type=align_type,
                                      num_cores=ncores,
                                      extra_params=extra_tree_params)
            click.echo('Finished building tree of all samples')


@load.command(name='snippy')
@click.pass_context
@click.argument('snippy_dir', type=click.Path(exists=True))
@click.option('--reference-file', help='Reference genome', required=True, type=click.Path(exists=True))
@click.option('--build-tree/--no-build-tree', default=False, help='Builds tree of all samples after loading')
@click.option('--align-type', help=f'The type of alignment to generate', default='core',
              type=click.Choice(CoreAlignmentService.ALIGN_TYPES))
@click.option('--extra-tree-params', help='Extra parameters to tree-building software',
              default=None)
def load_snippy(ctx, snippy_dir: Path, reference_file: Path, build_tree: bool, align_type: str, extra_tree_params: str):
    ncores = ctx.obj['ncores']
    project = get_project_exit_on_error(ctx)
    data_index_connection = project.create_connection()

    snippy_dir = Path(snippy_dir)
    reference_file = Path(reference_file)
    click.echo(f'Loading {snippy_dir}')
    sample_dirs = [snippy_dir / d for d in listdir(snippy_dir) if path.isdir(snippy_dir / d)]

    with TemporaryDirectory() as preprocess_dir:
        if ncores > 1:
            file_processor = MultipleProcessSampleFilesProcessor(preprocess_dir=Path(preprocess_dir),
                                                                 processing_cores=ncores)
        else:
            file_processor = NullSampleFilesProcessor.instance()

        data_package = NucleotideSampleDataPackage.create_from_snippy(sample_dirs,
                                                                      sample_files_processor=file_processor)

        load_variants_common(data_index_connection=data_index_connection, ncores=ncores, data_package=data_package,
                             reference_file=reference_file,
                             input=snippy_dir, build_tree=build_tree, align_type=align_type,
                             extra_tree_params=extra_tree_params)


@load.command(name='vcf')
@click.pass_context
@click.argument('vcf_fofns', type=click.Path(exists=True))
@click.option('--reference-file', help='Reference genome', required=True, type=click.Path(exists=True))
@click.option('--build-tree/--no-build-tree', default=False, help='Builds tree of all samples after loading')
@click.option('--align-type', help=f'The type of alignment to generate', default='core',
              type=click.Choice(CoreAlignmentService.ALIGN_TYPES))
@click.option('--extra-tree-params', help='Extra parameters to tree-building software',
              default=None)
def load_vcf(ctx, vcf_fofns: str, reference_file: str, build_tree: bool, align_type: str, extra_tree_params: str):
    ncores = ctx.obj['ncores']
    vcf_fofns = Path(vcf_fofns)
    reference_file = Path(reference_file)
    project = get_project_exit_on_error(ctx)
    data_index_connection = project.create_connection()

    click.echo(f'Loading files listed in {vcf_fofns}')
    sample_vcf_map = {}
    mask_files_map = {}
    files_df = pd.read_csv(vcf_fofns, sep='\t')
    for index, row in files_df.iterrows():
        if row['Sample'] in sample_vcf_map:
            raise Exception(f'Error, duplicate samples {row["Sample"]} in file {vcf_fofns}')

        sample_vcf_map[row['Sample']] = row['VCF']
        if not pd.isna(row['Mask File']):
            mask_files_map[row['Sample']] = row['Mask File']

    with TemporaryDirectory() as preprocess_dir:
        if ncores > 1:
            file_processor = MultipleProcessSampleFilesProcessor(preprocess_dir=Path(preprocess_dir),
                                                                 processing_cores=ncores)
        else:
            file_processor = NullSampleFilesProcessor.instance()

        data_package = NucleotideSampleDataPackage.create_from_sequence_masks(sample_vcf_map=sample_vcf_map,
                                                                              masked_genomic_files_map=mask_files_map,
                                                                              sample_files_processor=file_processor)

        load_variants_common(data_index_connection=data_index_connection, ncores=ncores, data_package=data_package,
                             reference_file=reference_file,
                             input=Path(vcf_fofns), build_tree=build_tree, align_type=align_type,
                             extra_tree_params=extra_tree_params)


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
    data_index_connection = get_project_exit_on_error(ctx).create_connection()
    items = [genome.name for genome in data_index_connection.reference_service.get_reference_genomes()]
    click.echo('\n'.join(items))


@list_data.command(name='samples')
@click.pass_context
def list_samples(ctx):
    data_index_connection = get_project_exit_on_error(ctx).create_connection()
    items = [sample.name for sample in data_index_connection.sample_service.get_samples()]
    click.echo('\n'.join(items))


def read_genomes_from_file(input_file: Path) -> List[Path]:
    with open(input_file, 'r') as fh:
        genome_paths = [Path(l.strip()) for l in fh.readlines()]
        return genome_paths

@main.command()
@click.option('--absolute/--no-absolute', help='Convert paths to absolute paths', required=False)
@click.option('--input-genomes-file',
              help='A file listing the genomes to process, one per line. This is an alternative'
                   ' to passing genomes as arguments on the command-line',
              type=click.Path(exists=True),
              required=False)
@click.argument('genomes', type=click.Path(exists=True), nargs=-1)
def input(absolute: bool, input_genomes_file: str, genomes: List[str]):
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
@click.option('--reference-file', help='Reference genome', required=True, type=click.Path(exists=True))
@click.option('--index/--no-index', help='Whether or not to load the processed files into the index or'
                                         ' just produce the VCFs from assemblies. --no-index implies --no-clean.',
              default=True)
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
@click.option('--kmer-size', help='Kmer size for indexing. List multiple for multiple kmer sizes in an index',
              default=[31], multiple=True, type=click.IntRange(min=1, max=201))
@click.option('--kmer-scaled', help='The scaled parameter to pass to sourmash. Defines how many kmers to keep in the '
                                    'sketch should be (i.e., a value of 1000 means to keep approx. 1/1000 kmers).',
              default=1000, type=click.IntRange(min=1))
@click.option('--batch-size', help='The maximum number of input files to process before dividing the Snakemake '
                                   'pipeline into batches. The number of jobs scheduled in each batch will be larger '
                                   'than this value.',
              default=2000, type=click.IntRange(min=1))
@click.option('--input-genomes-file',
              help='A file listing the genomes to process, one per line. This is an alternative'
                   ' to passing genomes as arguments on the command-line',
              type=click.Path(exists=True),
              required=False)
@click.argument('genomes', type=click.Path(exists=True), nargs=-1)
def analysis(ctx, reference_file: str, index: bool, clean: bool, build_tree: bool, align_type: str,
             extra_tree_params: str, use_conda: bool,
             include_mlst: bool, include_kmer: bool, ignore_snpeff: bool, kmer_size: List[int], kmer_scaled: int,
             batch_size: int,
             input_genomes_file: str, genomes: List[str]):
    data_index_connection = get_project_exit_on_error(ctx).create_connection()
    kmer_service = data_index_connection.kmer_service

    if not index:
        logger.debug('--no-index is enabled so setting --no-clean')
        clean = False

    if input_genomes_file is not None:
        genome_paths = read_genomes_from_file(Path(input_genomes_file))
    else:
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
                                                  snakemake_input_batch_size=batch_size)

    logger.info(f'Automatically structuring {len(genome_paths)} input files into assemblies/reads')
    sample_files = pipeline_executor.create_input_sample_files(genome_paths)

    logger.info(f'Processing {len(genome_paths)} genomes to identify mutations')
    results = pipeline_executor.execute(sample_files=sample_files,
                                        reference_file=Path(reference_file),
                                        ncores=ctx.obj['ncores'])

    processed_files_fofn = results.get_file('gdi-fofn')

    if index:
        try:
            logger.info(f'Indexing processed VCF files defined in [{processed_files_fofn}]')
            ctx.invoke(load_vcf, vcf_fofns=str(processed_files_fofn), reference_file=reference_file,
                       build_tree=build_tree, align_type=align_type, extra_tree_params=extra_tree_params)
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
        click.echo(f'Processed files found in: {processed_files_fofn}')

        if include_mlst:
            mlst_file = results.get_file('mlst')
            click.echo(f'MLST results found in: {mlst_file}')


@main.group()
@click.pass_context
def export(ctx):
    pass


@export.command(name='tree')
@click.pass_context
@click.argument('name', nargs=-1)
@click.option('--ascii/--no-ascii', help='Export as ASCII figure')
def export_tree(ctx, name: List[str], ascii: bool):
    data_index_connection = get_project_exit_on_error(ctx).create_connection()
    if len(name) == 0:
        logger.warning('No reference genome names passed, will not export tree')

    for ref_name in name:
        reference = data_index_connection.reference_service.find_reference_genome(ref_name)
        if ascii:
            click.echo(str(reference.tree))
        else:
            click.echo(reference.tree.write())


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
    data_index_connection = get_project_exit_on_error(ctx).create_connection()
    alignment_service = data_index_connection.alignment_service
    reference_service = data_index_connection.reference_service
    sample_service = data_index_connection.sample_service

    if not reference_service.exists_reference_genome(reference_name):
        logger.error(f'Reference genome [{reference_name}] does not exist')
        sys.exit(1)

    found_samples = set(sample_service.which_exists(sample))

    if len(sample) > 0 and found_samples != set(sample):
        logger.error(f'Samples {set(sample) - found_samples} do not exist')
        sys.exit(1)

    alignment_data = alignment_service.construct_alignment(reference_name=reference_name,
                                                           samples=sample,
                                                           align_type=align_type,
                                                           include_reference=True)

    with open(output_file, 'w') as f:
        AlignIO.write(alignment_data, f, 'fasta')
        click.echo(f'Wrote alignment to [{output_file}]')


@build.command()
@click.pass_context
@click.option('--output-file', help='Output file', required=True, type=click.Path())
@click.option('--reference-name', help='Reference genome name', type=str, required=True)
@click.option('--align-type', help=f'The type of alignment to use for generating the tree', default='core',
              type=click.Choice(CoreAlignmentService.ALIGN_TYPES))
@click.option('--tree-build-type', help=f'The type of tree building software', default='iqtree',
              type=click.Choice(TreeService.TREE_BUILD_TYPES))
@click.option('--sample', help='Sample to include in tree (can list more than one).',
              multiple=True, type=str)
@click.option('--extra-params', help='Extra parameters to tree-building software',
              default=None)
def tree(ctx, output_file: Path, reference_name: str, align_type: str,
         tree_build_type: str, sample: List[str], extra_params: str):
    data_index_connection = get_project_exit_on_error(ctx).create_connection()
    alignment_service = data_index_connection.alignment_service
    tree_service = data_index_connection.tree_service
    reference_service = ctx.obj['data_index_connection'].reference_service
    sample_service = ctx.obj['data_index_connection'].sample_service
    ncores = ctx.obj['ncores']

    if not reference_service.exists_reference_genome(reference_name):
        logger.error(f'Reference genome [{reference_name}] does not exist')
        sys.exit(1)

    found_samples = set(sample_service.which_exists(sample))

    if len(sample) > 0 and found_samples != set(sample):
        logger.error(f'Samples {set(sample) - found_samples} do not exist')
        sys.exit(1)

    if align_type == 'full' and tree_build_type == 'fasttree':
        logger.error(f'align_type=[{align_type}] is not supported for tree_build_type=[{tree_build_type}]')
        sys.exit(1)

    alignment_data = alignment_service.construct_alignment(reference_name=reference_name,
                                                           samples=sample,
                                                           align_type=align_type,
                                                           include_reference=True)

    log_file = f'{output_file}.log'

    tree_data, out = tree_service.build_tree(alignment_data, tree_build_type=tree_build_type,
                                             num_cores=ncores, align_type=align_type, extra_params=extra_params)
    tree_data.write(outfile=output_file)
    click.echo(f'Wrote tree to [{output_file}]')
    with open(log_file, 'w') as log:
        log.write(out)
        click.echo(f'Wrote log file to [{log_file}]')


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


@query.command(name='sample-mutation')
@click.pass_context
@click.argument('name', nargs=-1)
def query_sample_mutation(ctx, name: List[str]):
    data_index_connection = get_project_exit_on_error(ctx).create_connection()
    mutation_query_service = data_index_connection.mutation_query_service
    match_df = mutation_query_service.find_matches(samples=name)
    match_df.to_csv(sys.stdout, sep='\t', index=False, float_format='%0.4g', na_rep='-')


@query.command(name='sample-kmer')
@click.pass_context
@click.argument('name', nargs=-1)
def query_sample_kmer(ctx, name: List[str]):
    data_index_connection = get_project_exit_on_error(ctx).create_connection()
    kmer_query_service = data_index_connection.kmer_query_service
    match_df = kmer_query_service.find_matches(samples=name)
    match_df.to_csv(sys.stdout, sep='\t', index=False, float_format='%0.4g', na_rep='-')


@query.command(name='mutation')
@click.pass_context
@click.argument('name', nargs=-1)
@click.option('--include-unknown/--no-include-unknown',
              help='Including results where it is unknown if the search term is present or not.',
              required=False)
@click.option('--summarize/--no-summarize', help='Print summary information on query')
def query_mutation(ctx, name: List[str], include_unknown: bool, summarize: bool):
    data_index_connection = get_project_exit_on_error(ctx).create_connection()
    mutation_query_service = data_index_connection.mutation_query_service

    features = [QueryFeatureMutationSPDI(n) for n in name]
    if not summarize:
        match_df = mutation_query_service.find_by_features(features, include_unknown=include_unknown)
    else:
        match_df = mutation_query_service.count_by_features(features, include_unknown=include_unknown)

    match_df.to_csv(sys.stdout, sep='\t', index=False, float_format='%0.4g', na_rep='-')


@query.command(name='mlst')
@click.pass_context
@click.argument('name', nargs=-1)
@click.option('--include-unknown/--no-include-unknown',
              help='Including results where it is unknown if the search term is present or not.',
              required=False)
@click.option('--summarize/--no-summarize', help='Print summary information on query')
def query_mlst(ctx, name: List[str], include_unknown: bool, summarize: bool):
    data_index_connection = get_project_exit_on_error(ctx).create_connection()
    mlst_query_service = data_index_connection.mlst_query_service

    features = [QueryFeatureMLST(n) for n in name]
    if not summarize:
        match_df = mlst_query_service.find_by_features(features, include_unknown=include_unknown)
    else:
        match_df = mlst_query_service.count_by_features(features, include_unknown=include_unknown)

    match_df.to_csv(sys.stdout, sep='\t', index=False, float_format='%0.4g', na_rep='-')


@main.group(name='db')
@click.pass_context
def db(ctx):
    pass


UNITS = ['B', 'KB', 'MB', 'GB']


@db.command(name='size')
@click.pass_context
@click.option('--unit', default='B', help='The unit to display data sizes as.', type=click.Choice(UNITS))
def db_size(ctx, unit):
    data_index_connection = get_project_exit_on_error(ctx).create_connection()
    size_df = data_index_connection.db_size(unit)
    size_df.to_csv(sys.stdout, sep='\t', index=False, float_format='%0.2f', na_rep='-')
