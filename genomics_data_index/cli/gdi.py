import logging
import multiprocessing
import sys
from functools import partial
from os import path, listdir, getcwd
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import List, cast

import click
import click_config_file
import coloredlogs
import pandas as pd
from Bio import AlignIO

import genomics_data_index.storage.service.FeatureService as FeatureService
from genomics_data_index.cli import yaml_config_provider
from genomics_data_index.configuration.Project import Project
from genomics_data_index.storage.index.KmerIndexer import KmerIndexerSourmash, KmerIndexManager
from genomics_data_index.storage.io.mlst.MLSTChewbbacaReader import MLSTChewbbacaReader
from genomics_data_index.storage.io.mlst.MLSTSampleDataPackage import MLSTSampleDataPackage
from genomics_data_index.storage.io.mlst.MLSTSistrReader import MLSTSistrReader
from genomics_data_index.storage.io.mlst.MLSTTSeemannFeaturesReader import MLSTTSeemannFeaturesReader
from genomics_data_index.storage.io.mutation.NucleotideSampleDataPackage import NucleotideSampleDataPackage
from genomics_data_index.storage.io.processor.MultipleProcessSampleFilesProcessor import \
    MultipleProcessSampleFilesProcessor
from genomics_data_index.storage.io.processor.NullSampleFilesProcessor import NullSampleFilesProcessor
from genomics_data_index.storage.model.QueryFeatureMLST import QueryFeatureMLST
from genomics_data_index.storage.model.QueryFeatureMutation import QueryFeatureMutation
from genomics_data_index.storage.service import EntityExistsError
from genomics_data_index.storage.service.CoreAlignmentService import CoreAlignmentService
from genomics_data_index.storage.service.MLSTService import MLSTService
from genomics_data_index.storage.service.SampleService import SampleService
from genomics_data_index.storage.service.TreeService import TreeService
from genomics_data_index.storage.service.VariationService import VariationService
from genomics_data_index.storage.util import get_genome_name

logger = logging.getLogger('genomics_data_index')
max_cores = multiprocessing.cpu_count()

LOG_LEVELS = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']

click.option = partial(click.option, show_default=True)


def setup_logging(log_level: str) -> None:
    if log_level == 'INFO' or log_level == 'DEBUG':
        log_format = '%(asctime)s %(levelname)s %(name)s.%(funcName)s,%(lineno)s: %(message)s'
    else:
        log_format = '%(asctime)s %(levelname)s: %(message)s'

    coloredlogs.install(level=log_level, fmt=log_format, logger=logger)


@click.group()
@click.pass_context
@click.option('--project-dir', help='A project directory containing the data and connection information.')
@click.option('--ncores', help='Number of cores for any parallel processing', default=max_cores,
              type=click.IntRange(min=1, max=max_cores))
@click.option('--log-level', default='WARNING', help='Sets the log level', type=click.Choice(LOG_LEVELS))
@click_config_file.configuration_option(provider=yaml_config_provider,
                                        config_file_name='config.yaml',
                                        implicit=True)
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
    project = Project(ctx.obj['project_dir'])
    ctx.obj['data_index_connection'] = project.create_connection()


def load_variants_common(ctx, data_package: NucleotideSampleDataPackage, reference_file, input, build_tree, align_type,
                         extra_tree_params: str):
    reference_service = ctx.obj['data_index_connection'].reference_service
    variation_service = cast(VariationService, ctx.obj['data_index_connection'].variation_service)
    sample_service = cast(SampleService, ctx.obj['data_index_connection'].sample_service)
    tree_service = ctx.obj['data_index_connection'].tree_service
    ncores = ctx.obj['ncores']

    try:
        reference_service.add_reference_genome(reference_file)
    except EntityExistsError as e:
        logger.warning(f'Reference genome [{reference_file}] already exists, will not load')

    samples_exist = sample_service.which_exists(list(data_package.sample_names()))
    if len(samples_exist) > 0:
        logger.error(f'Samples {samples_exist} already exist, will not load any variants')
    else:
        reference_name = get_genome_name(reference_file)

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

        load_variants_common(ctx=ctx, data_package=data_package, reference_file=reference_file,
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
def load_vcf(ctx, vcf_fofns: Path, reference_file: Path, build_tree: bool, align_type: str, extra_tree_params: str):
    reference_file = Path(reference_file)
    ncores = ctx.obj['ncores']

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

        load_variants_common(ctx=ctx, data_package=data_package, reference_file=reference_file,
                             input=Path(vcf_fofns), build_tree=build_tree, align_type=align_type,
                             extra_tree_params=extra_tree_params)


@load.command(name='kmer')
@click.pass_context
@click.argument('kmer_fofns', type=click.Path(exists=True))
@click.option('--kmer-size', help='Kmer size for indexing. List multiple for multiple kmer sizes in an index',
              default=[31], multiple=True, type=click.IntRange(min=1, max=201))
def load_kmer(ctx, kmer_fofns, kmer_size):
    filesystem_storage = ctx.obj['data_index_connection'].filesystem_storage
    kmer_service = ctx.obj['data_index_connection'].kmer_service

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
    mlst_service = cast(MLSTService, ctx.obj['data_index_connection'].mlst_service)
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
    mlst_service = cast(MLSTService, ctx.obj['data_index_connection'].mlst_service)
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
    mlst_service = cast(MLSTService, ctx.obj['data_index_connection'].mlst_service)
    for file in mlst_file:
        click.echo(f'Loading MLST results from [{str(file)}] under scheme [{scheme_name}]')
        data_package = MLSTSampleDataPackage(MLSTChewbbacaReader(mlst_file=file, scheme=scheme_name))
        mlst_service.insert(data_package=data_package, feature_scope_name=FeatureService.AUTO_SCOPE)


@main.group(name='list')
@click.pass_context
def list_data(ctx):
    project = Project(ctx.obj['project_dir'])
    ctx.obj['data_index_connection'] = project.create_connection()


@list_data.command(name='genomes')
@click.pass_context
def list_genomes(ctx):
    items = [genome.name for genome in ctx.obj['data_index_connection'].reference_service.get_reference_genomes()]
    click.echo('\n'.join(items))


@list_data.command(name='samples')
@click.pass_context
def list_samples(ctx):
    items = [sample.name for sample in ctx.obj['data_index_connection'].sample_service.get_samples()]
    click.echo('\n'.join(items))


@main.group()
@click.pass_context
def export(ctx):
    project = Project(ctx.obj['project_dir'])
    ctx.obj['data_index_connection'] = project.create_connection()


@export.command(name='tree')
@click.pass_context
@click.argument('name', nargs=-1)
@click.option('--ascii/--no-ascii', help='Export as ASCII figure')
def export_tree(ctx, name: List[str], ascii: bool):
    if len(name) == 0:
        logger.warning('No reference genome names passed, will not export tree')

    for ref_name in name:
        reference = ctx.obj['data_index_connection'].reference_service.find_reference_genome(ref_name)
        if ascii:
            click.echo(str(reference.tree))
        else:
            click.echo(reference.tree.write())


@main.group()
@click.pass_context
def build(ctx):
    project = Project(ctx.obj['project_dir'])
    ctx.obj['data_index_connection'] = project.create_connection()


@build.command()
@click.pass_context
@click.option('--output-file', help='Output file', required=True, type=click.Path())
@click.option('--reference-name', help='Reference genome name', required=True, type=str)
@click.option('--align-type', help=f'The type of alignment to generate', default='core',
              type=click.Choice(CoreAlignmentService.ALIGN_TYPES))
@click.option('--sample', help='Sample to include in alignment (can list more than one).',
              multiple=True, type=str)
def alignment(ctx, output_file: Path, reference_name: str, align_type: str, sample: List[str]):
    alignment_service = ctx.obj['data_index_connection'].alignment_service
    reference_service = ctx.obj['data_index_connection'].reference_service
    sample_service = ctx.obj['data_index_connection'].sample_service

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
    alignment_service = ctx.obj['data_index_connection'].alignment_service
    tree_service = ctx.obj['data_index_connection'].tree_service
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
    project = Project(ctx.obj['project_dir'])
    ctx.obj['data_index_connection'] = project.create_connection()


@rebuild.command(name='tree')
@click.pass_context
@click.argument('reference', type=str, nargs=-1)
@click.option('--align-type', help=f'The type of alignment to use for generating the tree', default='core',
              type=click.Choice(CoreAlignmentService.ALIGN_TYPES))
@click.option('--extra-params', help='Extra parameters to tree-building software',
              default=None)
def rebuild_tree(ctx, reference: List[str], align_type: str, extra_params: str):
    tree_service = ctx.obj['data_index_connection'].tree_service
    reference_service = ctx.obj['data_index_connection'].reference_service
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
    project = Project(ctx.obj['project_dir'])
    ctx.obj['data_index_connection'] = project.create_connection()


@query.command(name='sample-mutation')
@click.pass_context
@click.argument('name', nargs=-1)
def query_sample_mutation(ctx, name: List[str]):
    mutation_query_service = ctx.obj['data_index_connection'].mutation_query_service
    match_df = mutation_query_service.find_matches(samples=name)
    match_df.to_csv(sys.stdout, sep='\t', index=False, float_format='%0.4g', na_rep='-')


@query.command(name='sample-kmer')
@click.pass_context
@click.argument('name', nargs=-1)
def query_sample_kmer(ctx, name: List[str]):
    kmer_query_service = ctx.obj['data_index_connection'].kmer_query_service
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
    mutation_query_service = ctx.obj['data_index_connection'].mutation_query_service

    features = [QueryFeatureMutation(n) for n in name]
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
    mlst_query_service = ctx.obj['data_index_connection'].mlst_query_service

    features = [QueryFeatureMLST(n) for n in name]
    if not summarize:
        match_df = mlst_query_service.find_by_features(features, include_unknown=include_unknown)
    else:
        match_df = mlst_query_service.count_by_features(features, include_unknown=include_unknown)

    match_df.to_csv(sys.stdout, sep='\t', index=False, float_format='%0.4g', na_rep='-')


@main.group(name='db')
@click.pass_context
def db(ctx):
    project = Project(ctx.obj['project_dir'])
    ctx.obj['data_index_connection'] = project.create_connection()


UNITS = ['B', 'KB', 'MB', 'GB']


@db.command(name='size')
@click.pass_context
@click.option('--unit', default='B', help='The unit to display data sizes as.', type=click.Choice(UNITS))
def db_size(ctx, unit):
    size_df = ctx.obj['data_index_connection'].db_size(unit)
    size_df.to_csv(sys.stdout, sep='\t', index=False, float_format='%0.2f', na_rep='-')
