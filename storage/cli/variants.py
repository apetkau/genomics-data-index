from os import path, listdir
from pathlib import Path
from typing import List
import logging
import sys
import multiprocessing

import click
import click_config_file
import coloredlogs
from Bio import AlignIO

from storage.cli import yaml_config_provider
from storage.variant.io.SnippyVariantsReader import SnippyVariantsReader
from storage.variant.service import DatabaseConnection, EntityExistsError
from storage.variant.service.CoreAlignmentService import CoreAlignmentService
from storage.variant.service.ReferenceService import ReferenceService
from storage.variant.service.SampleSequenceService import SampleSequenceService
from storage.variant.service.SampleService import SampleService
from storage.variant.service.TreeService import TreeService
from storage.variant.service.VariationService import VariationService
from storage.variant.service.SampleQueryService import SampleQueryService
from storage.variant.util import get_genome_name

logger = logging.getLogger('storage')
num_cores = multiprocessing.cpu_count()


@click.group()
@click.pass_context
@click.option('--database-connection', help='A connection string for the database.')
@click.option('--seqrepo-dir', help='The root directory for the seqrepo reference storage.',
              type=click.Path())
@click.option('--verbose/--no-verbose', default=False, help='Turn up verbosity of command')
@click_config_file.configuration_option(provider=yaml_config_provider,
                                        config_file_name='config.yaml',
                                        implicit=True)
def main(ctx, database_connection, seqrepo_dir, verbose):
    ctx.ensure_object(dict)

    if verbose:
        coloredlogs.install(level='DEBUG',
                            fmt='%(asctime)s %(levelname)s %(name)s.%(funcName)s,%(lineno)s: %(message)s',
                            logger=logger)
    else:
        coloredlogs.install(level='WARNING', fmt='%(asctime)s %(levelname)s: %(message)s',
                            logger=logger)

    logger.info(f'Connecting to database {database_connection}')
    database = DatabaseConnection(database_connection)

    logger.info(f'Use seqrepo directory {seqrepo_dir}')
    reference_service = ReferenceService(database, seqrepo_dir)

    sample_service = SampleService(database)
    sample_sequence_service = SampleSequenceService(database)
    variation_service = VariationService(database_connection=database,
                                         reference_service=reference_service,
                                         sample_service=sample_service)
    alignment_service = CoreAlignmentService(database=database,
                                             reference_service=reference_service,
                                             variation_service=variation_service,
                                             sample_sequence_service=sample_sequence_service)
    tree_service = TreeService(database, reference_service, alignment_service)
    sample_query_service = SampleQueryService(tree_service=tree_service,
                                              reference_service=reference_service)

    ctx.obj['database'] = database
    ctx.obj['reference_service'] = reference_service
    ctx.obj['variation_service'] = variation_service
    ctx.obj['alignment_service'] = alignment_service
    ctx.obj['tree_service'] = tree_service
    ctx.obj['sample_service'] = sample_service
    ctx.obj['sample_query_service'] = sample_query_service


@main.command()
@click.pass_context
@click.argument('snippy_dir', type=click.Path(exists=True))
@click.option('--reference-file', help='Reference genome', required=True, type=click.Path(exists=True))
@click.option('--build-tree/--no-build-tree', default=False, help='Builds tree of all samples after loading')
@click.option('--threads', help='Threads for building tree', default=1,
              type=click.IntRange(min=1, max=num_cores))
def load(ctx, snippy_dir: Path, reference_file: Path, build_tree: bool, threads: int):
    snippy_dir = Path(snippy_dir)
    reference_file = Path(reference_file)
    click.echo(f'Loading {snippy_dir}')
    sample_dirs = [snippy_dir / d for d in listdir(snippy_dir) if path.isdir(snippy_dir / d)]
    variants_reader = SnippyVariantsReader(sample_dirs)

    reference_service = ctx.obj['reference_service']
    variation_service = ctx.obj['variation_service']
    sample_service = ctx.obj['sample_service']
    tree_service = ctx.obj['tree_service']

    try:
        reference_service.add_reference_genome(reference_file)
    except EntityExistsError as e:
        logger.warning(f'Reference genome [{reference_file}] already exists, will not load')

    samples_exist = sample_service.which_exists(variants_reader.samples_list())
    if len(samples_exist) > 0:
        logger.error(f'Samples {samples_exist} already exist, will not load any variants')
    else:
        var_df = variants_reader.get_variants_table()
        core_masks = variants_reader.get_core_masks()

        reference_name = get_genome_name(reference_file)

        variation_service.insert_variants(var_df=var_df,
                                          reference_name=reference_name,
                                          core_masks=core_masks)
        click.echo(f'Loaded variants from [{snippy_dir}] into database')

        if build_tree:
            tree_service.rebuild_tree(reference_name=reference_name,
                                      num_cores=threads)
            click.echo('Finished building tree of all samples')


LIST_TYPES = ['genomes', 'samples']


@main.command(name='list')
@click.pass_context
@click.option('--type', 'data_type', required=True, help='Type of data to list',
              type=click.Choice(LIST_TYPES))
def list_data(ctx, data_type):
    if data_type == 'genomes':
        items = {genome.name for genome in ctx.obj['reference_service'].get_reference_genomes()}
    elif data_type == 'samples':
        items = {sample.name for sample in ctx.obj['sample_service'].get_samples()}
    else:
        raise Exception(f'Unknown data_type=[{data_type}]')

    click.echo('\n'.join(items))


@main.command()
@click.pass_context
@click.option('--output-file', help='Output file', required=True, type=click.Path())
@click.option('--reference-name', help='Reference genome name', required=True, type=str)
@click.option('--align-type', help=f'The type of alignment to generate', default='core',
              type=click.Choice(CoreAlignmentService.ALIGN_TYPES))
@click.option('--sample', help='Sample to include in alignment (can list more than one).',
              multiple=True, type=str)
def alignment(ctx, output_file: Path, reference_name: str, align_type: str, sample: List[str]):
    alignment_service = ctx.obj['alignment_service']
    reference_service = ctx.obj['reference_service']
    sample_service = ctx.obj['sample_service']

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


@main.command()
@click.pass_context
@click.option('--output-file', help='Output file', required=True, type=click.Path())
@click.option('--reference-name', help='Reference genome name', type=str, required=True)
@click.option('--align-type', help=f'The type of alignment to use for generating the tree', default='core',
              type=click.Choice(CoreAlignmentService.ALIGN_TYPES))
@click.option('--tree-build-type', help=f'The type of tree building software', default='iqtree',
              type=click.Choice(TreeService.TREE_BUILD_TYPES))
@click.option('--sample', help='Sample to include in tree (can list more than one).',
              multiple=True, type=str)
@click.option('--threads', help='Threads for building tree', default=1,
              type=click.IntRange(min=1, max=num_cores))
@click.option('--extra-params', help='Extra parameters to tree-building software',
              default=None)
def tree(ctx, output_file: Path, reference_name: str, align_type: str,
         tree_build_type: str, sample: List[str], threads: int, extra_params: str):
    alignment_service = ctx.obj['alignment_service']
    tree_service = ctx.obj['tree_service']
    reference_service = ctx.obj['reference_service']
    sample_service = ctx.obj['sample_service']

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
                                             num_cores=threads, align_type=align_type, extra_params=extra_params)
    tree_data.write(outfile=output_file)
    click.echo(f'Wrote tree to [{output_file}]')
    with open(log_file, 'w') as log:
        log.write(out)
        click.echo(f'Wrote log file to [{log_file}]')


QUERY_TYPES = ['sample']


@main.command()
@click.pass_context
@click.argument('name')
@click.option('--type', 'query_type', help='Query type',
              required=True, type=click.Choice(QUERY_TYPES))
def query(ctx, name, query_type):
    sample_query_service = ctx.obj['sample_query_service']

    match_df = None
    if query_type == 'sample':
        match_df = sample_query_service.find_matches(sample_name=name)
    else:
        logger.error(f'Invalid query_type=[{query_type}]')
        sys.exit(1)

    match_df.to_csv(sys.stdout, sep='\t', index=False)
