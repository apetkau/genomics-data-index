import logging
import multiprocessing
import sys
from os import path, listdir
from pathlib import Path
from typing import List

import click
import click_config_file
import coloredlogs
import pandas as pd
from Bio import AlignIO

from storage.FilesystemStorage import FilesystemStorage
from storage.cli import yaml_config_provider
from storage.variant.index.KmerIndexer import KmerIndexerSourmash, KmerIndexManager
from storage.variant.io.SnippyVariantsReader import SnippyVariantsReader
from storage.variant.io.VcfVariantsReader import VcfVariantsReader
from storage.variant.service import DatabaseConnection, EntityExistsError
from storage.variant.service.CoreAlignmentService import CoreAlignmentService
from storage.variant.service.KmerService import KmerService
from storage.variant.service.MutationQueryService import MutationQueryService, QueryFeatureMutation
from storage.variant.service.KmerQueryService import KmerQueryService
from storage.variant.service.ReferenceService import ReferenceService
from storage.variant.service.SampleService import SampleService
from storage.variant.service.TreeService import TreeService
from storage.variant.service.VariationService import VariationService
from storage.variant.util import get_genome_name

logger = logging.getLogger('storage')
num_cores = multiprocessing.cpu_count()


@click.group()
@click.pass_context
@click.option('--database-connection', help='A connection string for the database.')
@click.option('--database-dir', help='The root directory for the database files.',
              type=click.Path())
@click.option('--verbose/--no-verbose', default=False, help='Turn up verbosity of command')
@click_config_file.configuration_option(provider=yaml_config_provider,
                                        config_file_name='config.yaml',
                                        implicit=True)
def main(ctx, database_connection, database_dir, verbose):
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
    filesystem_storage = FilesystemStorage(Path(database_dir))

    logger.info(f'Use database directory {database_dir}')
    reference_service = ReferenceService(database, filesystem_storage.reference_dir)

    sample_service = SampleService(database)
    variation_service = VariationService(database_connection=database,
                                         variation_dir=filesystem_storage.variation_dir,
                                         reference_service=reference_service,
                                         sample_service=sample_service)
    alignment_service = CoreAlignmentService(database=database,
                                             reference_service=reference_service,
                                             sample_service=sample_service,
                                             variation_service=variation_service)
    tree_service = TreeService(database, reference_service, alignment_service)
    mutation_query_service = MutationQueryService(reference_service=reference_service,
                                                  sample_service=sample_service,
                                                  tree_service=tree_service)

    kmer_service = KmerService(database_connection=database,
                               sample_service=sample_service)

    kmer_query_service = KmerQueryService(sample_service=sample_service)

    ctx.obj['database'] = database
    ctx.obj['filesystem_storage'] = filesystem_storage
    ctx.obj['kmer_service'] = kmer_service
    ctx.obj['reference_service'] = reference_service
    ctx.obj['variation_service'] = variation_service
    ctx.obj['alignment_service'] = alignment_service
    ctx.obj['tree_service'] = tree_service
    ctx.obj['sample_service'] = sample_service
    ctx.obj['mutation_query_service'] = mutation_query_service
    ctx.obj['kmer_query_service'] = kmer_query_service


@main.group()
@click.pass_context
def load(ctx):
    pass


def load_variants_common(ctx, variants_reader, reference_file, input, build_tree, align_type, threads,
                         extra_tree_params: str):
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
        reference_name = get_genome_name(reference_file)

        variation_service.insert_variants(reference_name=reference_name,
                                          variants_reader=variants_reader)
        click.echo(f'Loaded variants from [{input}] into database')

        if build_tree:
            tree_service.rebuild_tree(reference_name=reference_name,
                                      align_type=align_type,
                                      num_cores=threads,
                                      extra_params=extra_tree_params)
            click.echo('Finished building tree of all samples')


@load.command(name='snippy')
@click.pass_context
@click.argument('snippy_dir', type=click.Path(exists=True))
@click.option('--reference-file', help='Reference genome', required=True, type=click.Path(exists=True))
@click.option('--build-tree/--no-build-tree', default=False, help='Builds tree of all samples after loading')
@click.option('--align-type', help=f'The type of alignment to generate', default='core',
              type=click.Choice(CoreAlignmentService.ALIGN_TYPES))
@click.option('--threads', help='Threads for building tree', default=1,
              type=click.IntRange(min=1, max=num_cores))
@click.option('--extra-tree-params', help='Extra parameters to tree-building software',
              default=None)
def load_snippy(ctx, snippy_dir: Path, reference_file: Path, build_tree: bool, align_type: str, threads: int,
                extra_tree_params: str):
    snippy_dir = Path(snippy_dir)
    reference_file = Path(reference_file)
    click.echo(f'Loading {snippy_dir}')
    sample_dirs = [snippy_dir / d for d in listdir(snippy_dir) if path.isdir(snippy_dir / d)]
    variants_reader = SnippyVariantsReader(sample_dirs)

    load_variants_common(ctx=ctx, variants_reader=variants_reader, reference_file=reference_file,
                         input=snippy_dir, build_tree=build_tree, align_type=align_type, threads=threads,
                         extra_tree_params=extra_tree_params)


@load.command(name='vcf')
@click.pass_context
@click.argument('vcf_fofns', type=click.Path(exists=True))
@click.option('--reference-file', help='Reference genome', required=True, type=click.Path(exists=True))
@click.option('--build-tree/--no-build-tree', default=False, help='Builds tree of all samples after loading')
@click.option('--align-type', help=f'The type of alignment to generate', default='core',
              type=click.Choice(CoreAlignmentService.ALIGN_TYPES))
@click.option('--threads', help='Threads for building tree', default=1,
              type=click.IntRange(min=1, max=num_cores))
@click.option('--extra-tree-params', help='Extra parameters to tree-building software',
              default=None)
def load_vcf(ctx, vcf_fofns: Path, reference_file: Path, build_tree: bool, align_type: str, threads: int,
             extra_tree_params: str):
    reference_file = Path(reference_file)

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

    variants_reader = VcfVariantsReader(sample_vcf_map=sample_vcf_map,
                                        masked_genomic_files_map=mask_files_map)

    load_variants_common(ctx=ctx, variants_reader=variants_reader, reference_file=reference_file,
                         input=Path(vcf_fofns), build_tree=build_tree, align_type=align_type, threads=threads,
                         extra_tree_params=extra_tree_params)


@load.command(name='kmer')
@click.pass_context
@click.argument('kmer_fofns', type=click.Path(exists=True))
@click.option('--kmer-size', help='Kmer size for indexing. List multiple for multiple kmer sizes in an index',
              default=[31], multiple=True, type=click.IntRange(min=1, max=201))
def load_kmer(ctx, kmer_fofns, kmer_size):
    filesystem_storage = ctx.obj['filesystem_storage']
    kmer_service = ctx.obj['kmer_service']

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


@main.group(name='list')
@click.pass_context
def list_data(ctx):
    pass


@list_data.command(name='genomes')
@click.pass_context
def list_genomes(ctx):
    items = [genome.name for genome in ctx.obj['reference_service'].get_reference_genomes()]
    click.echo('\n'.join(items))


@list_data.command(name='samples')
@click.pass_context
def list_samples(ctx):
    items = [sample.name for sample in ctx.obj['sample_service'].get_samples()]
    click.echo('\n'.join(items))


@main.group()
@click.pass_context
def export(ctx):
    pass


@export.command(name='tree')
@click.pass_context
@click.argument('name', nargs=-1)
@click.option('--ascii/--no-ascii', help='Export as ASCII figure')
def export_tree(ctx, name: List[str], ascii: bool):
    if len(name) == 0:
        logger.warning('No reference genome names passed, will not export tree')

    for ref_name in name:
        reference = ctx.obj['reference_service'].find_reference_genome(ref_name)
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


@main.group()
@click.pass_context
def rebuild(ctx):
    pass


@rebuild.command(name='tree')
@click.pass_context
@click.argument('reference', type=str, nargs=-1)
@click.option('--align-type', help=f'The type of alignment to use for generating the tree', default='core',
              type=click.Choice(CoreAlignmentService.ALIGN_TYPES))
@click.option('--threads', help='Threads for building tree', default=1,
              type=click.IntRange(min=1, max=num_cores))
@click.option('--extra-params', help='Extra parameters to tree-building software',
              default=None)
def rebuild_tree(ctx, reference: List[str], align_type: str, threads: int, extra_params: str):
    tree_service = ctx.obj['tree_service']
    reference_service = ctx.obj['reference_service']

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
                                  num_cores=threads,
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
    mutation_query_service = ctx.obj['mutation_query_service']
    match_df = mutation_query_service.find_matches(samples=name)
    match_df.to_csv(sys.stdout, sep='\t', index=False, float_format='%0.4g', na_rep='-')


@query.command(name='sample-kmer')
@click.pass_context
@click.argument('name', nargs=-1)
def query_sample_kmer(ctx, name: List[str]):
    kmer_query_service = ctx.obj['kmer_query_service']
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
    mutation_query_service = ctx.obj['mutation_query_service']

    features = [QueryFeatureMutation(n) for n in name]
    if not summarize:
        match_df = mutation_query_service.find_by_features(features, include_unknown=include_unknown)
    else:
        match_df = mutation_query_service.count_by_features(features, include_unknown=include_unknown)

    match_df.to_csv(sys.stdout, sep='\t', index=False, float_format='%0.4g', na_rep='-')


@main.group(name='db')
@click.pass_context
def db(ctx):
    pass


@db.command(name='size')
@click.pass_context
@click.option('--kb', help='Print in KB', is_flag=True)
@click.option('--mb', help='Print in MB', is_flag=True)
@click.option('--gb', help='Print in GB', is_flag=True)
def db_size(ctx, kb, mb, gb):
    database = ctx.obj['database']
    filesystem_storage = ctx.obj['filesystem_storage']

    factor = 1
    unit = 'B'
    if kb:
        factor = 1024
        unit = 'KB'
    elif mb:
        factor = 1024 ** 2
        unit = 'MB'
    elif gb:
        factor = 1024 ** 3
        unit = 'GB'

    filesystem_df = filesystem_storage.get_storage_size()
    database_df = database.get_database_size()

    size_df = pd.concat([filesystem_df, database_df])
    total_data_size = size_df['Data Size'].sum()
    total_index_size = size_df['Index Size'].sum()
    total_items = size_df['Number of Items'].sum()
    total_row = pd.DataFrame([['Total', pd.NA, pd.NA, total_data_size, total_index_size, total_items]],
                             columns=['Type', 'Name', 'Division', 'Data Size', 'Index Size', 'Number of Items'])
    size_df = pd.concat([size_df, total_row])

    # Reorder columns
    size_df = size_df[['Type', 'Name', 'Division', 'Data Size', 'Index Size', 'Number of Items']]

    size_df['Data Size'] = size_df['Data Size'] / factor
    size_df['Index Size'] = size_df['Index Size'] / factor
    size_df = size_df.rename({'Data Size': f'Data Size ({unit})',
                    'Index Size': f'Index Size ({unit})'}, axis='columns')

    size_df.to_csv(sys.stdout, sep='\t', index=False, float_format='%0.2f', na_rep='-')