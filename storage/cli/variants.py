from os import path, listdir
from pathlib import Path
from typing import List

import click
from Bio import AlignIO

from storage.variant.VariantsReader import SnippyVariantsReader
from storage.variant.service import DatabaseConnection
from storage.variant.service.CoreAlignmentService import CoreAlignmentService
from storage.variant.service.ReferenceService import ReferenceService
from storage.variant.service.VariationService import VariationService


@click.group()
@click.pass_context
@click.option('--database-connection', help='A connection string for the database.')
@click.option('--seqrepo-dir', help='The root directory for the seqrepo reference storage.',
              type=click.Path())
def main(ctx, database_connection, seqrepo_dir):
    ctx.ensure_object(dict)

    click.echo(f'Connecting to database {database_connection}')
    database = DatabaseConnection(database_connection)

    click.echo(f'Use seqrepo directory {seqrepo_dir}')
    reference_service = ReferenceService(database, seqrepo_dir)

    variation_service = VariationService(database, reference_service)

    ctx.obj['database'] = database
    ctx.obj['reference_service'] = reference_service
    ctx.obj['variation_service'] = variation_service


@main.command()
@click.pass_context
@click.argument('snippy_dir', type=click.Path(exists=True))
@click.option('--reference-file', help='Reference genome', type=click.Path(exists=True))
def load(ctx, snippy_dir: Path, reference_file: Path):
    snippy_dir = Path(snippy_dir)
    reference_file = Path(reference_file)
    click.echo(f'Loading {snippy_dir}')
    sample_dirs = [snippy_dir / d for d in listdir(snippy_dir) if path.isdir(snippy_dir / d)]
    variants_reader = SnippyVariantsReader(sample_dirs)

    reference_service = ctx.obj['reference_service']
    variation_service = ctx.obj['variation_service']

    reference_service.add_reference_genome(reference_file)
    var_df = variants_reader.get_variants_table()
    core_masks = variants_reader.get_core_masks()

    variation_service.insert_variants(var_df=var_df,
                                      reference_name='genome',
                                      core_masks=core_masks)
    click.echo(f'Loaded variants from [{snippy_dir}] into database')


@main.command()
@click.pass_context
@click.option('--output-file', help='Output file', type=click.Path())
@click.option('--reference-name', help='Reference genome name', type=str)
@click.option('--sample', help='Sample to include in alignment (can list more than one).', multiple=True, type=str)
def alignment(ctx, output_file: Path, reference_name: str, sample: List[str]):
    database = ctx.obj['database']
    reference_service = ctx.obj['reference_service']

    alignment_service = CoreAlignmentService(database, reference_service)

    alignment_data = alignment_service.construct_alignment(reference_name=reference_name,
                                                           samples=sample,
                                                           include_reference=True)

    with open(output_file, 'w') as f:
        AlignIO.write(alignment_data, f, 'fasta')