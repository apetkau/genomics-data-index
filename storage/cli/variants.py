import click
from os import path, listdir
from pathlib import Path

from storage.variant.service import DatabaseConnection
from storage.variant.service.ReferenceService import ReferenceService
from storage.variant.service.VariationService import VariationService
from storage.variant.VariantsReader import SnippyVariantsReader

database = DatabaseConnection('sqlite:///:memory:')
reference_service = None


@click.group()
@click.option('--database-connection', help='A connection string for the database.')
@click.option('--seqrepo-dir', help='The root directory for the seqrepo reference storage.',
              type=click.Path())
def main(database_connection, seqrepo_dir):
    global database
    global reference_service

    click.echo(f'Connecting to database {database_connection}')
    database = DatabaseConnection(database_connection)

    click.echo(f'Use seqrepo directory {seqrepo_dir}')
    reference_service = ReferenceService(database, seqrepo_dir)


@main.command()
@click.argument('snippy_dir', type=click.Path(exists=True))
@click.option('--reference-file', help='Reference genome', type=click.Path(exists=True))
def load(snippy_dir, reference_file):
    snippy_dir = Path(snippy_dir)
    reference_file = Path(reference_file)
    click.echo(f'Loading {snippy_dir}')
    sample_dirs = [snippy_dir / d for d in listdir(snippy_dir) if path.isdir(snippy_dir / d)]
    variants_reader = SnippyVariantsReader(sample_dirs)

    reference_service.add_reference_genome(reference_file)
    var_df = variants_reader.get_variants_table()
    core_masks = variants_reader.get_core_masks()
    print(database)

    var_service = VariationService(database, reference_service)
    var_service.insert_variants(var_df=var_df,
                                reference_name='genome',
                                core_masks=core_masks)


