import gzip
from functools import partial
from mimetypes import guess_type
from os.path import basename, splitext
from pathlib import Path
from typing import Tuple, List
import logging
from os import path
import shutil

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from jinja2 import Environment, PackageLoader, select_autoescape


logger = logging.getLogger(__name__)

jinja_env = Environment(
    loader=PackageLoader(__name__),
    autoescape=select_autoescape()
)


class SequenceFile:

    def __init__(self, file: Path):
        self._file = file

    def parse_sequence_file(self) -> Tuple[str, List[SeqRecord]]:
        # Code for handling gzipped/non-gzipped from https://stackoverflow.com/a/52839332
        encoding = guess_type(str(self._file))[1]  # uses file extension
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

        ref_name = self.get_genome_name()
        ref_extension = self.get_genome_extension_minus_compression().lower()

        logger.debug(f'Sequence file [{self._file}] has extension (minus compression) [{ref_extension}]')
        with _open(self._file) as f:
            if ref_extension == '.gb' or ref_extension == '.gbk':
                logger.debug(f'Parsing sequence file [{self._file}] as a genbank file')
                sequences = list(SeqIO.parse(f, 'genbank'))
            else:
                logger.debug(f'Parsing sequence file [{self._file}] as a fasta file')
                sequences = list(SeqIO.parse(f, 'fasta'))

            return ref_name, sequences

    def get_genome_extension_minus_compression(self):
        encoding = guess_type(str(self._file))[1]

        if encoding == 'gzip':
            ref_extension = splitext(basename(self._file).rstrip('.gz'))[1]
        else:
            ref_extension = splitext(basename(self._file))[1]

        return ref_extension

    def get_genome_name(self) -> str:
        """
        Gets the genome name (filename minus extension). Accounts for gzipped/non-gzipped files.
        :param file: The file.
        :return: The genome name from the file.
        """
        encoding = guess_type(str(self._file))[1]

        if encoding == 'gzip':
            ref_name = splitext(basename(self._file).rstrip('.gz'))[0]
        else:
            ref_name = splitext(basename(self._file))[0]

        return ref_name

    def create_snpeff_database(self, database_dir: Path) -> Path:
        logger.debug(f'Setting up snpeff database in [{database_dir}]')

        snpeff_config_path = database_dir / 'snpEff.config'
        snpeff_database_dir = database_dir / 'db'

        snpeff_config_template = jinja_env.get_template('snpEff.config')
        with open(snpeff_config_path, 'w') as snpeff_config:
            snpeff_config.write(snpeff_config_template.render(data_dir=str(snpeff_database_dir)))

        logger.debug(f'Writing snpeff config file [{snpeff_config_path}]')

        return snpeff_config_path
