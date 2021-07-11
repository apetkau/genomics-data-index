import gzip
import logging
from functools import partial
from mimetypes import guess_type
from os import mkdir, symlink, path
from os.path import basename, splitext
from pathlib import Path
from typing import Tuple, List

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from jinja2 import Environment, FileSystemLoader, select_autoescape

from genomics_data_index.storage.io.mutation.SnpEffDatabase import SnpEffDatabase
from genomics_data_index.storage.util import execute_commands

logger = logging.getLogger(__name__)

template_dir = Path(path.dirname(__file__)) / 'templates'

jinja_env = Environment(
    loader=FileSystemLoader(template_dir),
    autoescape=select_autoescape()
)


class SequenceFile:

    def __init__(self, file: Path):
        self._file = file

    def parse_sequence_file(self) -> Tuple[str, List[SeqRecord]]:
        # Code for handling gzipped/non-gzipped from https://stackoverflow.com/a/52839332
        _open = partial(gzip.open, mode='rt') if self.is_gzip() else open

        ref_name = self.get_genome_name()

        with _open(self._file) as f:
            if self.is_genbank():
                logger.debug(f'Parsing sequence file [{self._file}] as a genbank file')
                sequences = list(SeqIO.parse(f, 'genbank'))
            else:
                logger.debug(f'Parsing sequence file [{self._file}] as a fasta file')
                sequences = list(SeqIO.parse(f, 'fasta'))

            return ref_name, sequences

    def write(self, output_file: Path, format: str) -> None:
        ref_name, sequences = self.parse_sequence_file()
        with open(output_file, 'w') as fh:
            SeqIO.write(sequences, fh, format)

    def get_genome_extension_minus_compression(self):
        if self.is_gzip():
            ref_extension = splitext(basename(self._file).rstrip('.gz'))[1]
        else:
            ref_extension = splitext(basename(self._file))[1]

        logger.debug(f'Sequence file [{self._file}] has extension (minus compression) [{ref_extension}]')

        return ref_extension

    def is_gzip(self) -> bool:
        encoding = guess_type(str(self._file))[1]  # uses file extension
        return encoding == 'gzip'

    def is_genbank(self) -> bool:
        ref_extension = self.get_genome_extension_minus_compression().lower()
        return ref_extension == '.gb' or ref_extension == '.gbk'

    def can_use_snpeff(self) -> bool:
        return self.is_genbank()

    def get_genome_name(self) -> str:
        """
        Gets the genome name (filename minus extension). Accounts for gzipped/non-gzipped files.
        :param file: The file.
        :return: The genome name from the file.
        """
        if self.is_gzip():
            ref_name = splitext(basename(self._file).rstrip('.gz'))[0]
        else:
            ref_name = splitext(basename(self._file))[0]

        return ref_name

    def _write_snpeff_config(self, snpeff_config_file: Path, snpeff_database_dir: Path,
                             codon_type: str) -> str:
        genome_name, records = self.parse_sequence_file()
        sequence_ids = [s.id for s in records]

        snpeff_config_template = jinja_env.get_template('snpEff.config')
        with open(snpeff_config_file, 'w') as snpeff_config:
            snpeff_config.write(snpeff_config_template.render(
                data_dir=str(snpeff_database_dir),
                reference_name=genome_name,
                sequence_ids=sequence_ids,
                codon_type=codon_type,
            ))

        logger.debug(f'Wrote snpeff config file [{snpeff_config_file}]')

        return genome_name

    def _setup_snpeff_files(self, snpeff_database_dir: Path, reference_name: str):
        if not snpeff_database_dir.exists():
            mkdir(snpeff_database_dir)

        reference_dir = snpeff_database_dir / reference_name
        mkdir(reference_dir)

        if self.is_gzip():
            genbank_path = reference_dir / 'genes.gbk.gz'
        else:
            genbank_path = reference_dir / 'genes.gbk'
        symlink(self._file, genbank_path)

    def create_snpeff_database(self, database_dir: Path, codon_type: str = 'Standard') -> SnpEffDatabase:
        if not self.is_genbank():
            raise Exception(f'Sequence file [{self._file}] is not a genbank file. '
                            f'Can only build snpeff databases for genbank files.')

        logger.debug(f'Setting up snpeff database in [{database_dir}]')

        snpeff_config_path = database_dir / 'snpEff.config'
        snpeff_database_dir = database_dir / 'db'

        genome_name = self._write_snpeff_config(
            snpeff_config_file=snpeff_config_path,
            snpeff_database_dir=snpeff_database_dir,
            codon_type=codon_type)

        self._setup_snpeff_files(snpeff_database_dir=snpeff_database_dir,
                                 reference_name=genome_name)

        command = ['snpEff', 'build', '-v', '-c', str(snpeff_config_path), '-genbank', genome_name]
        execute_commands([command])

        return SnpEffDatabase(snpeff_config=snpeff_config_path,
                              database_dir=database_dir,
                              genome_name=genome_name)
