from __future__ import annotations

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
from genomics_data_index.storage.util import execute_commands, TRACE_LEVEL

logger = logging.getLogger(__name__)

template_dir = Path(path.dirname(__file__)) / 'templates'

jinja_env = Environment(
    loader=FileSystemLoader(template_dir),
    autoescape=select_autoescape()
)


class SequenceFile:

    def __init__(self, file: Path):
        self._file = file

    @property
    def file(self):
        return self._file

    def parse_sequence_file(self) -> Tuple[str, List[SeqRecord]]:
        # Code for handling gzipped/non-gzipped from https://stackoverflow.com/a/52839332
        _open = partial(gzip.open, mode='rt') if self.is_gzip() else open

        ref_name = self.get_genome_name()

        with _open(self._file) as f:
            if self.is_genbank():
                logger.log(TRACE_LEVEL, f'Parsing sequence file [{self._file}] as a genbank file')
                sequences = list(SeqIO.parse(f, 'genbank'))
            else:
                logger.log(TRACE_LEVEL, f'Parsing sequence file [{self._file}] as a fasta file')
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

        logger.log(TRACE_LEVEL, f'Sequence file [{self._file}] has extension (minus compression) [{ref_extension}]')

        return ref_extension

    def is_gzip(self) -> bool:
        encoding = guess_type(str(self._file))[1]  # uses file extension
        return encoding == 'gzip'

    def is_genbank(self) -> bool:
        ref_extension = self.get_genome_extension_minus_compression().lower()
        return ref_extension == '.gb' or ref_extension == '.gbk'

    def is_fasta(self) -> bool:
        ref_extension = self.get_genome_extension_minus_compression().lower()
        return ref_extension == '.fasta' or ref_extension == '.fa' or ref_extension == '.fna'

    def is_assembly(self) -> bool:
        return self.is_fasta() or self.is_genbank()

    def is_reads(self) -> bool:
        extension = self.get_genome_extension_minus_compression().lower()
        return extension == '.fastq' or extension == '.fq'

    def name_differences(self, other_file: SequenceFile) -> List[str]:
        """
        Gets the differences from this sequence file name and the passed file name.
        I assume both file names are of the same length (i.e., no insertion of gaps anywhere or
        complicated alignment of strings). This is intended to be used to figure out which of a pair
        of fastq files is the forward and which is the reverse.
        :param other_file: The other sequence file.
        :return: A list of strings containing mismatches between this file name and the other file name.
        """
        name1 = self._file.name
        name2 = other_file._file.name

        if len(name1) != len(name2):
            raise Exception(f'name1=[{name1}] is not the same length as name2=[{name2}]')
        else:
            differences = []
            current_difference = None
            previous_difference_index = -1
            for i in range(len(name1)):
                c1 = name1[i]
                c2 = name2[i]
                if c1 == c2:
                    if current_difference is not None:
                        differences.append(current_difference)
                        current_difference = None
                else:
                    if previous_difference_index == (i - 1):
                        if current_difference is None:
                            current_difference = c1
                            previous_difference_index = i
                        else:
                            current_difference += c1
                            previous_difference_index = i
                    elif current_difference is None:
                        current_difference = c1
                        previous_difference_index = i
                    else:
                        raise Exception(f'previous_difference_index=[{previous_difference_index}], i=[{i}], '
                                        f'but current_difference is not None')

            # Append last bit of differences
            if current_difference is not None:
                differences.append(current_difference)

            return differences

    def can_use_snpeff(self) -> bool:
        return self.is_genbank()

    def get_genome_name(self, exclude_paired_end_indicators: bool = False) -> str:
        """
        Gets the genome name (filename minus extension). Accounts for gzipped/non-gzipped files.
        :param file: The file.
        :param exclude_paired_end_indicators: Exclude paired-end indicators when getting genome name (e.g., _R1, _1, etc).
                                              this only applies if the sequence file is a reads file.
        :return: The genome name from the file.
        """
        if self.is_gzip():
            name_minus_extension = splitext(basename(self._file).rstrip('.gz'))[0]
        else:
            name_minus_extension = splitext(basename(self._file))[0]

        if self.is_reads() and exclude_paired_end_indicators:
            name_minus_pair_indicator = name_minus_extension
            for pair_indicator in ['_1', '_2', '_R1', '_R2', '_f', '_r']:
                last_index = name_minus_extension.rfind(pair_indicator)
                if last_index != -1:
                    name_minus_pair_indicator = name_minus_extension[:last_index]
                    break
            return name_minus_pair_indicator
        else:
            return name_minus_extension

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
