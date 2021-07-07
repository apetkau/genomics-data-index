import gzip
from functools import partial
from mimetypes import guess_type
from os.path import basename, splitext
from pathlib import Path
from typing import Tuple, List

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class SequenceFile:

    def __init__(self, file: Path):
        self._file = file

    def parse_sequence_file(self) -> Tuple[str, List[SeqRecord]]:
        # Code for handling gzipped/non-gzipped from https://stackoverflow.com/a/52839332
        encoding = guess_type(str(self._file))[1]  # uses file extension
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

        ref_name = self.get_genome_name()

        with _open(self._file) as f:
            sequences = list(SeqIO.parse(f, 'fasta'))

            return ref_name, sequences

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
