import gzip
from functools import partial
from mimetypes import guess_type
from os.path import basename, splitext
from pathlib import Path

from Bio import SeqIO

from storage.variant.model import Reference
from storage.variant.model import ReferenceSequence
from storage.variant.service import DatabaseConnection


class ReferenceService:

    def __init__(self, database_connection: DatabaseConnection):
        self._connection = database_connection

    def create_reference_genome(self, reference_file: Path):
        ref_length = 0
        ref_contigs = {}
        ref_name = basename(reference_file)

        # Code for handling gzipped/non-gzipped from https://stackoverflow.com/a/52839332
        encoding = guess_type(str(reference_file))[1]  # uses file extension
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

        if encoding == 'gzip':
            ref_name = splitext(basename(reference_file).rstrip('.gz'))[0]
        else:
            ref_name = splitext(basename(reference_file))

        with _open(reference_file) as f:
            for record in SeqIO.parse(f, 'fasta'):
                ref_contigs[record.id] = ReferenceSequence(
                    sequence_name=record.id, sequence_length=len(record.seq))
                ref_length += len(record.seq)

        reference = Reference(name=ref_name, length=ref_length, sequences=list(ref_contigs.values()))

        self._connection.get_session().add(reference)
        self._connection.get_session().commit()

    def find_reference_genome(self, name: str):
        return self._connection.get_session().query(Reference).filter_by(name=name).one()
