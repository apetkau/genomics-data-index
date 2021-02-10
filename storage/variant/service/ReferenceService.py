import gzip
from functools import partial
from mimetypes import guess_type
from pathlib import Path
from typing import Tuple, List

import ga4gh.vrs.dataproxy as dataproxy
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from biocommons.seqrepo import SeqRepo

from storage.variant.model import Reference
from storage.variant.model import ReferenceSequence
from storage.variant.service import DatabaseConnection, EntityExistsError
from storage.variant.util import get_genome_name


class ReferenceService:

    def __init__(self, database_connection: DatabaseConnection, seq_repo_dir: Path,
                 seq_repo_namespace: str = 'storage'):
        self._connection = database_connection
        self._seq_repo_namespace = seq_repo_namespace

        if seq_repo_dir is not None:
            self._seq_repo_updatable = SeqRepo(seq_repo_dir, writeable=True)
            self._seq_repo_proxy = dataproxy.SeqRepoDataProxy(SeqRepo(seq_repo_dir))

    def _parse_sequence_file(self, sequence_file: Path) -> Tuple[str, List[SeqRecord]]:
        # Code for handling gzipped/non-gzipped from https://stackoverflow.com/a/52839332
        encoding = guess_type(str(sequence_file))[1]  # uses file extension
        _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

        ref_name = get_genome_name(sequence_file)

        with _open(sequence_file) as f:
            sequences = list(SeqIO.parse(f, 'fasta'))

            return (ref_name, sequences)

    def add_reference_genome(self, genome_file: Path):
        (genome_name, sequences) = self._parse_sequence_file(genome_file)

        if self.exists_reference_genome(genome_name):
            raise EntityExistsError(f'Reference genome [{genome_name}] already exists in database')
        else:
            for record in sequences:
                self._seq_repo_updatable.store(str(record.seq),
                                               [{'namespace': self._seq_repo_namespace, 'alias': record.id}])
            self._seq_repo_updatable.commit()

            self._create_reference_genome_db(genome_file)

    def get_reference_contigs(self, reference_name: str):
        reference = self.find_reference_genome(reference_name)
        return {s.sequence_name: s for s in reference.sequences}

    def get_sequence(self, sequence_name: str):
        namespace = self._seq_repo_namespace
        seq_string = self._seq_repo_proxy.get_sequence(f'{namespace}:{sequence_name}')
        return SeqRecord(seq_string, id=sequence_name)

    def _create_reference_genome_db(self, reference_file: Path):
        ref_length = 0
        ref_contigs = {}

        (ref_name, sequences) = self._parse_sequence_file(reference_file)
        for record in sequences:
            ref_contigs[record.id] = ReferenceSequence(
                sequence_name=record.id, sequence_length=len(record.seq))
            ref_length += len(record.seq)

        reference = Reference(name=ref_name, length=ref_length, sequences=list(ref_contigs.values()))

        self._connection.get_session().add(reference)
        self._connection.get_session().commit()

    def find_reference_genome(self, name: str):
        return self._connection.get_session().query(Reference).filter_by(name=name).one()

    def exists_reference_genome(self, name: str):
        return self._connection.get_session().query(Reference.id).filter_by(name=name).scalar() is not None
