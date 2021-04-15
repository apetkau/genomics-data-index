from pathlib import Path
from typing import List, Dict

import ga4gh.vrs.dataproxy as dataproxy
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biocommons.seqrepo import SeqRepo
from ete3 import Tree

from storage.variant.model.db import Reference, SampleNucleotideVariation, ReferenceSequence, Sample
from storage.variant.service import DatabaseConnection, EntityExistsError
from storage.variant.util import parse_sequence_file


class ReferenceService:

    def __init__(self, database_connection: DatabaseConnection, seq_repo_dir: Path,
                 seq_repo_namespace: str = 'storage'):
        self._connection = database_connection
        self._seq_repo_namespace = seq_repo_namespace

        if seq_repo_dir is not None:
            self._seq_repo_updatable = SeqRepo(seq_repo_dir, writeable=True)
            self._seq_repo_proxy = dataproxy.SeqRepoDataProxy(SeqRepo(seq_repo_dir))

    def add_reference_genome(self, genome_file: Path):
        (genome_name, sequences) = parse_sequence_file(genome_file)

        if self.exists_reference_genome(genome_name):
            raise EntityExistsError(f'Reference genome [{genome_name}] already exists in database')
        else:
            for record in sequences:
                self._seq_repo_updatable.store(str(record.seq),
                                               [{'namespace': self._seq_repo_namespace, 'alias': record.id}])
            self._seq_repo_updatable.commit()

            self._create_reference_genome_db(genome_file)

    def get_reference_sequences(self, reference_name: str) -> Dict[str, ReferenceSequence]:
        reference = self.find_reference_genome(reference_name)
        return {s.sequence_name: s for s in reference.sequences}

    def get_sequence(self, sequence_name: str) -> SeqRecord:
        namespace = self._seq_repo_namespace
        seq_string = self._seq_repo_proxy.get_sequence(f'{namespace}:{sequence_name}')
        return SeqRecord(Seq(seq_string), id=sequence_name)

    def get_reference_genome_records(self, reference_name: str) -> List[SeqRecord]:
        reference = self.find_reference_genome(reference_name)
        return [self.get_sequence(sequence.sequence_name) for sequence in reference.sequences]

    def _create_reference_genome_db(self, reference_file: Path):
        ref_length = 0
        ref_contigs = {}

        (ref_name, sequences) = parse_sequence_file(reference_file)
        for record in sequences:
            ref_contigs[record.id] = ReferenceSequence(
                sequence_name=record.id, sequence_length=len(record.seq))
            ref_length += len(record.seq)

        reference = Reference(name=ref_name, length=ref_length, sequences=list(ref_contigs.values()))

        self._connection.get_session().add(reference)
        self._connection.get_session().commit()

    def update_tree(self, reference_name: str, tree: Tree, alignment_length: int):
        if alignment_length is None or alignment_length <= 0:
            raise Exception(f'Invalid alignment_length=[{alignment_length}]')

        reference = self.find_reference_genome(reference_name)
        reference.tree = tree
        reference.tree_alignment_length = alignment_length
        self._connection.get_session().commit()

    def find_reference_genome(self, name: str):
        return self._connection.get_session().query(Reference).filter_by(name=name).one()

    def exists_reference_genome(self, name: str):
        return self._connection.get_session().query(Reference.id).filter_by(name=name).scalar() is not None

    def get_reference_genomes(self) -> List[Reference]:
        return self._connection.get_session().query(Reference).all()

    def find_references_for_sample(self, sample_name: str) -> List[Reference]:
        return self._connection.get_session().query(Reference) \
            .join(SampleNucleotideVariation) \
            .join(Sample) \
            .filter(Sample.name == sample_name) \
            .all()

    def find_reference_for_sequence(self, sequence_name: str) -> Reference:
        return self._connection.get_session().query(Reference) \
            .join(Reference.sequences) \
            .filter(ReferenceSequence.sequence_name == sequence_name) \
            .one()
