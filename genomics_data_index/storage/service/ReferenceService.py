from pathlib import Path
from typing import List, Dict, Iterable

import ga4gh.vrs.dataproxy as dataproxy
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biocommons.seqrepo import SeqRepo
from ete3 import Tree
from sqlalchemy.exc import NoResultFound

from genomics_data_index.storage.io.mutation.SequenceFile import SequenceFile
from genomics_data_index.storage.model.QueryFeatureMutationSPDI import QueryFeatureMutationSPDI
from genomics_data_index.storage.model.db import Reference, SampleNucleotideVariation, ReferenceSequence, Sample
from genomics_data_index.storage.service import DatabaseConnection, EntityExistsError


class ReferenceService:
    MUTATION_ID_TYPES = ['spdi_ref']

    def __init__(self, database_connection: DatabaseConnection, seq_repo_dir: Path,
                 seq_repo_namespace: str = 'genomics_data_index'):
        self._connection = database_connection
        self._seq_repo_namespace = seq_repo_namespace

        if seq_repo_dir is not None:
            self._seq_repo_updatable = SeqRepo(seq_repo_dir, writeable=True)
            self._seq_repo_proxy = dataproxy.SeqRepoDataProxy(SeqRepo(seq_repo_dir))

    def add_reference_genome(self, genome_file: Path):
        (genome_name, sequences) = SequenceFile(genome_file).parse_sequence_file()

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

        (ref_name, sequences) = SequenceFile(reference_file).parse_sequence_file()
        for record in sequences:
            ref_contigs[record.id] = ReferenceSequence(
                sequence_name=record.id, sequence_length=len(record.seq))
            ref_length += len(record.seq)

        reference = Reference(name=ref_name, length=ref_length, sequences=list(ref_contigs.values()))

        self._connection.get_session().add(reference)
        self._connection.get_session().commit()

    def update_tree(self, reference_name: str, tree: Tree, alignment_length: int) -> None:
        if alignment_length is None or alignment_length <= 0:
            raise Exception(f'Invalid alignment_length=[{alignment_length}]')

        reference = self.find_reference_genome(reference_name)
        reference.tree = tree
        reference.tree_alignment_length = alignment_length
        self._connection.get_session().commit()

    def find_reference_genome(self, name: str) -> Reference:
        try:
            return self._connection.get_session().query(Reference).filter_by(name=name).one()
        except NoResultFound as e:
            raise EntityExistsError(f'No reference genome with name=[{name}]')

    def exists_reference_genome(self, name: str) -> bool:
        return self._connection.get_session().query(Reference.id).filter_by(name=name).scalar() is not None

    def get_reference_genomes(self) -> List[Reference]:
        return self._connection.get_session().query(Reference).all()

    def count_reference_genomes(self) -> int:
        return self._connection.get_session().query(Reference).count()

    def translate_spdi(self, spdi_ids: Iterable, to: str = 'spdi_ref') -> Dict[str, str]:
        if to not in self.MUTATION_ID_TYPES:
            raise Exception(f'to={to} must be one of {self.MUTATION_ID_TYPES}')

        # Convert to query feature mutations so I can get spdi parts
        spdi_features = {QueryFeatureMutationSPDI(i) for i in spdi_ids}

        # Get set of sequence names since it's likely all (or most) of the input spdi_ids are on the same sequence
        sequence_names = {f.sequence for f in spdi_features}

        # Create map of seq name to sequences so I don't have to re-generate them for every spdi_id
        namespace = self._seq_repo_namespace
        name_sequence_map = {n: self._seq_repo_proxy.get_sequence(f'{namespace}:{n}') for n in sequence_names}

        spdi_ids_map = {}
        for f in spdi_features:
            reference_sequence_str = name_sequence_map[f.sequence]
            spdi_del_str = reference_sequence_str[f.start0:f.stop0]
            spdi_ref = f'{f.sequence}:{f.position}:{spdi_del_str}:{f.insertion}'
            spdi_ids_map[f.id] = spdi_ref

        return spdi_ids_map

    def find_references_for_sample(self, sample_name: str) -> List[Reference]:
        return self._connection.get_session().query(Reference) \
            .join(Reference.sample_nucleotide_variation) \
            .join(SampleNucleotideVariation.sample) \
            .filter(Sample.name == sample_name) \
            .all()

    def find_reference_for_sequence(self, sequence_name: str) -> Reference:
        return self._connection.get_session().query(Reference) \
            .join(Reference.sequences) \
            .filter(ReferenceSequence.sequence_name == sequence_name) \
            .one()
