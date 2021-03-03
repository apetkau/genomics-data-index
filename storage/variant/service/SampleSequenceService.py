from typing import List, Dict

from storage.variant.model import Sample, SampleSequence, Reference, ReferenceSequence
from storage.variant.service import DatabaseConnection


class SampleSequenceService:

    def __init__(self, database_connection: DatabaseConnection):
        self._connection = database_connection

    def get_sample_sequences(self, reference_name: str, samples: List[str]) -> Dict[str, List[SampleSequence]]:
        sample_sequences = self._connection.get_session().query(SampleSequence) \
            .join(Sample) \
            .join(ReferenceSequence) \
            .join(Reference) \
            .filter(Sample.name.in_(samples)) \
            .filter(Reference.name == reference_name) \
            .all()

        sample_sequences_map = {}
        for ss in sample_sequences:
            if ss.sequence.sequence_name in sample_sequences_map:
                sample_sequences_map[ss.sequence.sequence_name].append(ss)
            else:
                sample_sequences_map[ss.sequence.sequence_name] = [ss]

        return sample_sequences_map

    def missing_in_sequence(self, sample_name: str, sequence_name: str, positions: List[int]) -> bool:
        """
        Checks if any of the given positions are missing in the sequence data.
        :param positions: The list of positions.
        :return: True if any positions are missing, False otherwise.
        """
        if positions is None or len(positions) == 0:
            return False

        sample_sequence = self._connection.get_session() \
            .query(SampleSequence) \
            .join(Sample) \
            .join(ReferenceSequence) \
            .filter(Sample.name == sample_name) \
            .filter(ReferenceSequence.sequence_name == sequence_name) \
            .one()

        core_mask = sample_sequence.core_mask
        return not all([pos in core_mask for pos in positions])
