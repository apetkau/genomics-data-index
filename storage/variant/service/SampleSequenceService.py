from typing import List, Dict

from storage.variant.service import DatabaseConnection
from storage.variant.model import Sample, SampleSequence, Reference, ReferenceSequence


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
