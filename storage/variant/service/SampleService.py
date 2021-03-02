from typing import List, Dict

from storage.variant.model import Sample, Reference, ReferenceSequence, VariationAllele, SampleSequence
from storage.variant.service import DatabaseConnection


class SampleService:

    def __init__(self, database_connection: DatabaseConnection):
        self._connection = database_connection

    def get_samples_with_variants(self, reference_name: str) -> List[Sample]:
        """
        Gets a list of all samples that have variants associated with the given reference genome name.
        :reference_name: The reference genome name.
        :return: A list of Samples with variants with respect to the reference genome name, empty list of no Samples.
        """
        samples = self._connection.get_session().query(Sample) \
            .join(Sample.variants) \
            .join(ReferenceSequence) \
            .join(Reference) \
            .filter(Reference.name == reference_name) \
            .all()
        return samples

    def get_samples_with_variants_on_sequence(self, sequence_name: str) -> List[Sample]:
        """
        Gets a list of all samples that have variants associated with the given sequence name.
        :sequence_name: The sequence name.
        :return: A list of Samples with variants with respect to the sequence name, empty list of no Samples.
        """
        samples = self._connection.get_session().query(Sample) \
            .join(Sample.variants) \
            .join(ReferenceSequence) \
            .filter(ReferenceSequence.sequence_name == sequence_name) \
            .all()
        return samples

    def get_samples_associated_with_sequence(self, sequence_name: str) -> List[Sample]:
        """
        Gets a list of all samples associated with a sequence name (whether they have variants or not).
        :sequence_name: The sequence name.
        :return: A list of Samples associated with the sequence name, empty list of no Samples.
        """
        samples = self._connection.get_session().query(Sample) \
            .join(Sample.sample_sequences) \
            .join(ReferenceSequence) \
            .filter(ReferenceSequence.sequence_name == sequence_name) \
            .all()
        return samples

    def count_samples_associated_with_sequence(self, sequence_name: str) -> int:
        return len(self.get_samples_associated_with_sequence(sequence_name))

    def get_samples(self) -> List[Sample]:
        return self._connection.get_session().query(Sample).all()

    def which_exists(self, sample_names: List[str]) -> List[str]:
        """
        Returns which of the given samples exist in the database.
        :param sample_names: The list of sample names.
        :return: A list of those passed sample names that exist in the database.
        """
        samples = self._connection.get_session().query(Sample) \
            .filter(Sample.name.in_(sample_names)) \
            .all()
        return [sample.name for sample in samples]

    def find_samples_by_variation_ids(self, variation_ids: List[str]) -> Dict[str, List[Sample]]:
        variants = self._connection.get_session().query(VariationAllele) \
            .filter(VariationAllele.id.in_(variation_ids)) \
            .all()

        return {v.id: v.samples for v in variants}
