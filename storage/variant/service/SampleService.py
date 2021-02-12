from typing import List

from storage.variant.model import Sample, Reference, ReferenceSequence
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
