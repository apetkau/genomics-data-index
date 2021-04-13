from typing import List, Dict, Set, Union

from storage.variant.SampleSet import SampleSet
from storage.variant.model import Sample, Reference, ReferenceSequence, NucleotideVariantsSamples
from storage.variant.model import SampleMLSTAlleles, MLSTScheme, MLSTAllelesSamples
from storage.variant.service import DatabaseConnection
from storage.variant.service.QueryService import QueryFeature


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
            .join(Sample.sample_nucleotide_variation) \
            .filter(Reference.name == reference_name) \
            .all()
        return samples

    def get_samples_with_mlst_alleles(self, scheme_name: str) -> List[Sample]:
        """
        Gets a list of all samples that have MLST alleles associated with the given scheme name.
        :scheme_name: The scheme name.
        :return: A list of Samples with MLST alleles with respect to the scheme name, empty list of no Samples.
        """
        samples = self._connection.get_session().query(Sample) \
            .join(Sample.sample_mlst_alleles) \
            .join(SampleMLSTAlleles.scheme) \
            .filter(MLSTScheme.name == scheme_name) \
            .all()
        return samples

    def get_samples_with_variants_on_sequence(self, sequence_name: str) -> List[Sample]:
        """
        Gets a list of all samples that have variants associated with the given sequence name.
        :sequence_name: The sequence name.
        :return: A list of Samples with variants with respect to the sequence name, empty list of no Samples.
        """
        samples = self._connection.get_session().query(Sample) \
            .join(Sample.sample_nucleotide_variation) \
            .join(Reference.sequences) \
            .filter(ReferenceSequence.sequence_name == sequence_name) \
            .all()
        return samples

    def get_samples_associated_with_reference(self, reference_name: str) -> List[Sample]:
        """
        Gets a list of all samples associated with a reference name.
        :reference_name: The reference name.
        :return: A list of Samples associated with the reference name or an empty list if no Samples.
        """
        samples = self._connection.get_session().query(Sample) \
            .join(Sample.sample_nucleotide_variation) \
            .filter(Reference.name == reference_name) \
            .all()
        return samples

    def count_samples_associated_with_reference(self, reference_name: str) -> int:
        return len(self.get_samples_associated_with_reference(reference_name))

    def count_samples_associated_with_mlst_scheme(self, scheme_name: str) -> int:
        return len(self.get_samples_with_mlst_alleles(scheme_name))

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

    def get_sample(self, sample_name: str) -> Sample:
        return self._connection.get_session().query(Sample) \
            .filter(Sample.name == sample_name) \
            .one()

    def exists(self, sample_name: str):
        return self._connection.get_session().query(Sample) \
                   .filter(Sample.name == sample_name).count() > 0

    def find_samples_by_ids(self, sample_ids: Union[List[int], SampleSet]) -> List[Sample]:
        if isinstance(sample_ids, SampleSet):
            sample_ids = list(sample_ids)

        return self._connection.get_session().query(Sample) \
            .filter(Sample.id.in_(sample_ids)) \
            .all()

    def _get_variants_samples_by_variation_ids(self, variation_ids: List[str]) -> List[NucleotideVariantsSamples]:
        variants = self._connection.get_session().query(NucleotideVariantsSamples) \
            .filter(NucleotideVariantsSamples._spdi.in_(variation_ids)) \
            .all()

        return variants

    def _get_mlst_samples_by_mlst_allele_ids(self, mlst_allele_ids: List[str]) -> List[MLSTAllelesSamples]:
        mlst_alleles = self._connection.get_session().query(MLSTAllelesSamples) \
            .filter(MLSTAllelesSamples._sla.in_(mlst_allele_ids)) \
            .all()

        return mlst_alleles

    def _get_feature_type(self, features: List[QueryFeature]) -> str:
        feature_types = {type(f).__name__ for f in features}

        if len(feature_types) != 1:
            raise Exception(f'Should only be one feature type but instead got: {feature_types}.')
        else:
            return feature_types.pop()

    def find_samples_by_features(self, features: List[QueryFeature]) -> Dict[str, List[Sample]]:
        feature_type = self._get_feature_type(features)

        feature_ids = list({f.id for f in features})

        if feature_type == 'QueryFeatureMutation':
            variants = self._get_variants_samples_by_variation_ids(feature_ids)

            return {v.spdi: self.find_samples_by_ids(v.sample_ids) for v in variants}
        elif feature_type == 'QueryFeatureMLST':
            mlst_alleles = self._get_mlst_samples_by_mlst_allele_ids(feature_ids)

            return {a.sla: self.find_samples_by_ids(a.sample_ids) for a in mlst_alleles}
        else:
            raise Exception(f'Invalid feature type {feature_type}')

    def count_samples_by_features(self, features: List[QueryFeature]) -> Dict[str, List[Sample]]:
        feature_type = self._get_feature_type(features)

        feature_ids = list({f.id for f in features})

        if feature_type == 'QueryFeatureMutation':
            variants = self._get_variants_samples_by_variation_ids(feature_ids)

            return {v.spdi: len(v.sample_ids) for v in variants}
        elif feature_type == 'QueryFeatureMLST':
            mlst_alleles = self._get_mlst_samples_by_mlst_allele_ids(feature_ids)

            allele_id_to_count =  {a.sla: len(a.sample_ids) for a in mlst_alleles}
            for f in features:
                if f.id not in allele_id_to_count:
                    allele_id_to_count[f.id] = 0
            return allele_id_to_count
        else:
            raise Exception(f'Invalid feature type {feature_type}')

    def find_sample_name_ids(self, sample_names: Set[str]) -> Dict[str, int]:
        """
        Given a list of sample names, returns a dictionary mapping the sample names to sample IDs.
        :param sample_names: The sample names to search.
        :return: A dictionary linking the sample names to IDs.
        """
        sample_tuples = self._connection.get_session().query(Sample.name, Sample.id) \
            .filter(Sample.name.in_(sample_names)) \
            .all()

        return dict(sample_tuples)
