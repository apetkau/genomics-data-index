from typing import List, Dict, Set, Union, cast, Type

import pandas as pd

from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.NucleotideMutationTranslater import NucleotideMutationTranslater
from genomics_data_index.storage.model.QueryFeature import QueryFeature
from genomics_data_index.storage.model.QueryFeatureHGVS import QueryFeatureHGVS
from genomics_data_index.storage.model.QueryFeatureMLST import QueryFeatureMLST
from genomics_data_index.storage.model.QueryFeatureMutation import QueryFeatureMutation
from genomics_data_index.storage.model.QueryFeatureMutationSPDI import QueryFeatureMutationSPDI
from genomics_data_index.storage.model.db import NucleotideVariantsSamples, Reference, ReferenceSequence, MLSTScheme, \
    SampleMLSTAlleles, MLSTAllelesSamples, Sample
from genomics_data_index.storage.model.db import SampleNucleotideVariation
from genomics_data_index.storage.service import DatabaseConnection


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
            .join(SampleNucleotideVariation.reference) \
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
            .join(SampleNucleotideVariation.reference) \
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
            .join(SampleNucleotideVariation.reference) \
            .filter(Reference.name == reference_name) \
            .all()
        return samples

    def get_samples_set_associated_with_reference(self, reference_name: str) -> SampleSet:
        """
        Gets a list of all samples associated with a reference name.
        :reference_name: The reference name.
        :return: A list of Samples associated with the reference name or an empty list if no Samples.
        """
        sample_ids = [i for i, in self._connection.get_session().query(Sample.id) \
            .join(Sample.sample_nucleotide_variation) \
            .join(SampleNucleotideVariation.reference) \
            .filter(Reference.name == reference_name) \
            .all()]

        return SampleSet(sample_ids=sample_ids)

    def create_dataframe_from_sample_set(self, sample_set: SampleSet,
                                         universe_set: SampleSet,
                                         exclude_absent: bool,
                                         queries_expression: str) -> pd.DataFrame:
        samples = self.find_samples_by_ids(sample_set)
        data = []
        for sample in samples:
            data.append([queries_expression, sample.name, sample.id, 'Present'])

        if not exclude_absent:
            complement_samples_set = self.find_samples_by_ids(universe_set.minus(sample_set))
            for sample in complement_samples_set:
                data.append([queries_expression, sample.name, sample.id, 'Absent'])

        results_df = pd.DataFrame(data=data, columns=['Query', 'Sample Name', 'Sample ID', 'Status'])
        return results_df

    def count_samples_associated_with_reference(self, reference_name: str) -> int:
        return self._connection.get_session().query(Sample) \
            .join(Sample.sample_nucleotide_variation) \
            .join(SampleNucleotideVariation.reference) \
            .filter(Reference.name == reference_name) \
            .count()

    def count_samples_associated_with_mlst_scheme(self, scheme_name: str) -> int:
        return len(self.get_samples_with_mlst_alleles(scheme_name))

    def get_samples(self) -> List[Sample]:
        return self._connection.get_session().query(Sample).all()

    def count_samples(self) -> int:
        return self._connection.get_session().query(Sample).count()

    def get_all_sample_ids(self) -> SampleSet:
        ids_list = [id for id, in self._connection.get_session().query(Sample.id).all()]
        return SampleSet(ids_list)

    def get_existing_samples_by_names(self, sample_names: List[str]) -> List[Sample]:
        return self._connection.get_session().query(Sample) \
            .filter(Sample.name.in_(sample_names)) \
            .all()

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

    def get_variants_samples_by_variation_features(self, features: List[QueryFeatureMutation]) -> Dict[
        str, NucleotideVariantsSamples]:
        standardized_features_to_input_feature = {}
        standardized_features_ids = set()
        standardized_feature_hgvs_c_ids = set()
        standardized_feature_hgvs_p_ids = set()
        for feature in features:
            if isinstance(feature, QueryFeatureMutationSPDI):
                dbf = NucleotideMutationTranslater.to_db_feature(feature)
                if dbf.id in standardized_features_to_input_feature:
                    standardized_features_to_input_feature[dbf.id].append(feature.id)
                else:
                    standardized_features_to_input_feature[dbf.id] = [feature.id]
                standardized_features_ids.add(dbf.id)
            elif isinstance(feature, QueryFeatureHGVS):
                if feature.is_nucleotide():
                    standardized_feature_hgvs_c_ids.add(feature.id)
                elif feature.is_protein():
                    standardized_feature_hgvs_p_ids.add(feature.id)
                else:
                    raise Exception(f'feature=[{feature}] is neither nucleotide or protein')
            else:
                raise Exception(f'Invalid type for feature=[{feature}]. '
                                f'Must be either {QueryFeatureMutationSPDI.__class__.__name__} or '
                                f'{QueryFeatureHGVS.__class__.__name__}')

        if len(standardized_features_ids) > 0:
            variants_spdi = self._connection.get_session().query(NucleotideVariantsSamples) \
                .filter(NucleotideVariantsSamples._spdi.in_(standardized_features_ids)) \
                .all()
        else:
            variants_spdi = []

        if len(standardized_feature_hgvs_c_ids) > 0:
            variants_hgvs_c = self._connection.get_session().query(NucleotideVariantsSamples) \
                .filter(NucleotideVariantsSamples._id_hgvs_c.in_(standardized_feature_hgvs_c_ids)) \
                .all()
        else:
            variants_hgvs_c = []

        if len(standardized_feature_hgvs_p_ids) > 0:
            variants_hgvs_p = self._connection.get_session().query(NucleotideVariantsSamples) \
                .filter(NucleotideVariantsSamples._id_hgvs_p.in_(standardized_feature_hgvs_p_ids)) \
                .all()
        else:
            variants_hgvs_p = []

        # Map back unstandardized IDs to the actual variant object
        # Use this because some features can have multiple identifiers for the same feature
        # (e.g., ref:10:A:T and ref:10:1:T). I want to make sure I map each passed id to the
        # same object (that is, in this example, I want to return a dictionary with two keys, one for each ID)
        unstandardized_variants = {}
        for v in variants_spdi:
            for vid in standardized_features_to_input_feature[v.spdi]:
                unstandardized_variants[vid] = v

        unstandardized_variants.update({v.id_hgvs_c: v for v in variants_hgvs_c})
        unstandardized_variants.update({v.id_hgvs_p: v for v in variants_hgvs_p})
        return unstandardized_variants

    def _get_mlst_samples_by_mlst_features(self, features: List[QueryFeatureMLST]) -> List[MLSTAllelesSamples]:
        feature_ids = list({f.id for f in features})
        mlst_alleles = self._connection.get_session().query(MLSTAllelesSamples) \
            .filter(MLSTAllelesSamples._sla.in_(feature_ids)) \
            .all()

        return mlst_alleles

    def _get_feature_type(self, features: List[QueryFeature]) -> Type[QueryFeature]:
        feature_types = {f.__class__ for f in features}

        if len(feature_types) != 1:
            raise Exception(f'Should only be one feature type but instead got: {feature_types}.')
        else:
            return feature_types.pop()

    def find_unknown_sample_sets_by_feature(self, features: List[QueryFeature]) -> Dict[str, SampleSet]:
        unknown_to_features_dict = {}
        unknown_features = []
        for feature in features:
            unknown_feature = feature.to_unknown()
            unknown_features.append(unknown_feature)
            unknown_to_features_dict[unknown_feature.id] = feature

        unknown_features_sets = self.find_sample_sets_by_features(unknown_features)

        return {unknown_to_features_dict[fid].id: unknown_features_sets[fid] for fid in unknown_features_sets}

    def find_sample_sets_by_features(self, features: List[QueryFeature]) -> Dict[str, SampleSet]:
        feature_type = self._get_feature_type(features)

        if issubclass(feature_type, QueryFeatureMutation):
            features = cast(List[QueryFeatureMutation], features)
            variants_dict = self.get_variants_samples_by_variation_features(features)

            return {id: variants_dict[id].sample_ids for id in variants_dict}
        elif issubclass(feature_type, QueryFeatureMLST):
            features = cast(List[QueryFeatureMLST], features)
            mlst_alleles = self._get_mlst_samples_by_mlst_features(features)

            return {a.sla: a.sample_ids for a in mlst_alleles}
        else:
            raise Exception(f'Invalid feature type {feature_type}')

    def find_samples_by_features(self, features: List[QueryFeature]) -> Dict[str, List[Sample]]:
        feature_type = self._get_feature_type(features)

        if issubclass(feature_type, QueryFeatureMutation):
            features = cast(List[QueryFeatureMutation], features)
            variants_dict = self.get_variants_samples_by_variation_features(features)

            return {id: self.find_samples_by_ids(variants_dict[id].sample_ids) for id in variants_dict}
        elif issubclass(feature_type, QueryFeatureMLST):
            features = cast(List[QueryFeatureMLST], features)
            mlst_alleles = self._get_mlst_samples_by_mlst_features(features)

            return {a.sla: self.find_samples_by_ids(a.sample_ids) for a in mlst_alleles}
        else:
            raise Exception(f'Invalid feature type {feature_type}')

    def count_samples_by_features(self, features: List[QueryFeature]) -> Dict[str, List[Sample]]:
        feature_type = self._get_feature_type(features)

        if issubclass(feature_type, QueryFeatureMutation):
            features = cast(List[QueryFeatureMutation], features)
            variants_dict = self.get_variants_samples_by_variation_features(features)

            return {id: len(variants_dict[id].sample_ids) for id in variants_dict}
        elif issubclass(feature_type, QueryFeatureMLST):
            features = cast(List[QueryFeatureMLST], features)
            mlst_alleles = self._get_mlst_samples_by_mlst_features(features)

            allele_id_to_count = {a.sla: len(a.sample_ids) for a in mlst_alleles}
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
