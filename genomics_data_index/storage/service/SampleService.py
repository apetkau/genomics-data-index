import logging
from typing import List, Dict, Set, Union, cast, Type

import pandas as pd

from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.NucleotideMutationTranslater import NucleotideMutationTranslater
from genomics_data_index.storage.model.QueryFeature import QueryFeature
from genomics_data_index.storage.model.QueryFeatureHGVS import QueryFeatureHGVS
from genomics_data_index.storage.model.QueryFeatureHGVSGN import QueryFeatureHGVSGN
from genomics_data_index.storage.model.QueryFeatureMLST import QueryFeatureMLST
from genomics_data_index.storage.model.QueryFeatureMutation import QueryFeatureMutation
from genomics_data_index.storage.model.QueryFeatureMutationSPDI import QueryFeatureMutationSPDI
from genomics_data_index.storage.model.db import NucleotideVariantsSamples, Reference, ReferenceSequence, MLSTScheme, \
    SampleMLSTAlleles, MLSTAllelesSamples, Sample
from genomics_data_index.storage.model.db import SampleNucleotideVariation
from genomics_data_index.storage.service import DatabaseConnection

logger = logging.getLogger(__name__)


class FeatureExplodeUnknownError(Exception):

    def __init__(self, msg: str):
        super().__init__(msg)


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

    def feature_explode_unknown(self, feature: QueryFeature) -> List[QueryFeature]:
        if isinstance(feature, QueryFeatureHGVSGN):
            features_spdi = self.find_features_spdi_for_hgvsgn(feature)

            if len(features_spdi) == 0:
                raise FeatureExplodeUnknownError(f'feature={feature} is of type HGVSGN but the corresponding SPDI '
                                                 f'feature does not exist in the database. Cannot convert to unknown '
                                                 f'SPDI representation.')
            else:
                unknown_features = []
                for feature in features_spdi:
                    unknown_features.extend(feature.to_unknown_explode())
                return unknown_features
        elif isinstance(feature, QueryFeatureHGVS):
            if feature.is_nucleotide():
                variants_hgvs = self._connection.get_session().query(NucleotideVariantsSamples) \
                    .filter(NucleotideVariantsSamples._id_hgvs_c == feature.id) \
                    .all()
            elif feature.is_protein():
                variants_hgvs = self._connection.get_session().query(NucleotideVariantsSamples) \
                    .filter(NucleotideVariantsSamples._id_hgvs_p == feature.id) \
                    .all()
            else:
                raise Exception(f'feature=[{feature}] is neither nucleotide or protein')

            if len(variants_hgvs) == 0:
                raise FeatureExplodeUnknownError(f'feature={feature} is of type HGVS but the corresponding SPDI '
                                                 f'feature does not exist in the database. Cannot convert to unknown '
                                                 f'SPDI representation.')
            else:
                unknown_features = []
                for variants_sample_obj in variants_hgvs:
                    unknown_features.extend(QueryFeatureMutationSPDI(variants_sample_obj.spdi).to_unknown_explode())

                return unknown_features
        else:
            return feature.to_unknown_explode()

    def find_features_spdi_for_hgvsgn(self, feature: QueryFeatureHGVSGN) -> List[QueryFeatureMutationSPDI]:
        if not isinstance(feature, QueryFeatureHGVSGN):
            raise Exception(f'Cannot handle feature={feature}. Not of type {QueryFeatureHGVSGN.__name__}')

        query = self._connection.get_session().query(NucleotideVariantsSamples).filter(
            NucleotideVariantsSamples.sequence == feature.sequence)

        if feature.has_gene():
            query = query.filter(NucleotideVariantsSamples.annotation_gene_name == feature.gene)

        if feature.is_nucleotide():
            query = query.filter(NucleotideVariantsSamples.annotation_hgvs_c == feature.variant)
        elif feature.is_protein():
            query = query.filter(NucleotideVariantsSamples.annotation_hgvs_p == feature.variant)
        else:
            raise Exception(f'feature={feature} is neither protein nor nucleotide')

        return [QueryFeatureMutationSPDI(s.spdi) for s in query.all()]

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

    def create_dataframe_from_sample_set(self, present_set: SampleSet,
                                         absent_set: SampleSet,
                                         unknown_set: SampleSet,
                                         queries_expression: str) -> pd.DataFrame:
        sample_sets_status_list = [(present_set, 'Present'), (absent_set, 'Absent'), (unknown_set, 'Unknown')]
        data = []
        for sample_status in sample_sets_status_list:
            sample_set = sample_status[0]
            status = sample_status[1]

            if not sample_set.is_empty():
                samples = self.find_samples_by_ids(sample_set)
                for sample in samples:
                    data.append([queries_expression, sample.name, sample.id, status])

        return pd.DataFrame(data=data, columns=['Query', 'Sample Name', 'Sample ID', 'Status'])

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
            elif isinstance(feature, QueryFeatureHGVSGN):
                logger.warning(f'feature=[{feature}] is a QueryFeatureHGVSGN and I do not handle it here.')
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
        feature_ids = list({f.id_no_prefix for f in features})
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

    def find_unknown_sample_sets_by_features(self, features: List[QueryFeature]) -> Dict[str, SampleSet]:
        unknown_to_features_dict = {}
        unknown_features = []
        for feature in features:
            try:
                unknown_features_exploded = self.feature_explode_unknown(feature)
                unknown_features.extend(unknown_features_exploded)
                for unknown_feature in unknown_features_exploded:
                    unknown_to_features_dict[unknown_feature.id] = feature
            except FeatureExplodeUnknownError as e:
                logger.warning(
                    f'Could not map feature={feature} to a set of unknown features. Will assume no unknowns exist.')

        if len(unknown_features) > 0:
            unknown_features_sets = self.find_sample_sets_by_features(unknown_features)
        else:
            unknown_features_sets = set()

        features_to_unknown_sample_sets = {}
        for uid in unknown_features_sets:
            fid = unknown_to_features_dict[uid].id
            sample_set = unknown_features_sets[uid]

            # If we've already set this sample set with the same feature,
            # We need to merge together the unknown sample sets
            # This can occur if, e.g., we have a large deletion and are iterating over each
            # Base in the deletion in turn (e.g., ref:10:ATT:A -> gets converted to
            # ['ref:10:A:?', 'ref:11:T:?', 'ref:12:T:?'], we need to merge unknown sample results
            # for each of these features in turn.
            if fid in features_to_unknown_sample_sets:
                previous_sample_set = features_to_unknown_sample_sets[fid]
                features_to_unknown_sample_sets[fid] = previous_sample_set.union(sample_set)
            else:
                features_to_unknown_sample_sets[fid] = sample_set

        return features_to_unknown_sample_sets

    def find_sample_sets_by_features(self, features: List[QueryFeature]) -> Dict[str, SampleSet]:
        feature_type = self._get_feature_type(features)

        if issubclass(feature_type, QueryFeatureHGVSGN):
            # In this case where I'm querying by gene name, first convert to SPDI features before lookup
            # TODO: it's not the most efficient to do this as a loop, but it's easier to implement right now
            hgvs_gn_id_to_sampleset = dict()
            for feature in features:
                feature = cast(QueryFeatureHGVSGN, feature)
                features_spdi = self.find_features_spdi_for_hgvsgn(feature)
                variants_dict = self.get_variants_samples_by_variation_features(features_spdi)
                variants_nuc_variants_samples = list(variants_dict.values())

                if len(variants_nuc_variants_samples) == 0:
                    samples_union = SampleSet.create_empty()
                else:
                    first_nuc_variant_samples = variants_nuc_variants_samples.pop()
                    samples_union = first_nuc_variant_samples.sample_ids

                    # Handle remaining, if any
                    for nuc_variant_samples in variants_nuc_variants_samples:
                        samples_union = samples_union.union(nuc_variant_samples.sample_ids)

                hgvs_gn_id_to_sampleset[feature.id] = samples_union
            return hgvs_gn_id_to_sampleset
        elif issubclass(feature_type, QueryFeatureMutation):
            features = cast(List[QueryFeatureMutation], features)
            variants_dict = self.get_variants_samples_by_variation_features(features)

            return {id: variants_dict[id].sample_ids for id in variants_dict}
        elif issubclass(feature_type, QueryFeatureMLST):
            features = cast(List[QueryFeatureMLST], features)
            mlst_alleles = self._get_mlst_samples_by_mlst_features(features)

            return {a.query_id: a.sample_ids for a in mlst_alleles}
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

            return {a.query_id: self.find_samples_by_ids(a.sample_ids) for a in mlst_alleles}
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

            allele_id_to_count = {a.query_id: len(a.sample_ids) for a in mlst_alleles}
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

    def get_sample_set_by_names(self, sample_names: Union[List[str], Set[str]],
                                ignore_not_found: bool = False) -> SampleSet:
        """
        Given a collection of sample names, get a SampleSet of the corresponding IDs.
        :param sample_names: The names to convert to an ID set.
        :param ignore_not_found: Whether or not to ignore sample names that were not found.
        :return: A SampleSet with all the corresponding samples by the passed names. If ignore_not_found is false,
                 raises an exception if some sample names have no ids.
        """
        if isinstance(sample_names, list):
            sample_names = set(sample_names)
        elif not isinstance(sample_names, set):
            raise Exception(f'Invalid type=[{type(sample_names)}] for passed sample_names. Must be list or set.')

        sample_ids_tuples = self._connection.get_session().query(Sample.id) \
            .filter(Sample.name.in_(sample_names)) \
            .all()
        sample_ids = {i for i, in sample_ids_tuples}
        sample_set = SampleSet(sample_ids=sample_ids)

        if ignore_not_found or len(sample_names) == len(sample_set):
            return sample_set
        else:
            # Find matching sample names to ids we did find for a nicer error message
            found_sample_names = {s.name for s in self.find_samples_by_ids(sample_set)}
            names_not_found = sample_names - found_sample_names
            if len(names_not_found) > 10:
                small_not_found = list(names_not_found)[:10]
                msg = f'[{", ".join(small_not_found)}, ...]'
            else:
                msg = f'{names_not_found}'

            raise Exception(f'Did not find an equal number of sample names and ids. '
                            f'Number sample_names={len(sample_names)}. Number returned sample_ids={len(sample_ids)}. '
                            f'Sample names with missing ids {msg}')
