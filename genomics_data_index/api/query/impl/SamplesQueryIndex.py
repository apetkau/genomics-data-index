from __future__ import annotations

import logging
from typing import Union, List, Set, Tuple, cast

import numpy as np
import pandas as pd
from ete3 import Tree

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.impl.DataFrameSamplesQuery import DataFrameSamplesQuery
from genomics_data_index.api.query.impl.QueriesCollection import QueriesCollection
from genomics_data_index.api.query.impl.TreeSamplesQueryFactory import TreeSamplesQueryFactory
from genomics_data_index.configuration.connector import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.QueryFeature import QueryFeature
from genomics_data_index.storage.model.QueryFeatureHGVS import QueryFeatureHGVS
from genomics_data_index.storage.model.QueryFeatureMLST import QueryFeatureMLST
from genomics_data_index.storage.model.QueryFeatureMutationSPDI import QueryFeatureMutationSPDI
from genomics_data_index.storage.service.KmerService import KmerService

logger = logging.getLogger(__name__)


class SamplesQueryIndex(SamplesQuery):
    """
    The main class implementing :py:class:`genomics_data_index.api.query.SamplesQuery`. This is used
    to store selections of sets of samples as well as operations to select subsets of samples.
    """

    HAS_KINDS = ['mutation', 'mutations', 'mlst']
    SUMMARY_FEATURES_KINDS = ['mutations']
    FEATURES_SELECTIONS = ['all', 'unique']
    ISIN_TYPES = ['sample', 'samples', 'distance', 'distances']
    ISA_TYPES = ['sample', 'samples']
    DISTANCES_UNITS = ['kmer_jaccard']

    def __init__(self, connection: DataIndexConnection,
                 universe_set: SampleSet,
                 sample_set: SampleSet,
                 unknown_set: SampleSet,
                 queries_collection: QueriesCollection = QueriesCollection.create_empty()):
        """
        Builds a new SamplesQueryIndex from the given information. In most normal operations SamplesQuery objects
        are not created directly but are instead created from an :py:class:`genomics_data_index.api.GenomicsDataIndex`
        object or from operations applied to a SamplesQuery.

        :param connection: A connection to a database containing samples.
        :param universe_set: The :py:class:`genomics_data_index.storage.SampleSet` representing a set of samples defining
                       the universe (used for e.g., complement() operations).
        :param sample_set: The :py:class:`genomics_data_index.storage.SampleSet` representing the set of selected samples.
        :param unknown_set: The :py:class:`genomics_data_index.storage.SampleSet` representing the set of unknown samples.
        :param queries_collection: A collection of strings representing the queries performed to arrive at this SamplesQuery.
        :return: A new SamplesQueryIndex object.
        """
        super().__init__()
        self._query_connection = connection
        self._universe_set = universe_set
        self._sample_set = sample_set
        self._unknown_set = unknown_set
        self._queries_collection = queries_collection

    @property
    def universe_set(self) -> SampleSet:
        return self._universe_set

    @property
    def sample_set(self) -> SampleSet:
        return self._sample_set

    @property
    def unknown_set(self) -> SampleSet:
        return self._unknown_set

    @property
    def _absent_set(self) -> SampleSet:
        return self._universe_set.minus([self._sample_set, self._unknown_set])

    @property
    def absent_set(self) -> SampleSet:
        return self._absent_set

    def reset_universe(self, include_unknown: bool = True) -> SamplesQuery:
        if include_unknown:
            universe_set = self._sample_set.union(self._unknown_set)
            unknown_set = self._unknown_set
        else:
            universe_set = self._sample_set
            unknown_set = SampleSet.create_empty()

        return self._create_from(sample_set=self._sample_set, universe_set=universe_set,
                                 unknown_set=unknown_set,
                                 queries_collection=self._queries_collection)

    def join(self, data_frame: pd.DataFrame, sample_ids_column: str = None,
             sample_names_column: str = None, default_isa_kind: str = 'sample',
             default_isa_column: str = None) -> SamplesQuery:
        if sample_ids_column is None and sample_names_column is None:
            raise Exception('At least one of sample_ids_column or sample_names_column must be set.')
        elif sample_ids_column is not None:
            return DataFrameSamplesQuery.create_with_sample_ids_column(sample_ids_column=sample_ids_column,
                                                                       data_frame=data_frame,
                                                                       wrapped_query=self,
                                                                       connection=self._query_connection,
                                                                       default_isa_kind=default_isa_kind,
                                                                       default_isa_column=default_isa_column)
        else:
            return DataFrameSamplesQuery.create_with_sample_names_column(sample_names_column=sample_names_column,
                                                                         data_frame=data_frame,
                                                                         wrapped_query=self,
                                                                         connection=self._query_connection,
                                                                         default_isa_kind=default_isa_kind,
                                                                         default_isa_column=default_isa_column)

    def join_tree(self, tree: Tree, kind='mutation', **kwargs) -> SamplesQuery:
        return TreeSamplesQueryFactory.instance().join_tree(tree=tree,
                                                            kind=kind,
                                                            database_connection=self._query_connection,
                                                            wrapped_query=self,
                                                            **kwargs)

    def intersect(self, sample_set: SampleSet, query_message: str = None) -> SamplesQuery:
        intersected_set = self._intersect_sample_set(sample_set)

        if query_message is None:
            query_message = f'intersect(samples={len(sample_set)}'

        queries_collection = self._queries_collection.append(query_message)
        return self._create_from(sample_set=intersected_set, universe_set=self._universe_set,
                                 unknown_set=self._unknown_set,
                                 queries_collection=queries_collection)

    def _found_in_self_and(self, found_in_other: SampleSet) -> SampleSet:
        return self._sample_set.intersection(found_in_other)

    def _found_in_self_or(self, found_in_other: SampleSet) -> SampleSet:
        return self._sample_set.union(found_in_other)

    def _unknown_in_self_or(self, unknown_in_other: SampleSet, found_either: SampleSet):
        return self._unknown_set.union(unknown_in_other).minus(found_either)

    def _unknown_in_self_and(self, found_in_other: SampleSet, unknown_in_other: SampleSet) -> SampleSet:
        """
        Given the above SampleSets, this returns the set of those unknown in self and the other sets
        (representing another query). This is defined using the help of Kleene's three-valued logic truth tables
        <https://en.wikipedia.org/wiki/Three-valued_logic#Kleene_and_Priest_logics>. Specifically, for two queries A
        and B with sets of samples consisting of either found (True), Unknown, or absent (False), this will return the
        sets of samples in the unknown state (U) in the truth table for AND(A,B).
        :param found_in_other: The SampleSet found in the other query.
        :param unknown_in_other: The SampleSet of unknowns in the other query.
        :return: The SampleSet of unknowns in self AND the other query.
        """
        # Each corresponds to one of three possible combinations in the three-valued truth-table for A and B
        # For those samples that should be in the unknown state (U) for "A AND B".
        unknown_and_unknown = self._unknown_set.intersection(unknown_in_other)
        unknown_and_found = self._unknown_set.intersection(found_in_other)
        found_and_unknown = self._sample_set.intersection(unknown_in_other)
        return unknown_and_unknown.union(unknown_and_found).union(found_and_unknown)

    def _intersect_sample_set(self, other: SampleSet) -> SampleSet:
        return self.sample_set.intersection(other)

    def _union_sample_set(self, other: SampleSet) -> SampleSet:
        return self.sample_set.union(other)

    def _union_unknown_set(self, other: SampleSet) -> SampleSet:
        return self.unknown_set.union(other)

    def _intersect_unknown_set(self, other: SampleSet) -> SampleSet:
        return self.unknown_set.intersection(other)

    def _get_has_kinds(self) -> List[str]:
        return self.HAS_KINDS

    def query_expression(self) -> str:
        return self._queries_collection.query_expression()

    def build_tree(self, kind: str, **kwargs) -> SamplesQuery:
        return TreeSamplesQueryFactory.instance().build_tree(kind=kind,
                                                             database_connection=self._query_connection,
                                                             wrapped_query=self, **kwargs)

    def toframe(self, include_present: bool = True, include_unknown: bool = False,
                include_absent: bool = False) -> pd.DataFrame:
        sample_service = self._query_connection.sample_service
        queries_expression = self._queries_collection.query_expression()

        present_set = self.sample_set if include_present else SampleSet.create_empty()
        absent_set = self.absent_set if include_absent else SampleSet.create_empty()
        unknown_set = self.unknown_set if include_unknown else SampleSet.create_empty()

        return sample_service.create_dataframe_from_sample_set(present_set=present_set,
                                                               absent_set=absent_set,
                                                               unknown_set=unknown_set,
                                                               queries_expression=queries_expression)

    def summary(self) -> pd.DataFrame:
        present = len(self)
        total = len(self.universe_set)
        unknown = len(self.unknown_set)
        absent = total - present - unknown
        per_present = (present / total) * 100
        per_absent = (absent / total) * 100
        per_unknown = (unknown / total) * 100

        return pd.DataFrame([{
            'Query': self._queries_collection.query_expression(),
            'Present': present,
            'Absent': absent,
            'Unknown': unknown,
            'Total': total,
            '% Present': per_present,
            '% Absent': per_absent,
            '% Unknown': per_unknown,
        }])

    def summary_features(self, kind: str = 'mutations', selection: str = 'all', **kwargs) -> pd.DataFrame:
        if kind == 'mutations':
            return self._summary_features_mutations(kind=kind, selection=selection, **kwargs)
        else:
            raise Exception(f'Unsupported value kind=[{kind}]. Must be one of {self.SUMMARY_FEATURES_KINDS}.')

    def _summary_features_mutations(self, kind: str, selection: str = 'all',
                                    ncores: int = 1,
                                    batch_size: int = 500,
                                    mutation_type: str = 'all',
                                    ignore_annotations: bool = False) -> pd.DataFrame:
        if selection not in self.FEATURES_SELECTIONS:
            raise Exception(f'selection=[{selection}] is unknown. Must be one of {self.FEATURES_SELECTIONS}')

        vs = self._query_connection.variation_service
        features_all_df = vs.count_mutations_in_sample_ids_dataframe(sample_ids=self._sample_set,
                                                                     ncores=ncores,
                                                                     batch_size=batch_size,
                                                                     mutation_type=mutation_type
                                                                     )
        features_all_df['Total'] = len(self)
        features_all_df['Percent'] = 100 * (features_all_df['Count'] / features_all_df['Total'])

        if selection == 'all':
            features_results_df = features_all_df
        elif selection == 'unique':
            features_complement_df = self.complement().summary_features(kind=kind, selection='all',
                                                                        ncores=ncores, batch_size=batch_size,
                                                                        mutation_type=mutation_type)
            features_merged_df = features_all_df.merge(features_complement_df, left_index=True, right_index=True,
                                                       how='left', indicator=True, suffixes=('_x', '_y'))
            features_merged_df = features_merged_df[features_merged_df['_merge'] == 'left_only'].rename({
                'Sequence_x': 'Sequence',
                'Position_x': 'Position',
                'Deletion_x': 'Deletion',
                'Insertion_x': 'Insertion',
                'Count_x': 'Count',
                'Total_x': 'Total',
                'Percent_x': 'Percent',
            }, axis='columns')
            features_results_df = features_merged_df[['Sequence', 'Position', 'Deletion', 'Insertion',
                                                      'Count', 'Total', 'Percent']]
        else:
            raise Exception(f'selection=[{selection}] is unknown. Must be one of {self.FEATURES_SELECTIONS}')

        if not ignore_annotations:
            features_results_df = vs.append_mutation_annotations(features_results_df)

        return features_results_df

    def tofeaturesset(self, kind: str = 'mutations', selection: str = 'all',
                      ncores: int = 1) -> Set[str]:
        return set(self.summary_features(kind=kind, selection=selection, ncores=ncores,
                                         ignore_annotations=True).index)

    def and_(self, other):
        if isinstance(other, SamplesQuery):
            present_set = self._found_in_self_and(other.sample_set)
            unknown_set = self._unknown_in_self_and(other.sample_set, other.unknown_set)
            universe_set = self.universe_set.union(other.universe_set)
            queries_collection = self._queries_collection.append(str(other))
            return self._create_from(present_set, universe_set=universe_set,
                                     unknown_set=unknown_set,
                                     queries_collection=queries_collection)
        else:
            raise Exception(f'Cannot perform an "and" on object {other}')

    def or_(self, other: SamplesQuery) -> SamplesQuery:
        if isinstance(other, SamplesQuery):
            present_set = self._found_in_self_or(other.sample_set)
            unknown_set = self._unknown_in_self_or(unknown_in_other=other.unknown_set, found_either=present_set)
            universe_set = self.universe_set.union(other.universe_set)

            queries_collection = self._queries_collection.append(f'OR({str(other)}')
            return self._create_from(present_set, universe_set=universe_set,
                                     unknown_set=unknown_set,
                                     queries_collection=queries_collection)
        else:
            raise Exception(f'Cannot perform an "or" on object {other}')

    def is_empty(self, include_unknown=False) -> bool:
        if include_unknown:
            return self.sample_set.is_empty() and self.unknown_set.is_empty()
        else:
            return self.sample_set.is_empty()

    def complement(self):
        complement_present_set = self.absent_set
        query_collection = self._queries_collection.append('complement')
        return self._create_from(sample_set=complement_present_set, universe_set=self._universe_set,
                                 unknown_set=self._unknown_set,
                                 queries_collection=query_collection)

    @property
    def tree(self):
        raise Exception(f'No tree exists for {self.__class__}.'
                        f' Perhaps you should try to run build_tree() first to build a tree.')

    def has_tree(self) -> bool:
        return False

    def __len__(self):
        return len(self.sample_set)

    def __str__(self) -> str:
        universe_length = len(self.universe_set)
        unknown_length = len(self.unknown_set)
        selected_length = len(self)
        percent_selected = (selected_length / universe_length) * 100
        percent_unknown = (unknown_length / universe_length) * 100
        return f'<{self.__class__.__name__}[selected={percent_selected:0.0f}% ' \
               f'({selected_length}/{universe_length}) samples, unknown={percent_unknown:0.0f}% ' \
               f'({unknown_length}/{universe_length}) samples]>'

    def tolist(self, names=True):
        if names:
            sample_service = self._query_connection.sample_service
            return [s.name for s in sample_service.find_samples_by_ids(self._sample_set)]
        else:
            return list(self.sample_set)

    def __repr__(self) -> str:
        return str(self)

    def hasa(self, property: Union[QueryFeature, str, pd.Series], kind='mutation') -> SamplesQuery:
        if isinstance(property, QueryFeature):
            query_feature = property
        elif isinstance(property, pd.Series):
            raise Exception(f'The query type {self.__class__.__name__} cannot support querying with respect to a '
                            f'dataframe. Perhaps you could try attaching a dataframe with join() first before querying.')
        elif kind is None:
            raise Exception(f'property=[{property}] is not of type QueryFeature so must set "kind" parameter')
        elif kind == 'mutation' or kind == 'mutations':
            if property.startswith('hgvs:'):
                query_feature = QueryFeatureHGVS.create_from_id(property)
            else:
                query_feature = QueryFeatureMutationSPDI(property)
        elif kind == 'mlst':
            query_feature = QueryFeatureMLST(property)
        else:
            raise Exception(f'kind={kind} is not recognized for {self}. Must be one of {self._get_has_kinds()}')

        found_hasa_set_dict = self._query_connection.sample_service.find_sample_sets_by_features([query_feature])
        unknown_hasa_set_dict = self._query_connection.sample_service.find_unknown_sample_sets_by_features([query_feature])

        if query_feature.id in found_hasa_set_dict:
            found_hasa_set = found_hasa_set_dict[query_feature.id]
            found_in_query = self._found_in_self_and(found_hasa_set)
        else:
            found_in_query = SampleSet.create_empty()

        if query_feature.id in unknown_hasa_set_dict:
            unknown_hasa_set = unknown_hasa_set_dict[query_feature.id]
        else:
            unknown_hasa_set = SampleSet.create_empty()

        unknown_in_query = self._unknown_in_self_and(found_in_other=found_in_query, unknown_in_other=unknown_hasa_set)

        # Universe remains the same in this case.
        universe_in_query = self._universe_set

        queries_collection = self._queries_collection.append(query_feature)
        return self._create_from(found_in_query, universe_set=universe_in_query,
                                 unknown_set=unknown_in_query,
                                 queries_collection=queries_collection)

    def _prepare_isin_query_message(self, query_message_prefix: str,
                                    query_msg_infix: str,
                                    additional_messages: str) -> str:
        return f"{query_message_prefix}({query_msg_infix}{additional_messages})"

    def _isin_samples(self, data: Union[str, List[str], pd.Series, SamplesQuery, SampleSet],
                      query_message_prefix: str) -> SamplesQuery:
        present_set, unknown_set, query_infix = self._get_present_unknown_sets_query_infix_from_data(data)
        query_message = self._prepare_isin_query_message(query_message_prefix=query_message_prefix,
                                                         query_msg_infix=query_infix,
                                                         additional_messages='')
        queries_collection = self._queries_collection.append(query_message)
        return self._create_from(sample_set=present_set, universe_set=self._universe_set,
                                 unknown_set=unknown_set,
                                 queries_collection=queries_collection)

    def _within_kmer_jaccard(self, data: Union[str, List[str], pd.Series, SamplesQuery, SampleSet],
                             distance: float, kmer_size: int = 31) -> SamplesQuery:
        additional_messages = f', dist={distance}, k={kmer_size}'
        sample_names, query_infix = self._get_sample_names_query_infix_from_data(data)
        query_message = self._prepare_isin_query_message(query_message_prefix='isin_kmer_jaccard',
                                                         query_msg_infix=query_infix,
                                                         additional_messages=additional_messages)

        kmer_service: KmerService = self._query_connection.kmer_service

        if len(sample_names) == 0:
            sample_set_matches = SampleSet.create_empty()
        else:
            sample_set_matches = kmer_service.find_matches_within(sample_names=list(sample_names),
                                                                  distance_threshold=distance,
                                                                  kmer_size=kmer_size,
                                                                  samples_universe=self._sample_set)
        queries_collection = self._queries_collection.append(query_message)
        return self._create_from(sample_set=sample_set_matches, universe_set=self._universe_set,
                                 unknown_set=self.unknown_set,
                                 queries_collection=queries_collection)

    def _within_distance(self, data: Union[str, List[str], pd.Series, SamplesQuery, SampleSet], distance: float,
                         units: str = 'kmer_jaccard', **kwargs) -> SamplesQuery:
        if units == 'kmer_jaccard':
            return self._within_kmer_jaccard(data=data, distance=distance, **kwargs)
        else:
            raise Exception(f'units=[{units}] is not supported. Must be one of {self._distance_units()}. '
                            f'For additional distance queries you perhaps need to build or attach a tree to '
                            f'the query.')

    def _get_present_unknown_sets_query_infix_from_data(self, data: Union[str, List[str], pd.Series, SamplesQuery, SampleSet]
                                              ) -> Tuple[SampleSet, SampleSet, str]:
        if isinstance(data, str) or isinstance(data, list):
            logger.debug(f'data=[{data}] contains sample names')
            if isinstance(data, str):
                query_msg = f"'{data}'"
                data = [data]
            else:
                query_msg = f'{data}'
            samples = self._query_connection.sample_service.get_existing_samples_by_names(data)
            sample_ids = {s.id for s in samples}
            present_set = SampleSet(sample_ids)
        elif isinstance(data, SamplesQuery) or isinstance(data, SampleSet):
            logger.debug(f'data=[{data}] contains sample ids')
            if isinstance(data, SamplesQuery):
                present_set = data.sample_set
                query_msg = f'{data}'
            else:
                present_set = data
                query_msg = f'set({len(data)} samples)'
        else:
            raise Exception(f'Unknown type for data=[{data}. Got type [{type(data)}]. '
                            'Must a string or list of strings (representing sample names) '
                            'or a set of sample IDs (as a SamplesQuery or SampleSet).')
        unknown_set = self.unknown_set.intersection(present_set)
        present_set = self.sample_set.intersection(present_set)
        return present_set, unknown_set, query_msg

    def _get_sample_names_query_infix_from_data(self, data: Union[str, List[str], pd.Series, SamplesQuery, SampleSet]
                                                ) -> Tuple[Set[str], str]:
        if isinstance(data, str) or isinstance(data, list):
            logger.debug(f'data={data} contains sample names')
            if isinstance(data, str):
                query_msg = f"'{data}'"
                data = [data]
            else:
                query_msg = f'{data}'
            sample_names = data
        elif isinstance(data, SamplesQuery) or isinstance(data, SampleSet):
            logger.debug(f'data=[{data}] contains sample ids')
            if isinstance(data, SamplesQuery):
                sample_set = data.sample_set
                query_msg = f'{data}'
            else:
                sample_set = data
                query_msg = f'set({len(data)} samples)'

            samples = self._query_connection.sample_service.find_samples_by_ids(sample_set)
            sample_names = {s.name for s in samples}
        else:
            raise Exception(f'Unknown type for data=[{data}. Got type [{type(data)}]. '
                            'Must a string or list of strings (representing sample names) '
                            'or a set of sample IDs (as a SamplesQuery or SampleSet).')
        return sample_names, query_msg

    def isin(self, data: Union[str, List[str], pd.Series, SamplesQuery, SampleSet], kind: str = 'samples',
             **kwargs) -> SamplesQuery:
        if kind == 'sample' or kind == 'samples':
            return self._isin_samples(data=data, query_message_prefix='isin_samples')
        elif kind == 'distance' or kind == 'distances':
            return self._within_distance(data=data, **kwargs)
        else:
            raise Exception(f'kind=[{kind}] is not supported. Must be one of {self.ISIN_TYPES}')

    def _isin_kinds(self) -> List[str]:
        return self.ISIN_TYPES

    def _isa_kinds(self) -> List[str]:
        return self.ISA_TYPES

    def _can_handle_isa_kind(self, kind: str) -> bool:
        return kind in self.ISA_TYPES

    def _can_handle_isin_kind(self, kind: str) -> bool:
        return kind in self.ISIN_TYPES

    def _distance_units(self) -> List[str]:
        return self.DISTANCES_UNITS

    def _can_handle_distance_units(self, units: str) -> bool:
        return units in self.DISTANCES_UNITS

    def isa(self, data: Union[str, List[str]], kind: str = 'sample', **kwargs) -> SamplesQuery:
        if kind == 'sample' or kind == 'samples':
            return self._isin_samples(data=data, query_message_prefix='isa_sample')
        else:
            raise Exception(f'kind=[{kind}] is not supported. Must be one of {self._isa_kinds()}')

    def to_distances(self, kind: str = 'kmer', **kwargs) -> Tuple[np.ndarray, List[str]]:
        if kind == 'kmer':
            return self._to_distances_kmer(**kwargs)
        else:
            raise Exception(f'kind=[{kind}] is not supported. Must be one of {self.DISTANCES_UNITS}')

    def _to_distances_kmer(self, kmer_size: int = 31, ncores: int = 1) -> Tuple[np.ndarray, List[str]]:
        return self._query_connection.kmer_service.get_distance_matrix(sample_ids=self._sample_set,
                                                                       kmer_size=kmer_size,
                                                                       ncores=ncores)

    def _create_from(self, sample_set: SampleSet, universe_set: SampleSet, unknown_set: SampleSet,
                     queries_collection: QueriesCollection) -> SamplesQuery:
        return SamplesQueryIndex(connection=self._query_connection,
                                 universe_set=universe_set,
                                 sample_set=sample_set,
                                 unknown_set=unknown_set,
                                 queries_collection=queries_collection)
