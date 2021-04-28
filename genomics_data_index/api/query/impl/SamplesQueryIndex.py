from __future__ import annotations

from typing import Union, List, Set

import pandas as pd

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.impl.DataFrameSamplesQuery import DataFrameSamplesQuery
from genomics_data_index.api.query.impl.QueriesCollection import QueriesCollection
from genomics_data_index.api.query.impl.TreeSamplesQueryFactory import TreeSamplesQueryFactory
from genomics_data_index.configuration.connector import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.QueryFeature import QueryFeature
from genomics_data_index.storage.model.QueryFeatureMLST import QueryFeatureMLST
from genomics_data_index.storage.model.QueryFeatureMutation import QueryFeatureMutation


class SamplesQueryIndex(SamplesQuery):
    HAS_KINDS = ['mutation', 'mutations', 'mlst']
    SUMMARY_FEATURES_KINDS = ['mutations']
    FEATURES_SELECTIONS = ['all', 'unique']
    ISIN_TYPES = ['names']
    ISA_TYPES = ['names']

    def __init__(self, connection: DataIndexConnection,
                 universe_set: SampleSet,
                 sample_set: SampleSet,
                 queries_collection: QueriesCollection = QueriesCollection.create_empty()):
        super().__init__()
        self._query_connection = connection
        self._universe_set = universe_set
        self._sample_set = sample_set
        self._queries_collection = queries_collection

    @property
    def universe_set(self) -> SampleSet:
        return self._universe_set

    @property
    def sample_set(self) -> SampleSet:
        return self._sample_set

    def reset_universe(self) -> SamplesQuery:
        return self._create_from(sample_set=self._sample_set, universe_set=self._sample_set,
                                 queries_collection=self._queries_collection)

    def join(self, data_frame: pd.DataFrame, sample_ids_column: str = None,
             sample_names_column: str = None, default_isa_kind: str = 'names',
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

    def intersect(self, sample_set: SampleSet, query_message: str = None) -> SamplesQuery:
        intersected_set = self._intersect_sample_set(sample_set)

        if query_message is None:
            query_message = f'intersect(samples={len(sample_set)}'

        queries_collection = self._queries_collection.append(query_message)
        return self._create_from(sample_set=intersected_set, universe_set=self._universe_set,
                                 queries_collection=queries_collection)

    def _intersect_sample_set(self, other: SampleSet) -> SampleSet:
        return self.sample_set.intersection(other)

    def _union_sample_set(self, other: SampleSet) -> SampleSet:
        return self.sample_set.union(other)

    def _get_has_kinds(self) -> List[str]:
        return self.HAS_KINDS

    def query_expression(self) -> str:
        return self._queries_collection.query_expression()

    def build_tree(self, kind: str, scope: str, **kwargs) -> SamplesQuery:
        return TreeSamplesQueryFactory.instance().build_tree(kind=kind, scope=scope,
                                                             database_connection=self._query_connection,
                                                             wrapped_query=self, **kwargs)

    def toframe(self, exclude_absent: bool = True) -> pd.DataFrame:
        sample_service = self._query_connection.sample_service
        queries_expression = self._queries_collection.query_expression()
        return sample_service.create_dataframe_from_sample_set(self.sample_set,
                                                               universe_set=self.universe_set,
                                                               exclude_absent=exclude_absent,
                                                               queries_expression=queries_expression)

    def summary(self) -> pd.DataFrame:
        present = len(self)
        total = len(self.universe_set)
        unknown = pd.NA
        absent = total - present
        per_present = (present / total) * 100
        per_absent = (absent / total) * 100
        per_unknown = pd.NA

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

    def summary_features(self, kind: str = 'mutations', **kwargs) -> pd.DataFrame:
        if kind == 'mutations':
            return self._summary_features_mutations(kind=kind, **kwargs)
        else:
            raise Exception(f'Unsupported value kind=[{kind}]. Must be one of {self.SUMMARY_FEATURES_KINDS}.')

    def _summary_features_mutations(self, kind: str, ncores: int = 1,
                                    batch_size: int = 500,
                                    mutation_type: str = 'all'):
        vs = self._query_connection.variation_service
        return vs.count_mutations_in_sample_ids_dataframe(sample_ids=self._sample_set,
                                                          ncores=ncores,
                                                          batch_size=batch_size,
                                                          mutation_type=mutation_type
                                                          )

    def tofeaturesset(self, kind: str = 'mutations', selection: str = 'all',
                      ncores: int = 1) -> Set[str]:
        if selection == 'all':
            return set(self.summary_features(kind=kind, ncores=ncores).index)
        elif selection == 'unique':
            features_set = set(self.summary_features(kind=kind, ncores=ncores).index)
            complement_features_set = set(self.complement().summary_features(kind=kind, ncores=ncores).index)
            return features_set - complement_features_set
        else:
            raise Exception(f'Unsupported selection=[{selection}]. Must be one of {self.FEATURES_SELECTIONS}.')

    def and_(self, other):
        if isinstance(other, SamplesQuery):
            intersect_set = self._intersect_sample_set(other.sample_set)
            queries_collection = self._queries_collection.append(str(other))
            return self._create_from(intersect_set, universe_set=self._universe_set,
                                     queries_collection=queries_collection)
        else:
            raise Exception(f'Cannot perform an "and" on object {other}')

    def or_(self, other: SamplesQuery) -> SamplesQuery:
        if isinstance(other, SamplesQuery):
            union_set = self._union_sample_set(other.sample_set)
            queries_collection = self._queries_collection.append(f'OR({str(other)}')
            return self._create_from(union_set, universe_set=self._universe_set,
                                     queries_collection=queries_collection)
        else:
            raise Exception(f'Cannot perform an "or" on object {other}')

    def is_empty(self) -> bool:
        return self.sample_set.is_empty()

    def complement(self):
        complement_set = self.universe_set.minus(self.sample_set)
        query_collection = self._queries_collection.append('complement')
        return self._create_from(sample_set=complement_set, universe_set=self._universe_set,
                                 queries_collection=query_collection)

    @property
    def tree(self):
        raise Exception(f'No tree exists for {self.__class__}.'
                        f' Perhaps you should try to run build_tree() first to build a tree.')

    def __and__(self, other):
        return self.and_(other)

    def __or__(self, other):
        return self.or_(other)

    def __len__(self):
        return len(self.sample_set)

    def __str__(self) -> str:
        universe_length = len(self.universe_set)
        percent_selected = (len(self) / universe_length) * 100
        return f'<{self.__class__.__name__}[{percent_selected:0.0f}% ({len(self)}/{universe_length}) samples]>'

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
            query_feature = QueryFeatureMutation(property)
        elif kind == 'mlst':
            query_feature = QueryFeatureMLST(property)
        else:
            raise Exception(f'kind={kind} is not recognized for {self}. Must be one of {self._get_has_kinds()}')

        found_set_dict = self._query_connection.sample_service.find_sample_sets_by_features([query_feature])

        if query_feature.id in found_set_dict:
            found_set = found_set_dict[query_feature.id]
            intersect_found = self._intersect_sample_set(found_set)
        else:
            intersect_found = SampleSet.create_empty()

        queries_collection = self._queries_collection.append(query_feature)
        return self._create_from(intersect_found, universe_set=self._universe_set,
                                 queries_collection=queries_collection)

    def _isin_names(self, sample_names: Union[str, List[str]], query_message_prefix: str) -> SamplesQuery:
        if isinstance(sample_names, str):
            query_message = f"{query_message_prefix}('{sample_names}')"
            sample_names = [sample_names]
        elif isinstance(sample_names, list):
            query_message = f"{query_message_prefix}({sample_names})"
        else:
            raise Exception(f'Unrecognized input type sample_names={sample_names}. Can only be a str or list')

        samples = self._query_connection.sample_service.get_existing_samples_by_names(sample_names)
        sample_ids = {s.id for s in samples}
        sample_set = SampleSet(sample_ids)
        queries_collection = self._queries_collection.append(query_message)
        return self._create_from(sample_set=sample_set, universe_set=self._universe_set,
                                 queries_collection=queries_collection)

    def isin(self, data: Union[str, List[str]], kind: str = 'names', **kwargs) -> SamplesQuery:
        if kind == 'names':
            return self._isin_names(sample_names=data, query_message_prefix='isin')
        else:
            raise Exception(f'kind=[{kind}] is not supported. Must be one of {self.ISIN_TYPES}')

    def isa(self, data: Union[str, List[str]], kind: str = 'names', **kwargs) -> SamplesQuery:
        if kind == 'names':
            return self._isin_names(sample_names=data, query_message_prefix='isa')
        else:
            raise Exception(f'kind=[{kind}] is not supported. Must be one of {self.ISA_TYPES}')

    def _create_from(self, sample_set: SampleSet, universe_set: SampleSet,
                     queries_collection: QueriesCollection) -> SamplesQuery:
        return SamplesQueryIndex(connection=self._query_connection,
                                 universe_set=universe_set,
                                 sample_set=sample_set,
                                 queries_collection=queries_collection)
