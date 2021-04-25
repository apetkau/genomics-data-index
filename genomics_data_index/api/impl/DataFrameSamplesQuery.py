from __future__ import annotations

from typing import Union, List

import pandas as pd

from genomics_data_index.api.SamplesQuery import SamplesQuery
from genomics_data_index.api.impl.TreeSamplesQuery import TreeSamplesQuery
from genomics_data_index.api.impl.WrappedSamplesQuery import WrappedSamplesQuery
from genomics_data_index.configuration.connector import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.model.QueryFeature import QueryFeature


class DataFrameSamplesQuery(WrappedSamplesQuery):
    HAS_KINDS = ['dataframe']

    def __init__(self, connection: DataIndexConnection, wrapped_query: SamplesQuery,
                 universe_set: SampleSet,
                 data_frame: pd.DataFrame,
                 sample_ids_col: str):
        super().__init__(connection=connection, wrapped_query=wrapped_query, universe_set=universe_set)
        self._sample_ids_col = sample_ids_col
        self._data_frame = data_frame

    def has(self, property: Union[QueryFeature, str, pd.Series], kind=None) -> SamplesQuery:
        if kind == 'dataframe':
            if isinstance(property, pd.Series) and property.dtype == bool:
                if property.index.equals(self._data_frame.index):
                    return self._handle_select_by_series(property)
                else:
                    raise Exception(f'Passed property=[series with length {len(property)}] '
                                    f'does not have same index as internal data frame (length={len(self._data_frame)})')
            else:
                raise Exception(f'property=[{property}] is wrong type for kind=[{kind}]. '
                                f'Must be a boolean pandas.Series')
        else:
            return self._wrap_create(self._wrapped_query.has(property=property, kind=kind))

    def _get_has_kinds(self) -> List[str]:
        return self._wrapped_query._get_has_kinds() + self.HAS_KINDS

    def _handle_select_by_series(self, series_selection: pd.Series) -> SamplesQuery:
        subset_df = self._data_frame[series_selection]
        subset_sample_set = SampleSet(subset_df[self._sample_ids_col].tolist())
        query_message = 'has(subset from series)'
        subset_query = self._wrapped_query.intersect(sample_set=subset_sample_set, query_message=query_message)
        return self._wrap_create(subset_query)

    def toframe(self, exclude_absent: bool = True) -> pd.DataFrame:
        samples_dataframe = super().toframe()
        samples_df_cols = list(samples_dataframe.columns)
        merged_df = self._data_frame.merge(samples_dataframe, how='inner', left_on=self._sample_ids_col,
                                           right_on='Sample ID')
        new_col_order = samples_df_cols + [col for col in list(merged_df.columns) if col not in samples_df_cols]
        return merged_df[new_col_order]

    def _wrap_create(self, wrapped_query: SamplesQuery) -> WrappedSamplesQuery:
        return DataFrameSamplesQuery(connection=self._query_connection,
                                     wrapped_query=wrapped_query,
                                     data_frame=self._data_frame,
                                     sample_ids_col=self._sample_ids_col,
                                     universe_set=self.universe_set)

    def build_tree(self, kind: str, scope: str, **kwargs) -> SamplesQuery:
        return TreeSamplesQuery.create(kind=kind, scope=scope, database_connection=self._query_connection,
                                       wrapped_query=self, **kwargs)

    @classmethod
    def create_with_sample_ids_column(self, sample_ids_column: str, data_frame: pd.DataFrame,
                                      wrapped_query: SamplesQuery, connection: DataIndexConnection,
                                      query_message: str = None) -> DataFrameSamplesQuery:
        sample_ids = data_frame[sample_ids_column].tolist()
        df_sample_set = SampleSet(sample_ids=sample_ids)
        universe_set = wrapped_query.universe_set.intersection(df_sample_set)

        if query_message is None:
            query_message = f'dataframe(ids_col=[{sample_ids_column}])'

        wrapped_query_intersect = wrapped_query.intersect(sample_set=df_sample_set,
                                                          query_message=query_message)

        return DataFrameSamplesQuery(connection=connection,
                                     wrapped_query=wrapped_query_intersect,
                                     universe_set=universe_set,
                                     data_frame=data_frame,
                                     sample_ids_col=sample_ids_column)

    def join(self, data_frame: pd.DataFrame, sample_ids_column: str = None,
             sample_names_column: str = None) -> SamplesQuery:
        raise Exception(f'Cannot join a new dataframe onto an existing data frame query: {self}')

    @classmethod
    def create_with_sample_names_column(self, sample_names_column: str, data_frame: pd.DataFrame,
                                        wrapped_query: SamplesQuery,
                                        connection: DataIndexConnection) -> DataFrameSamplesQuery:
        sample_names = set(data_frame[sample_names_column].tolist())
        sample_ids_column = 'Sample ID'

        sample_name_ids = connection.sample_service.find_sample_name_ids(sample_names)

        # Only attempt once to rename sample IDs column if it already exists
        if sample_ids_column in data_frame:
            sample_ids_column = sample_ids_column + '_gdi'
            if sample_ids_column in data_frame:
                raise Exception(f'Column to be used for sample_ids [{sample_ids_column}] already in data frame')

        sample_ids_series = pd.Series(sample_name_ids, name=sample_ids_column)
        data_frame = data_frame.merge(sample_ids_series, left_on=sample_names_column, right_index=True)

        query_message = f'dataframe(names_col=[{sample_names_column}])'

        return DataFrameSamplesQuery.create_with_sample_ids_column(sample_ids_column=sample_ids_column,
                                                                   data_frame=data_frame,
                                                                   wrapped_query=wrapped_query,
                                                                   connection=connection,
                                                                   query_message=query_message)
