from __future__ import annotations

from typing import Union, List

import pandas as pd

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.impl.TreeSamplesQueryFactory import TreeSamplesQueryFactory
from genomics_data_index.api.query.impl.WrappedSamplesQuery import WrappedSamplesQuery
from genomics_data_index.configuration.connector import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet


class DataFrameSamplesQuery(WrappedSamplesQuery):
    ISIN_KINDS = ['dataframe']
    ISA_KINDS = ['dataframe']

    def __init__(self, connection: DataIndexConnection, wrapped_query: SamplesQuery,
                 universe_set: SampleSet,
                 data_frame: pd.DataFrame,
                 sample_ids_col: str,
                 default_isa_kind: str,
                 default_isa_column: str):
        super().__init__(connection=connection, wrapped_query=wrapped_query, universe_set=universe_set)
        self._sample_ids_col = sample_ids_col
        self._data_frame = data_frame
        self._default_isa_kind = default_isa_kind
        self._default_isa_column = default_isa_column

    def _isin_internal(self, data: Union[str, List[str], pd.Series, SamplesQuery, SampleSet], kind: str, **kwargs) -> SamplesQuery:
        if kind == 'dataframe':
            if isinstance(data, pd.Series) and data.dtype == bool:
                if data.index.equals(self._data_frame.index):
                    return self._handle_select_by_series(data, query_message='isin(subset from series)')
                else:
                    raise Exception(f'Passed data=[series with length {len(data)}] '
                                    f'does not have same index as internal data frame (length={len(self._data_frame)})')
            else:
                raise Exception(f'data=[{data}] is wrong type for kind=[{kind}]. '
                                f'Must be a boolean pandas.Series')
        else:
            raise Exception(f'kind=[{kind}] is not supported. Must be one of {self._isin_kinds()}')

    def _isin_kinds(self) -> List[str]:
        return list(set(super()._isin_kinds() + self.ISIN_KINDS))

    def _can_handle_isin_kind(self, kind: str) -> bool:
        return kind in self.ISIN_KINDS

    def _can_handle_isa_kind(self, kind: str) -> bool:
        return kind in self.ISA_KINDS

    def isa(self, data: Union[str, List[str], pd.Series, SamplesQuery, SampleSet], kind: str = None, **kwargs) -> SamplesQuery:
        if kind is None:
            if self._default_isa_kind is None:
                kind = 'sample'
            else:
                kind = self._default_isa_kind

        return super().isa(data=data, kind=kind, **kwargs)

    def _isa_internal(self, data: Union[str, List[str], pd.Series, SamplesQuery, SampleSet], kind: str, isa_column: str = None,
                      regex: bool = False) -> SamplesQuery:
        if kind == 'dataframe':
            if isa_column is None and self._default_isa_column is not None:
                isa_column = self._default_isa_column
            elif isa_column is None:
                raise Exception(f'No defined isa_column, cannot execute isa for kind={kind}')

            if regex:
                no_na_select = (~self._data_frame[isa_column].isna()) & (
                    self._data_frame[isa_column].str.contains(data))
                return self._handle_select_by_series(no_na_select,
                                                     query_message=f"isa('{isa_column}' contains '{data}')")
            else:
                return self._handle_select_by_series(self._data_frame[isa_column] == data,
                                                     query_message=f"isa('{isa_column}' is '{data}')")
        else:
            raise Exception(f'Invalid kind={kind}. Must be one of {self._isa_kinds()}')

    def _isa_kinds(self):
        return super()._isa_kinds() + self.ISA_KINDS

    def _handle_select_by_series(self, series_selection: pd.Series, query_message: str) -> SamplesQuery:
        subset_df = self._data_frame[series_selection]
        subset_sample_set = SampleSet(subset_df[self._sample_ids_col].tolist())
        subset_query = self._wrapped_query.intersect(sample_set=subset_sample_set, query_message=query_message)
        return self._wrap_create(subset_query)

    def toframe(self, exclude_absent: bool = True) -> pd.DataFrame:
        samples_dataframe = super().toframe()
        samples_df_cols = list(samples_dataframe.columns)
        merged_df = self._data_frame.merge(samples_dataframe, how='inner', left_on=self._sample_ids_col,
                                           right_on='Sample ID')
        new_col_order = samples_df_cols + [col for col in list(merged_df.columns) if col not in samples_df_cols]
        return merged_df[new_col_order]

    def reset_universe(self) -> SamplesQuery:
        new_wrapped_query = self._wrapped_query.reset_universe()
        return self._wrap_create(new_wrapped_query, universe_set=new_wrapped_query.universe_set)

    def _wrap_create(self, wrapped_query: SamplesQuery, universe_set: SampleSet = None) -> WrappedSamplesQuery:
        if universe_set is None:
            universe_set = self.universe_set

        return DataFrameSamplesQuery(connection=self._query_connection,
                                     wrapped_query=wrapped_query,
                                     data_frame=self._data_frame,
                                     sample_ids_col=self._sample_ids_col,
                                     universe_set=universe_set,
                                     default_isa_kind=self._default_isa_kind,
                                     default_isa_column=self._default_isa_column)

    def build_tree(self, kind: str, **kwargs) -> SamplesQuery:
        return TreeSamplesQueryFactory.instance().build_tree(kind=kind,
                                                             database_connection=self._query_connection,
                                                             wrapped_query=self, **kwargs)

    def join(self, data_frame: pd.DataFrame, sample_ids_column: str = None,
             sample_names_column: str = None, default_isa_kind: str = 'names',
             default_isa_column: str = None) -> SamplesQuery:
        raise Exception(f'Cannot join a new dataframe onto an existing data frame query: {self}')

    @classmethod
    def create_with_sample_ids_column(self, sample_ids_column: str, data_frame: pd.DataFrame,
                                      wrapped_query: SamplesQuery, connection: DataIndexConnection,
                                      query_message: str = None,
                                      default_isa_kind: str = None,
                                      default_isa_column: str = None) -> DataFrameSamplesQuery:
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
                                     sample_ids_col=sample_ids_column,
                                     default_isa_column=default_isa_column,
                                     default_isa_kind=default_isa_kind)

    @classmethod
    def create_with_sample_names_column(self, sample_names_column: str, data_frame: pd.DataFrame,
                                        wrapped_query: SamplesQuery,
                                        connection: DataIndexConnection,
                                        default_isa_kind: str = None,
                                        default_isa_column: str = None) -> DataFrameSamplesQuery:
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
                                                                   query_message=query_message,
                                                                   default_isa_column=default_isa_column,
                                                                   default_isa_kind=default_isa_kind)
