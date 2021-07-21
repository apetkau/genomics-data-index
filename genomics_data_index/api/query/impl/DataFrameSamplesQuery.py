from __future__ import annotations

from typing import Union, List

import pandas as pd
from ete3 import Tree

from genomics_data_index.api.query.SamplesQuery import SamplesQuery
from genomics_data_index.api.query.impl.TreeSamplesQueryFactory import TreeSamplesQueryFactory
from genomics_data_index.api.query.impl.WrappedSamplesQuery import WrappedSamplesQuery
from genomics_data_index.configuration.connector import DataIndexConnection
from genomics_data_index.storage.SampleSet import SampleSet


class DataFrameSamplesQuery(WrappedSamplesQuery):
    """
    Defines a SamplesQuery which is joined to a DataFrame. This allows for queries that make use of data in the DataFrame.
    """

    ISIN_KINDS = ['dataframe']
    ISA_KINDS = ['dataframe']

    def __init__(self, connection: DataIndexConnection, wrapped_query: SamplesQuery,
                 data_frame: pd.DataFrame,
                 sample_ids_col: str,
                 default_isa_kind: str,
                 default_isa_column: str):
        """
        Builds a new DataFrameSamplesQuery from the given information. In most normal operations DataFrameSamplesQuery objects
        are not created directly but are instead created from an :py:class:`genomics_data_index.api.GenomicsDataIndex`
        object or from the join operation applied to a SamplesQuery.

        :param connection: A connection to a database containing samples.
        :param wrapped_query: The SamplesQuery to wrap around/decorate.
        :param data_frame: The DataFrame to join to this query. This dataframe must have a column which matches rows to
                           the Sample IDs stored in database defined by the connection object.
        :param sample_ids_col: The name of the column in the DataFrame defining the Sample ID for a particular row.
        :param default_isa_kind: If A is a DataFrameSamplesQuery defines the default `kind` parameter to be used
                                 when running A.isa(data, kind='...').
        :param default_isa_column: If A is a DataFrameSamplesQuery defines the default 'isa_column' parameter
                                   to be used when running A.isa(data, kind='...', isa_column='...').
        :return: A new DataFrameSamplesQuery object.
        """
        super().__init__(connection=connection, wrapped_query=wrapped_query)
        self._sample_ids_col = sample_ids_col
        self._data_frame = data_frame
        self._default_isa_kind = default_isa_kind
        self._default_isa_column = default_isa_column

    def _isin_internal(self, data: Union[str, List[str], pd.Series, SamplesQuery, SampleSet], kind: str,
                       **kwargs) -> SamplesQuery:
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

    def isa(self, data: Union[str, List[str], pd.Series, SamplesQuery, SampleSet], kind: str = None,
            **kwargs) -> SamplesQuery:
        if kind is None:
            if self._default_isa_kind is None:
                kind = 'sample'
            else:
                kind = self._default_isa_kind

        return super().isa(data=data, kind=kind, **kwargs)

    def _isa_internal(self, data: Union[str, List[str], pd.Series, SamplesQuery, SampleSet], kind: str,
                      isa_column: str = None,
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

    def toframe(self, include_present: bool = True, include_unknown: bool = False,
                include_absent: bool = False) -> pd.DataFrame:
        samples_dataframe = super().toframe(include_present=include_present,
                                            include_unknown=include_unknown,
                                            include_absent=include_absent)
        samples_df_cols = list(samples_dataframe.columns)
        merged_df = self._data_frame.merge(samples_dataframe, how='inner', left_on=self._sample_ids_col,
                                           right_on='Sample ID')
        new_col_order = samples_df_cols + [col for col in list(merged_df.columns) if col not in samples_df_cols]
        return merged_df[new_col_order]

    def reset_universe(self, include_unknown: bool = True) -> SamplesQuery:
        new_wrapped_query = self._wrapped_query.reset_universe(include_unknown=include_unknown)
        return self._wrap_create(new_wrapped_query)

    def _wrap_create(self, wrapped_query: SamplesQuery) -> WrappedSamplesQuery:
        return DataFrameSamplesQuery(connection=self._query_connection,
                                     wrapped_query=wrapped_query,
                                     data_frame=self._data_frame,
                                     sample_ids_col=self._sample_ids_col,
                                     default_isa_kind=self._default_isa_kind,
                                     default_isa_column=self._default_isa_column)

    def build_tree(self, kind: str, **kwargs) -> SamplesQuery:
        return TreeSamplesQueryFactory.instance().build_tree(kind=kind,
                                                             database_connection=self._query_connection,
                                                             wrapped_query=self, **kwargs)

    def join_tree(self, tree: Tree, kind='mutation', **kwargs) -> SamplesQuery:
        return TreeSamplesQueryFactory.instance().join_tree(tree=tree,
                                                            kind=kind,
                                                            database_connection=self._query_connection,
                                                            wrapped_query=self,
                                                            **kwargs)

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
        """
        Builds a new DataFrameSamplesQuery from the given information matching by a column defining the Sample IDs.
        In most normal operations DataFrameSamplesQuery objects are not created directly but are instead created from
        an :py:class:`genomics_data_index.api.GenomicsDataIndex` object or from the join operation applied to a SamplesQuery.

        :param sample_ids_col: The name of the column in the DataFrame defining the Sample ID for a particular row.
        :param data_frame: The DataFrame to join to this query. This dataframe must have a column which matches rows to
                   the Sample IDs stored in database defined by the connection object.
        :param wrapped_query: The SamplesQuery to wrap around/decorate.
        :param connection: A connection to a database containing samples.
        :param query_message: The query_message to use when creating this new DataFrameSamplesQuery object.
        :param default_isa_kind: If A is a DataFrameSamplesQuery defines the default `kind` parameter to be used
                                 when running A.isa(data, kind='...').
        :param default_isa_column: If A is a DataFrameSamplesQuery defines the default 'isa_column' parameter
                                   to be used when running A.isa(data, kind='...', isa_column='...').
        :return: A new DataFrameSamplesQuery object.
        """
        sample_ids = data_frame[sample_ids_column].tolist()
        df_sample_set = SampleSet(sample_ids=sample_ids)
        universe_set = wrapped_query.universe_set.intersection(df_sample_set)

        if query_message is None:
            query_message = f'dataframe(ids_col=[{sample_ids_column}])'

        wrapped_query_intersect = wrapped_query.set_universe(universe_set, query_message=query_message)

        return DataFrameSamplesQuery(connection=connection,
                                     wrapped_query=wrapped_query_intersect,
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
        """
        Builds a new DataFrameSamplesQuery from the given information matching by a column defining the Sample names.
        In most normal operations DataFrameSamplesQuery objects are not created directly but are instead created from
        an :py:class:`genomics_data_index.api.GenomicsDataIndex` object or from the join operation applied to a SamplesQuery.

        :param sample_names_column: The name of the column in the DataFrame defining the Sample name for a particular row.
                                    This will internally be matched up to Sample IDs in the database defined by the connection
                                    object and used to create a new column in the DataFrame with the Sample ID.
        :param data_frame: The DataFrame to join to this query. This dataframe must have a column which matches rows to
                   the Sample names stored in database defined by the connection object.
        :param wrapped_query: The SamplesQuery to wrap around/decorate.
        :param connection: A connection to a database containing samples.
        :param default_isa_kind: If A is a DataFrameSamplesQuery defines the default `kind` parameter to be used
                                 when running A.isa(data, kind='...').
        :param default_isa_column: If A is a DataFrameSamplesQuery defines the default 'isa_column' parameter
                                   to be used when running A.isa(data, kind='...', isa_column='...').
        :return: A new DataFrameSamplesQuery object.
        """
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
