from __future__ import annotations

import pandas as pd

from storage.api.SamplesQuery import SamplesQuery
from storage.api.impl.QueriesCollection import QueriesCollection
from storage.api.impl.SamplesQueryIndex import SamplesQueryIndex
from storage.connector.DataIndexConnection import DataIndexConnection
from storage.variant.SampleSet import SampleSet


class DataFrameSamplesQuery(SamplesQueryIndex):

    def __init__(self, data_frame: pd.DataFrame,
                 sample_ids_col: str,
                 connection: DataIndexConnection,
                 queries_collection: QueriesCollection = QueriesCollection.create_empty(),
                 sample_set: SampleSet = SamplesQuery.ALL_SAMPLES):
        super().__init__(connection=connection, sample_set=sample_set,
                         queries_collection=queries_collection)
        self._sample_ids_col = sample_ids_col
        self._data_frame = data_frame

    def toframe(self) -> pd.DataFrame:
        samples_dataframe = super().toframe()
        merged_df = self._data_frame.merge(samples_dataframe, how='inner', left_on=self._sample_ids_col,
                                           right_on='Sample ID')
        return merged_df

    def _create_from(self, connection: DataIndexConnection, sample_set: SampleSet,
                     queries_collection: QueriesCollection) -> SamplesQuery:
        return DataFrameSamplesQuery(data_frame=self._data_frame,
                                     sample_ids_col=self._sample_ids_col,
                                     connection=connection,
                                     sample_set=sample_set,
                                     queries_collection=queries_collection)

    @classmethod
    def create_with_sample_ids_column(self, sample_ids_column: str, data_frame: pd.DataFrame,
                                      connection: DataIndexConnection) -> DataFrameSamplesQuery:
        sample_ids = data_frame[sample_ids_column].tolist()
        sample_set = SampleSet(sample_ids=sample_ids)

        return DataFrameSamplesQuery(data_frame=data_frame,
                                     sample_ids_col=sample_ids_column,
                                     connection=connection,
                                     sample_set=sample_set)

    @classmethod
    def create_with_sample_names_column(self, sample_names_column: str, data_frame: pd.DataFrame,
                                        connection: DataIndexConnection) -> DataFrameSamplesQuery:
        sample_names = set(data_frame[sample_names_column].tolist())
        sample_ids_column = 'Sample ID'

        sample_name_ids = connection.sample_service.find_sample_name_ids(sample_names)
        sample_set = SampleSet(sample_ids=sample_name_ids.values())

        # Only attempt once to rename sample IDs column if it already exists
        if sample_ids_column in data_frame:
            sample_ids_column = sample_ids_column + '_database'
            if sample_ids_column in data_frame:
                raise Exception(f'Column to be used for sample_ids [{sample_ids_column}] already in data frame')

        sample_name_series = pd.Series(sample_name_ids, name=sample_ids_column)
        data_frame = data_frame.merge(sample_name_series, left_on=sample_names_column, right_index=True)

        return DataFrameSamplesQuery(data_frame=data_frame,
                                     sample_ids_col=sample_ids_column,
                                     connection=connection,
                                     sample_set=sample_set)
