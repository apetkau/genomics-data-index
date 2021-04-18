from __future__ import annotations

import pandas as pd

from storage.api.SamplesQuery import SamplesQuery
from storage.api.impl.SamplesQueryIndex import SamplesQueryIndex
from storage.variant.SampleSet import SampleSet
from storage.connector.DataIndexConnection import DataIndexConnection


class DataFrameSamplesQuery(SamplesQueryIndex):

    def __init__(self, data_frame: pd.DataFrame,
                 sample_ids_col: str,
                 connection: DataIndexConnection,
                 sample_set: SampleSet = SamplesQuery.ALL_SAMPLES):
        super().__init__(connection=connection, sample_set=sample_set)
        self._sample_ids_col = sample_ids_col
        self._data_frame = data_frame

    def dataframe(self) -> pd.DataFrame:
        samples_dataframe = super().dataframe()
        merged_df = self._data_frame.merge(samples_dataframe, how='inner', left_on=self._sample_ids_col,
                               right_on='Sample ID')
        return merged_df

    def _create_from(self, connection: DataIndexConnection, sample_set: SampleSet) -> SamplesQuery:
        return DataFrameSamplesQuery(data_frame=self._data_frame,
                                     sample_ids_col=self._sample_ids_col,
                                     connection=connection,
                                     sample_set=sample_set)

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
        sample_names = data_frame[sample_names_column].tolist()
        sample_ids_col = 'Sample ID'

        sample_name_ids = {} # do lookup based on sample names
        raise Exception('Not implemented')
        sample_set = SampleSet(sample_ids=sample_name_ids.values())

        if sample_ids_col in data_frame:
            raise Exception(f'Column to be used for sample_ids [{sample_ids_col}] already in data frame')

        data_frame = data_frame.merge()

        return DataFrameSamplesQuery(data_frame=data_frame,
                                     sample_ids_col=sample_ids_col,
                                     connection=connection,
                                     sample_set=sample_set)
