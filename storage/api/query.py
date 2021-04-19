from __future__ import annotations

from pathlib import Path

import pandas as pd

from storage.api.SamplesQuery import SamplesQuery
from storage.api.impl.DataFrameSamplesQuery import DataFrameSamplesQuery
from storage.api.impl.SamplesQueryIndex import SamplesQueryIndex
from storage.connector.DataIndexConnection import DataIndexConnection


def connect(database_connection: str, database_dir: Path) -> DataIndexConnection:
    return DataIndexConnection.connect(database_connection=database_connection,
                                       database_dir=database_dir)


def query(connection: DataIndexConnection,
          data_frame: pd.DataFrame = None,
          sample_ids_column=None, sample_names_column=None) -> SamplesQuery:
    if data_frame is not None:
        if sample_ids_column is None and sample_names_column is None:
            raise Exception('If querying with a data_frame, then one of sample_names_column or sample_ids_column must '
                            'be set')
        elif sample_ids_column is not None:
            return DataFrameSamplesQuery.create_with_sample_ids_column(sample_ids_column,
                                                                       data_frame=data_frame,
                                                                       connection=connection)
        else:
            return DataFrameSamplesQuery.create_with_sample_names_column(sample_names_column,
                                                                         data_frame=data_frame,
                                                                         connection=connection)
    else:
        return SamplesQueryIndex(connection=connection)
