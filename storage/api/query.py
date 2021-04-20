from __future__ import annotations

from pathlib import Path

import pandas as pd

from storage.api.SamplesQuery import SamplesQuery
from storage.api.impl.TreeSamplesQuery import TreeSamplesQuery
from storage.api.impl.DataFrameSamplesQuery import DataFrameSamplesQuery
from storage.api.impl.SamplesQueryIndex import SamplesQueryIndex
from storage.connector.DataIndexConnection import DataIndexConnection

QUERY_UNIVERSE = ['all', 'mutations', 'dataframe']


def connect(database_connection: str, database_dir: Path) -> DataIndexConnection:
    return DataIndexConnection.connect(database_connection=database_connection,
                                       database_dir=database_dir)


def query(connection: DataIndexConnection, universe: str = 'all', **kwargs) -> SamplesQuery:
    if universe == 'all':
        return _query_all_samples(connection)
    elif universe == 'mutations':
        return _query_reference(connection=connection, **kwargs)
    elif universe == 'dataframe':
        return _query_data_frame(connection=connection, **kwargs)
    else:
        raise Exception(f'Invalid universe=[{universe}]. Must be one of {QUERY_UNIVERSE}')


def _query_all_samples(connection: DataIndexConnection):
    all_samples = connection.sample_service.get_all_sample_ids()
    return SamplesQueryIndex(connection=connection, sample_set=all_samples, universe_set=all_samples)


def _query_reference(connection: DataIndexConnection, reference_name: str):
    reference_samples = connection.sample_service.get_samples_associated_with_reference(reference_name)
    reference_genome = connection.reference_service.find_reference_genome(reference_name)

    sample_query = SamplesQueryIndex(connection=connection, sample_set=reference_samples,
                                     universe_set=reference_samples)

    if reference_genome.has_tree():
        sample_query = TreeSamplesQuery(connection=connection, wrapped_query=sample_query,
                                        tree=reference_genome.tree,
                                        alignment_length=reference_genome.tree_alignment_length)

    return sample_query


def _query_data_frame(connection: DataIndexConnection,
                      data_frame: pd.DataFrame = None,
                      sample_ids_column=None,
                      sample_names_column=None
                      ):

    if data_frame is None:
        raise Exception('data_frame must be set when querying with universe="dataframe"')
    if sample_ids_column is None and sample_names_column is None:
        raise Exception('If querying with universe="dataframe", then one of sample_names_column or sample_ids_column '
                        'must be set')
    elif sample_ids_column is not None:
        all_samples = connection.sample_service.get_all_sample_ids()
        return DataFrameSamplesQuery.create_with_sample_ids_column(sample_ids_column,
                                                                   data_frame=data_frame,
                                                                   database_sample_set=all_samples,
                                                                   connection=connection)
    else:
        all_samples = connection.sample_service.get_all_sample_ids()
        return DataFrameSamplesQuery.create_with_sample_names_column(sample_names_column,
                                                                     data_frame=data_frame,
                                                                     database_sample_set=all_samples,
                                                                     connection=connection)
