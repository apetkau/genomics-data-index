from __future__ import annotations
from pathlib import Path

from storage.connector.DataIndexConnection import DataIndexConnection
from storage.api.SamplesQuery import SamplesQuery
from storage.api.impl.SamplesQueryIndex import SamplesQueryIndex


def connect(database_connection: str, database_dir: Path) -> DataIndexConnection:
    return DataIndexConnection.connect(database_connection=database_connection,
                                       database_dir=database_dir)


def query(connection: DataIndexConnection) -> SamplesQuery:
    return SamplesQueryIndex(connection=connection)
