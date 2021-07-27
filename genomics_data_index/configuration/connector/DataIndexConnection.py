from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from genomics_data_index.configuration.connector.FilesystemStorage import FilesystemStorage
from genomics_data_index.storage.model.db.DatabasePathTranslator import DatabasePathTranslator
from genomics_data_index.storage.service import DatabaseConnection
from genomics_data_index.storage.service.CoreAlignmentService import CoreAlignmentService
from genomics_data_index.storage.service.KmerService import KmerService
from genomics_data_index.storage.service.MLSTService import MLSTService
from genomics_data_index.storage.service.ReferenceService import ReferenceService
from genomics_data_index.storage.service.SampleService import SampleService
from genomics_data_index.storage.service.TreeService import TreeService
from genomics_data_index.storage.service.VariationService import VariationService

logger = logging.getLogger(__name__)


class DataIndexConnection:

    def __init__(self, reference_service: ReferenceService, sample_service: SampleService,
                 variation_service: VariationService, alignment_service: CoreAlignmentService,
                 tree_service: TreeService, kmer_service: KmerService,
                 mlst_service: MLSTService, filesystem_storage: FilesystemStorage,
                 database_connection: DatabaseConnection):
        self._reference_service = reference_service
        self._sample_service = sample_service
        self._variation_service = variation_service
        self._alignment_service = alignment_service
        self._tree_service = tree_service
        self._kmer_service = kmer_service
        self._mlst_service = mlst_service
        self._filesystem_storage = filesystem_storage
        self._database_connection = database_connection

    @property
    def database(self):
        return self._database_connection

    @property
    def reference_service(self):
        return self._reference_service

    @property
    def sample_service(self):
        return self._sample_service

    @property
    def variation_service(self):
        return self._variation_service

    @property
    def alignment_service(self):
        return self._alignment_service

    @property
    def tree_service(self):
        return self._tree_service

    @property
    def kmer_service(self):
        return self._kmer_service

    @property
    def mlst_service(self):
        return self._mlst_service

    @property
    def filesystem_storage(self):
        return self._filesystem_storage

    def db_size(self, unit: str = 'B') -> pd.DataFrame:
        if unit == 'B':
            factor = 1
        elif unit == 'KB':
            factor = 1024
        elif unit == 'MB':
            factor = 1024 ** 2
        elif unit == 'GB':
            factor = 1024 ** 3
        else:
            raise Exception(f'Unknown unit=[{unit}]')

        filesystem_df = self.filesystem_storage.get_storage_size()
        database_df = self._database_connection.get_database_size()

        size_df = pd.concat([filesystem_df, database_df])
        total_data_size = size_df['Data Size'].sum()
        total_index_size = size_df['Index Size'].sum()
        total_items = size_df['Number of Items'].sum()
        total_row = pd.DataFrame([['Total', pd.NA, pd.NA, total_data_size, total_index_size, total_items]],
                                 columns=['Type', 'Name', 'Division', 'Data Size', 'Index Size', 'Number of Items'])
        size_df = pd.concat([size_df, total_row])

        # Reorder columns
        size_df = size_df[['Type', 'Name', 'Division', 'Data Size', 'Index Size', 'Number of Items']]

        size_df['Data Size'] = size_df['Data Size'] / factor
        size_df['Index Size'] = size_df['Index Size'] / factor
        size_df = size_df.rename({'Data Size': f'Data Size ({unit})',
                                  'Index Size': f'Index Size ({unit})'}, axis='columns')
        return size_df

    @classmethod
    def connect(cls, database_connection: str, database_dir: Path,
                index_unknown_missing: bool = True) -> DataIndexConnection:
        filesystem_storage = FilesystemStorage(Path(database_dir))
        dpt = DatabasePathTranslator(filesystem_storage.root_dir)
        logger.debug(f'Using database directory {database_dir}')

        logger.debug(f'Connecting to database {database_connection}')
        database = DatabaseConnection(connection_string=database_connection,
                                      database_path_translator=dpt)

        sql_select_limit = 500

        reference_service = ReferenceService(database, filesystem_storage.reference_dir)
        sample_service = SampleService(database)
        variation_service = VariationService(database_connection=database,
                                             variation_dir=filesystem_storage.variation_dir,
                                             reference_service=reference_service,
                                             sample_service=sample_service,
                                             index_unknown_missing=index_unknown_missing,
                                             sql_select_limit=sql_select_limit)

        alignment_service = CoreAlignmentService(database=database,
                                                 reference_service=reference_service,
                                                 sample_service=sample_service,
                                                 variation_service=variation_service)

        tree_service = TreeService(database, reference_service, alignment_service)

        kmer_service = KmerService(database_connection=database,
                                   sample_service=sample_service,
                                   features_dir=filesystem_storage.kmer_dir)

        mlst_service = MLSTService(database_connection=database, sample_service=sample_service,
                                   mlst_dir=filesystem_storage.mlst_dir)

        return DataIndexConnection(reference_service=reference_service, sample_service=sample_service,
                                   variation_service=variation_service, alignment_service=alignment_service,
                                   tree_service=tree_service, kmer_service=kmer_service, mlst_service=mlst_service,
                                   filesystem_storage=filesystem_storage, database_connection=database)
