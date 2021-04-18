from __future__ import annotations

import logging
from pathlib import Path

from storage.FilesystemStorage import FilesystemStorage
from storage.variant.service import DatabaseConnection
from storage.variant.service.CoreAlignmentService import CoreAlignmentService
from storage.variant.service.KmerQueryService import KmerQueryService
from storage.variant.service.KmerService import KmerService
from storage.variant.service.MLSTQueryService import MLSTQueryService
from storage.variant.service.MLSTService import MLSTService
from storage.variant.service.MutationQueryService import MutationQueryService
from storage.variant.service.ReferenceService import ReferenceService
from storage.variant.service.SampleService import SampleService
from storage.variant.service.TreeService import TreeService
from storage.variant.service.VariationService import VariationService

logger = logging.getLogger(__name__)


class DataIndexConnection:

    def __init__(self, reference_service: ReferenceService, sample_service: SampleService,
                 variation_service: VariationService, alignment_service: CoreAlignmentService,
                 tree_service: TreeService, mutation_query_service: MutationQueryService,
                 kmer_service: KmerService, kmer_query_service: KmerQueryService,
                 mlst_service: MLSTService, mlst_query_service: MLSTQueryService):
        self._reference_service = reference_service
        self._sample_service = sample_service
        self._variation_service = variation_service
        self._alignment_service = alignment_service
        self._tree_service = tree_service
        self._mutation_query_service = mutation_query_service
        self._kmer_service = kmer_service
        self._kmer_query_service = kmer_query_service
        self._mlst_service = mlst_service
        self._mlst_query_service = mlst_query_service

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
    def mutation_query_service(self):
        return self._mutation_query_service

    @property
    def kmer_service(self):
        return self._kmer_service

    @property
    def kmer_query_service(self):
        return self._kmer_query_service

    @property
    def mlst_service(self):
        return self._mlst_service

    @property
    def mlst_query_service(self):
        return self._mlst_query_service

    @classmethod
    def connect(cls, database_connection: str, database_dir: Path) -> DataIndexConnection:
        logger.debug(f'Connecting to database {database_connection}')
        database = DatabaseConnection(database_connection)
        filesystem_storage = FilesystemStorage(Path(database_dir))

        logger.info(f'Using database directory {database_dir}')
        reference_service = ReferenceService(database, filesystem_storage.reference_dir)
        sample_service = SampleService(database)
        variation_service = VariationService(database_connection=database,
                                             variation_dir=filesystem_storage.variation_dir,
                                             reference_service=reference_service,
                                             sample_service=sample_service)

        alignment_service = CoreAlignmentService(database=database,
                                                 reference_service=reference_service,
                                                 sample_service=sample_service,
                                                 variation_service=variation_service)

        tree_service = TreeService(database, reference_service, alignment_service)

        mutation_query_service = MutationQueryService(reference_service=reference_service,
                                                      sample_service=sample_service,
                                                      tree_service=tree_service)

        kmer_service = KmerService(database_connection=database,
                                   sample_service=sample_service)

        kmer_query_service = KmerQueryService(sample_service=sample_service)

        mlst_service = MLSTService(database_connection=database, sample_service=sample_service,
                                   mlst_dir=filesystem_storage.mlst_dir)
        mlst_query_service = MLSTQueryService(sample_service=sample_service,
                                              mlst_service=mlst_service)

        return DataIndexConnection(reference_service=reference_service, sample_service=sample_service,
                                   variation_service=variation_service, alignment_service=alignment_service,
                                   tree_service=tree_service, mutation_query_service=mutation_query_service,
                                   kmer_service=kmer_service, kmer_query_service=kmer_query_service,
                                   mlst_service=mlst_service, mlst_query_service=mlst_query_service)
