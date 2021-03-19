from pathlib import Path

from storage.variant.service import DatabaseConnection
from storage.variant.service.SampleService import SampleService
from storage.variant.model import Sample
from storage.variant.model import SampleKmerIndex


class KmerService():

    def __init__(self, database_connection: DatabaseConnection, sample_service: SampleService):
        self._database = database_connection
        self._sample_service = sample_service

    def has_kmer_index(self, sample_name: str) -> bool:
        if self._sample_service.exists(sample_name):
            sample = self._sample_service.get_sample(sample_name)
            return sample.sample_kmer_index is not None
        else:
            return False

    def insert_kmer_index(self, sample_name: str, kmer_index_path: Path):
        kmer_index = SampleKmerIndex(kmer_index_path=kmer_index_path)

        if self._sample_service.exists(sample_name):
            sample = self._sample_service.get_sample(sample_name)
        else:
            sample = Sample(name=sample_name)
            self._database.get_session().add(sample)
        sample.sample_kmer_index = kmer_index
        self._database.get_session().commit()