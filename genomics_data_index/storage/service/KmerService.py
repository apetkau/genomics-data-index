from pathlib import Path
import shutil

from genomics_data_index.storage.model.db import Sample, SampleKmerIndex
from genomics_data_index.storage.service import DatabaseConnection
from genomics_data_index.storage.service.SampleService import SampleService


class KmerService:

    def __init__(self, database_connection: DatabaseConnection, features_dir: Path, sample_service: SampleService):
        self._database = database_connection
        self._sample_service = sample_service
        self._features_dir = features_dir

    def has_kmer_index(self, sample_name: str) -> bool:
        if self._sample_service.exists(sample_name):
            sample = self._sample_service.get_sample(sample_name)
            return sample.sample_kmer_index is not None
        else:
            return False

    def insert_kmer_index(self, sample_name: str, kmer_index_path: Path):

        if self._sample_service.exists(sample_name):
            sample = self._sample_service.get_sample(sample_name)
        else:
            sample = Sample(name=sample_name)
            self._database.get_session().add(sample)

        kmer_path_internal = self._features_dir / kmer_index_path.name
        shutil.copy(kmer_index_path, kmer_path_internal)
        kmer_index = SampleKmerIndex(sample=sample, kmer_index_path=kmer_path_internal)
        self._database.get_session().add(kmer_index)
        self._database.get_session().commit()