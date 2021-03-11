from pathlib import Path

from storage.variant.service import DatabaseConnection
from storage.variant.model import Sample
from storage.variant.model import SampleKmerIndex


class KmerService():

    def __init__(self, database: DatabaseConnection):
        self._database = database

    def insert_kmer_index(self, sample: Sample, kmer_index_path: Path):
        kmer_index = SampleKmerIndex(kmer_index_path=kmer_index_path)
        kmer_index.sample = sample
        self._database.get_session().add(kmer_index)
        self._database.get_session().commit()
