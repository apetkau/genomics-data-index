import shutil
from pathlib import Path
from typing import List
import pandas as pd

from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.index.KmerSearchManager import KmerSearchManagerSourmash
from genomics_data_index.storage.model.db import Sample, SampleKmerIndex
from genomics_data_index.storage.service import DatabaseConnection
from genomics_data_index.storage.service.SampleService import SampleService


class KmerService:

    def __init__(self, database_connection: DatabaseConnection, features_dir: Path, sample_service: SampleService):
        self._database = database_connection
        self._sample_service = sample_service
        self._features_dir = features_dir
        self._sourmash_search = KmerSearchManagerSourmash()

    def find_matches_within(self, sample_names: List[str], distance_threshold: float) -> SampleSet:
        all_samples = self._sample_service.get_samples()
        kmer_index_paths = [s.sample_kmer_index.kmer_index_path for s in all_samples if s.sample_kmer_index is not None]
        kmer_size = 31

        matches_df = pd.DataFrame(data=[], columns=[
            'Query',
            'Match',
            'Similarity',
            'Distance',
        ])

        for sample_name in sample_names:
            query_sample = self._sample_service.get_sample(sample_name)
            results_df = self._sourmash_search.search(kmer_size, query_sample.sample_kmer_index.kmer_index_path,
                                                      kmer_index_paths)
            results_df['Distance'] = 1 - results_df['similarity']
            results_df['Query'] = sample_name
            results_df = results_df.rename({
                'name': 'Match',
                'similarity': 'Similarity',
            }, axis='columns')
            results_df = results_df[['Query', 'Match', 'Similarity', 'Distance']]

            matches_df = pd.concat([matches_df, results_df])

        if distance_threshold is not None:
            matches_df = matches_df[matches_df['Distance'] <= distance_threshold]

        sample_name_ids = self._sample_service.find_sample_name_ids(set(matches_df['Match'].tolist()))
        matches_set = SampleSet(sample_name_ids.values())

        return matches_set

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
