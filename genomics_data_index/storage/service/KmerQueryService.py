from pathlib import Path
from typing import List, Dict

import pandas as pd

from genomics_data_index.storage.index.KmerSearchManager import KmerSearchManagerSourmash
from genomics_data_index.storage.service.QueryService import QueryService
from genomics_data_index.storage.service.SampleService import SampleService


class KmerQueryService(QueryService):

    def __init__(self, sample_service: SampleService):
        super().__init__()
        self._sample_service = sample_service
        self._sourmash_search = KmerSearchManagerSourmash()

    def _find_matches_genome_files_internal(self, sample_reads: Dict[str, List[Path]],
                                            distance_threshold: float = None) -> pd.DataFrame:
        raise Exception('Not implemented')

    def _find_matches_internal(self, sample_names: List[str], distance_threshold: float):
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
            matches_df = matches_df.loc[:, matches_df['Distance'] <= distance_threshold]

        return matches_df.sort_values(['Query', 'Distance'], ascending=True)

    def _pairwise_distance_internal(self, samples: List[str]) -> pd.DataFrame:
        raise Exception('Not implemented')

    def _differences_between_genomes_internal(self, sample1: str, sample2: str):
        raise Exception('Not implemented')

    def get_data_type(self) -> str:
        return 'kmer'
