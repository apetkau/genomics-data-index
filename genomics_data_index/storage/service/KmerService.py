import logging
import shutil
from pathlib import Path
from typing import List, Tuple, Union

import numpy as np
import pandas as pd

from genomics_data_index.storage.SampleSet import SampleSet
from genomics_data_index.storage.index.KmerSearchManager import KmerSearchManagerSourmash
from genomics_data_index.storage.model.db import Sample, SampleKmerIndex
from genomics_data_index.storage.service import DatabaseConnection
from genomics_data_index.storage.service.SampleService import SampleService

logger = logging.getLogger(__name__)


class KmerService:
    FIND_MATCHES_MERGE_TYPES = ['union']

    def __init__(self, database_connection: DatabaseConnection, features_dir: Path, sample_service: SampleService):
        self._database = database_connection
        self._sample_service = sample_service
        self._features_dir = features_dir
        self._sourmash_search = KmerSearchManagerSourmash()

    def find_matches_within(self, sample_names: List[str],
                            kmer_size: int, distance_threshold: float,
                            results_merge_type: str = 'union',
                            samples_universe: SampleSet = None) -> SampleSet:
        """
        Find samples within a particular distance of a list of samples. This is based on kmer signatures/sketches.

        :param sample_names: The list of sample names to search for matches.
        :param kmer_size: The kmer size to use for searching through signatures/sketches.
        :param distance_threshold: A number from 0 to 1, with 0 indicating the distance threshold to include matches to
                                   to the samples listed in 'sample_names'.
        :param results_merge_type: Defines how to combine results when passing multiple 'sample_names'. A type of
                                   'union' means that matches will be the union of all samples matching anything in
                                   'sample_names'. Currently only 'union' is supported (this parameter is here to
                                   make a bit more clear how results are combined until I implement additional ways of
                                   merging results).
        :param samples_universe: The universe of samples to search through. Can be used to restrict which samples we will
                           consider for matches. Set to 'None' to search for matches in all samples in the system.
        :return: A SampleSet of the matches.
        """
        if results_merge_type != 'union':
            raise Exception(f'results_merge_type=[{results_merge_type}] is not supported. '
                            f'Only {self.FIND_MATCHES_MERGE_TYPES} are supported.')

        if samples_universe is None:
            sample_universe_objects = self._sample_service.get_samples()
        else:
            sample_universe_objects = self._sample_service.find_samples_by_ids(sample_ids=samples_universe)

        kmer_index_paths = [s.sample_kmer_index.kmer_index_path for s in sample_universe_objects if
                            s.sample_kmer_index is not None]

        if len(kmer_index_paths) < len(sample_universe_objects):
            logger.debug(f'Not all samples (number={len(sample_universe_objects)} have associated kmer signatures '
                         f'(number={len(kmer_index_paths)}). These will be excluded from the search.')

        if len(kmer_index_paths) == 0:
            return SampleSet.create_empty()
        else:
            similarity_threshold = 1 - distance_threshold

            matches_df = pd.DataFrame(data=[], columns=[
                'Query',
                'Match',
                'Similarity',
                'Distance',
            ])

            for sample_name in sample_names:
                query_sample = self._sample_service.get_sample(sample_name)
                results_df = self._sourmash_search.search(kmer_size=kmer_size,
                                                          similarity_threshold=similarity_threshold,
                                                          query_file=query_sample.sample_kmer_index.kmer_index_path,
                                                          search_files=kmer_index_paths)
                results_df['Distance'] = 1 - results_df['similarity']
                results_df['Query'] = sample_name
                results_df = results_df.rename({
                    'name': 'Match',
                    'similarity': 'Similarity',
                }, axis='columns')
                results_df = results_df[['Query', 'Match', 'Similarity', 'Distance']]

                matches_df = pd.concat([matches_df, results_df])

            sample_name_ids = self._sample_service.find_sample_name_ids(set(matches_df['Match'].tolist()))
            matches_set = SampleSet(sample_name_ids.values())

            return matches_set

    def get_distance_matrix(self, sample_ids: Union[List[int], SampleSet], kmer_size: int,
                            ncores: int = 1) -> Tuple[
        np.ndarray, List[str]]:
        if isinstance(sample_ids, list):
            sample_ids = SampleSet(sample_ids)

        sourmash_search_multicore = KmerSearchManagerSourmash(ncores=ncores)

        samples = self._sample_service.find_samples_by_ids(sample_ids)
        kmer_index_paths = [s.sample_kmer_index.kmer_index_path for s in samples if
                            s.sample_kmer_index is not None]

        if len(kmer_index_paths) < len(samples):
            raise Exception(f'Not all samples (number={len(samples)}) have associated kmer signatures '
                            f'(number={len(kmer_index_paths)}).')

        return sourmash_search_multicore.distances(kmer_size=kmer_size, signature_files=kmer_index_paths)

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
