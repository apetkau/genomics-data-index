import logging
import tempfile
import time
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
from sourmash.fig import load_matrix_and_labels

from genomics_data_index.storage.util import execute_commands

logger = logging.getLogger(__name__)


class KmerSearchManagerSourmash:

    def __init__(self, ncores: int = 1):
        self._ncores = ncores

    def search(self, kmer_size: int, similarity_threshold: float, query_file: Path,
               search_files: List[Path]) -> pd.DataFrame:
        start_time = time.time()
        logger.debug(f'Start search for matches to [{query_file}].')

        with tempfile.TemporaryDirectory() as tmp:
            out_file = Path(tmp, 'output.csv')
            command = ['sourmash', 'search', '-k', str(kmer_size), '--threshold', str(similarity_threshold),
                       '-o', str(out_file),
                       str(query_file)]
            command.extend([str(f) for f in search_files])

            execute_commands([command])
            end_time = time.time()
            logger.debug(f'Finished search for matches. Took {end_time - start_time:0.2f} seconds')

            return pd.read_csv(out_file)

    def distances(self, kmer_size: int, signature_files: List[Path]) -> Tuple[np.ndarray, List[str]]:
        sigs_to_file_threshold = 1000

        start_time = time.time()
        logger.debug(f'Start creating pairwise distance matrix.')

        with tempfile.TemporaryDirectory() as tmp:
            if len(signature_files) > sigs_to_file_threshold:
                signatures_list_file = Path(tmp, 'signatures')
                with open(signatures_list_file, 'w') as fh:
                    for file in signature_files:
                        fh.write(f'{file}\n')

                files_to_search = ['--from-file', str(signatures_list_file)]
            else:
                files_to_search = [str(f) for f in signature_files]

            out_file = Path(tmp, 'output-matrix')
            command = ['sourmash', 'compare', '-k', str(kmer_size), '--output', str(out_file), '-p', str(self._ncores)]
            command.extend(files_to_search)

            execute_commands([command])
            end_time = time.time()
            logger.debug(f'Finished search for matches. Took {end_time - start_time:0.2f} seconds')

            similarity_matrix, labels = load_matrix_and_labels(str(out_file))
            distance_matrix = 1 - similarity_matrix

            return distance_matrix, labels
